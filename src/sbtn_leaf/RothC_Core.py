## File: rothc_core.py
"""
Core RothC model routines: pool initialization, decomposition (with CO2 tracking), optional spin-up, and simulation
"""

# -----------------------------------------------------------------------------
# MODULES
# -----------------------------------------------------------------------------
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional
import polars as pl

# -----------------------------------------------------------------------------
# Constants used by the Tillage Rate Modifier (TRM) calculation
# -----------------------------------------------------------------------------
_TRM_COEFFICIENTS = np.array(
    [
        [1.54, 1.71, 1.54, 0.72],  # DPM
        [0.35, 0.35, 2.15, 0.97],  # RPM
        [1.42, 0.38, 2.38, 0.99],  # BIO
        [0.42, 0.87, 2.93, 0.94],  # HUM
    ],
    dtype=float,
)

# -----------------------------------------------------------------------------
# Data structure for carbon pools (including CO2 emissions)
# -----------------------------------------------------------------------------
@dataclass
class CarbonPools:
    DPM: float  # Decomposable Plant Material (t C/ha)
    RPM: float  # Resistant Plant Material (t C/ha)
    BIO: float  # Microbial Biomass (t C/ha)
    HUM: float  # Humified Organic Matter (t C/ha)
    IOM: float  # Inert Organic Matter (t C/ha)
    SOC: float  # Total Soil Organic Carbon (t C/ha)
    CO2: float  # Cumulative CO2-C emissions (t C/ha)

# -----------------------------------------------------------------------------
# Initialize pools from total SOC and clay content (no initial CO2)
# -----------------------------------------------------------------------------
def initialize_pools(soc_initial: float, clay: float) -> CarbonPools:
    IOM_initial = 0.049 * soc_initial**1.139
    RPM_initial = (0.1847 * soc_initial + 0.1555) * (clay + 1.275)**(-0.1158)
    HUM_initial = (0.7148 * soc_initial) * (clay + 0.3421)**(0.0184)
    BIO_initial = (0.014  * soc_initial + 0.0075) * (clay + 8.8473)**(0.0567)
    DPM_initial = soc_initial - IOM_initial - RPM_initial - HUM_initial - BIO_initial
    
    return CarbonPools(
        DPM=DPM_initial,
        RPM=RPM_initial,
        BIO=BIO_initial,
        HUM=HUM_initial,
        IOM=IOM_initial,
        SOC=soc_initial,
        CO2=0.0
    )

# -----------------------------------------------------------------------------
# Rate-modifying factors
# -----------------------------------------------------------------------------
def RMF_Tmp(temp: float):
    """
    Temperature modifying factor, vectorized over an array of temp (°C).
    temp < -5 → 0.0, else 47.91/(exp(106.06/(temp+18.27))+1).
    """
    tmp = np.asarray(temp, dtype=float)
    rmf = 47.91 / (np.exp(106.06 / (tmp + 18.27)) + 1.0)
    # zero out below -5
    return np.where(tmp < -5.0, 0.0, rmf)

def RMF_Moist(rain, evap, clay, depth, pc, swc):
    """Moisture modifying factor and updated soil water content
    The maximum Soil Moisture Deficit is calculated as  -(20.0 + 1.3 (%clay) - 0.01 (%clay)**2) for 23 cm. If calculating for a different depth, divide by 23 and multiply for th actual depth in cm.

    evap must be given in open-pan evaporation, and the model later multiplies it by 0.75. 


    """
    rain = np.asarray(rain, dtype=float)
    evap = np.asarray(evap, dtype=float)
    clay = np.asarray(clay, dtype=float)
    pc   = np.asarray(pc, dtype=int)
    swc  = np.asarray(swc, dtype=float)

    # Compute SMD limits
    SMDMax = -(20 + 1.3*clay - 0.01*clay**2)
    SMDAdj = SMDMax * depth / 23.0
    SMD1bar = 0.444 * SMDAdj
    SMDBare = 0.556 * SMDAdj

    # water deficit
    df = rain - 0.75 * evap

    # update SWC differently for bare vs vegetated
    swc_bare = np.clip(swc + df, SMDAdj, 0.0)
    swc_veg  = np.clip(swc + df, SMDBare, swc)
    new_swc  = np.where(pc == 1, swc_bare, swc_veg)

    # RMF (b)
    above = new_swc > SMD1bar
    rmf = np.where(
        above,
        1.0,
        0.2 + (1.0 - 0.2) * ((SMDAdj - new_swc) / (SMDAdj - SMD1bar))
    )
    return rmf, new_swc


def RMF_PC(pc: int):
    """
    Plant‐cover modifying factor, vectorized over array of pc (0/1).
    1 if bare (pc==0), 0.6 otherwise.
    """
    arr = np.asarray(pc)
    return np.where(arr == 0, 1.0, 0.6)


def RMF_TRM(sand, SOC):
    """Calculate the tillage rate modifiers (TRM) for conservation tillage.

    The implementation follows Hyun & Yoo (2024) and mirrors the previous
    behaviour of this routine, but relies on pre-computed coefficient tables so
    repeated allocation/initialisation is avoided.
    """

    sand = np.asarray(sand, dtype=float)
    SOC = np.asarray(SOC, dtype=float)

    # Classify each cell into one of the four decision-tree "nodes".  This
    # retains the original step-by-step logic even though it results in only
    # nodes 3 and 4 being reachable; the goal of the refactor is parity.
    SN1 = sand > 37.6
    high_soc = (SOC > 75.7) & SN1
    nodes = np.where(high_soc, 1, 2)
    high_sand = (sand > 35.0) & (~SN1)
    nodes = np.where(high_sand, 3, 4)
    nodes_idx = nodes.astype(int) - 1

    # Look up the TRM coefficients for every pool simultaneously.
    trm_values = np.take(_TRM_COEFFICIENTS, nodes_idx, axis=1)
    TRM_DPM, TRM_RPM, TRM_BIO, TRM_HUM = trm_values

    return TRM_DPM, TRM_RPM, TRM_BIO, TRM_HUM


# -----------------------------------------------------------------------------
# Decomposition step with CO2 emission calculation
# -----------------------------------------------------------------------------
def decomp(
    time_fact: float,
    pools: CarbonPools,
    rate_m: float,
    clay: float,
    plant_c_inp: float,
    fym_inp: float,
    dpm_rpm: float
) -> CarbonPools:
    """
    Update carbon pools and track CO2 emitted this timestep.
    """
    # unpack current
    DPM_0, RPM_0, BIO_0, HUM_0, IOM_v, CO2_0 = (
        pools.DPM, pools.RPM, pools.BIO, pools.HUM, pools.IOM, pools.CO2
    )
    k = dict(DPM=10.0, RPM=0.3, BIO=0.66, HUM=0.02)
    tstep = 1.0 / time_fact

    # decay existing pools
    DPM_1 = DPM_0 * np.exp(-rate_m * k['DPM'] * tstep)
    RPM_1 = RPM_0 * np.exp(-rate_m * k['RPM'] * tstep)
    BIO_1 = BIO_0 * np.exp(-rate_m * k['BIO'] * tstep)
    HUM_1 = HUM_0 * np.exp(-rate_m * k['HUM'] * tstep)

    # losses in each pool
    losses = {
        'DPM': DPM_0 - DPM_1,
        'RPM': RPM_0 - RPM_1,
        'BIO': BIO_0 - BIO_1,
        'HUM': HUM_0 - HUM_1
    }

    # CO2/BIO+HUM partition ratio
    x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay))
    resp_frac = x / (x + 1)

    # compute CO2 emission
    total_co2 = (
        losses['DPM'] * resp_frac +
        losses['RPM'] * resp_frac +
        losses['BIO'] * resp_frac +
        losses['HUM'] * resp_frac
    )

    # remaining partition to BIO/HUM
    def BIO_HUM_part(val): return val * (0.46 / (x + 1)), val * (0.54 / (x + 1))
    DPM_BIO, DPM_HUM = BIO_HUM_part(losses['DPM'])
    RPM_BIO, RPM_HUM = BIO_HUM_part(losses['RPM'])
    BIO_BIO, BIO_HUM = BIO_HUM_part(losses['BIO'])
    HUM_BIO, HUM_HUM = BIO_HUM_part(losses['HUM'])

    # assemble new pools before adding external Carbon inputs
    new_DPM = DPM_1
    new_RPM = RPM_1
    new_BIO = BIO_1 + DPM_BIO + RPM_BIO + BIO_BIO + HUM_BIO
    new_HUM = HUM_1 + DPM_HUM + RPM_HUM + BIO_HUM + HUM_HUM

    # add external inputs
    # Divides plant carbon input between DPM and RPM
    pi_dpm = dpm_rpm / (dpm_rpm + 1.0) * plant_c_inp
    pi_rpm = plant_c_inp - pi_dpm
    # Divides farm manure into dpm, rpm, and hum
    f_dpm = 0.49 * fym_inp
    f_rpm = 0.49 * fym_inp
    f_hum = 0.02 * fym_inp

    new_DPM += pi_dpm + f_dpm
    new_RPM += pi_rpm + f_rpm
    new_HUM += f_hum

    new_SOC = new_DPM + new_RPM + new_BIO + new_HUM + IOM_v
    new_CO2 = CO2_0 + total_co2

    return CarbonPools(
        DPM=new_DPM,
        RPM=new_RPM,
        BIO=new_BIO,
        HUM=new_HUM,
        IOM=IOM_v,
        SOC=new_SOC,
        CO2=new_CO2
    )

# -----------------------------------------------------------------------------
# One-month RothC step
# -----------------------------------------------------------------------------
def onemonth_step_rothc(
    pools: CarbonPools,
    clay: float,
    depth: float,
    temp: float,
    rain: float,
    evap: float,
    cover: int,
    dpm_rpm: float,
    c_inp: float,
    fym: float,
    swc: float
) -> Tuple[CarbonPools, float]:
    """Advance pools one month"""
    rm_tmp = RMF_Tmp(temp)
    rm_moist, new_swc = RMF_Moist(rain, evap, clay, depth, cover, swc)
    rm_pc = RMF_PC(cover)
    rate_m = rm_tmp * rm_moist * rm_pc
    new_pools = decomp(12.0, pools, rate_m, clay, c_inp, fym, dpm_rpm)
    return new_pools, new_swc

# -----------------------------------------------------------------------------
# Equilibrium spin-up logging (formerly run_equilibrium)
# -----------------------------------------------------------------------------
def run_equilibrium(
    pools: CarbonPools,
    clay: float,
    depth: float,
    tmp12: np.ndarray,
    rain12: np.ndarray,
    evap12: np.ndarray,
    pc12: np.ndarray,
    dpm_rpm: float,
    c12: Optional[np.ndarray] = None,
    f12: Optional[np.ndarray] = None,
    tol: float = 1e-6,
    max_cycles: int = 10000
) -> Tuple[
    CarbonPools, float,
    np.ndarray, np.ndarray,
    np.ndarray, np.ndarray
]:
    """
    Iterate yearly cycles until SOC change < tol.

    Variables:
      - cycle: automatically increments via 'for cycle in range(max_cycles)'
      - prev_soc: initialized to current pools.SOC at start of each cycle

    Returns:
      - pools at equilibrium
      - final SWC
      - soc_monthly_log: array of shape (cycles,12)
      - co2_monthly_log: array of shape (cycles,12)
      - soc_annual_log: array of shape (cycles,)
      - co2_annual_log: array of shape (cycles,)
    """
    
    # default external inputs = 0
    if c12 is None:
        c12 = np.zeros_like(tmp12, dtype=float)
    if f12 is None:
        f12 = np.zeros_like(tmp12, dtype=float)

    swc = 0.0
    prev_soc = pools.SOC

    # logs: list of cycles × 12
    soc_monthly_log = []
    co2_monthly_log = []
    soc_annual_log = []
    co2_annual_log = []

    # Loop cycles; 'cycle' index auto-updates each iteration
    for cycle in range(max_cycles):
        # Store SOC at beginning of cycle for convergence check
        prev_soc = pools.SOC
        # Prepare containers for this cycle
        monthly_soc = np.empty(12, dtype=float)
        monthly_co2 = np.empty(12, dtype=float)
        swc_cycle = swc

        # Run 12-month loop
        for k in range(12):
            # Track CO2 before step
            prev_co2 = pools.CO2
            pools, swc_cycle = onemonth_step_rothc(
                pools, clay, depth,
                tmp12[k], rain12[k], evap12[k], int(pc12[k]),
                dpm_rpm, c12[k], f12[k], swc_cycle
            )
            # Log end-of-month SOC and CO2 emitted that month
            monthly_soc[k] = pools.SOC
            monthly_co2[k] = pools.CO2 - prev_co2

        # Annual metrics for this cycle
        annual_soc = monthly_soc[-1]          # December SOC
        annual_co2 = monthly_co2.sum()        # sum of monthly CO2

        # Append logs
        soc_monthly_log.append(monthly_soc)
        co2_monthly_log.append(monthly_co2)
        soc_annual_log.append(annual_soc)
        co2_annual_log.append(annual_co2)

        # Update SWC for next cycle
        swc = swc_cycle

        # Check convergence: SOC change between cycles
        if abs(pools.SOC - prev_soc) < tol:
            # Converged; break out of cycles
            break

    return (
        pools,
        swc,
        np.array(soc_monthly_log),
        np.array(co2_monthly_log),
        np.array(soc_annual_log),
        np.array(co2_annual_log)
    )

# -----------------------------------------------------------------------------
# Run full simulation, optional equilibrium, returning annual SOC and CO2
# -----------------------------------------------------------------------------
def run_simulation(
    clay: float,
    depth: float,
    iom: float,
    soc0: float,
    tmp: np.ndarray,
    rain: np.ndarray,
    evap: np.ndarray,
    pc: np.ndarray,
    dpm_rpm: float,
    n_years: int,
    do_equilibrium: bool = False,
    c_inp: Optional[np.ndarray] = None,
    fym: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Run RothC for n_years; return (soc_annual, co2_annual) arrays"""
    tmp = np.asarray(tmp, dtype=float)
    rain = np.asarray(rain, dtype=float)
    evap = np.asarray(evap, dtype=float)
    pc = np.asarray(pc, dtype=int)
    if c_inp is None:
        c_inp = np.zeros_like(tmp, dtype=float)
    else:
        c_inp = np.asarray(c_inp, dtype=float)
    if fym is None:
        fym = np.zeros_like(tmp, dtype=float)
    else:
        fym = np.asarray(fym, dtype=float)

    pools = initialize_pools(soc0, clay)
    pools.IOM = iom
    swc = 0.0

    # optional spin-up
    if do_equilibrium:
        _, swc, _, _, _, _ = run_equilibrium(
            pools, clay, depth,
            tmp[:12], rain[:12], evap[:12], pc[:12],
            dpm_rpm, c_inp[:12], fym[:12]
        )

    # run full simulation
    months = n_years * 12
    soc_monthly = np.empty(months, dtype=float)
    co2_monthly = np.empty(months, dtype=float)
    for t in range(months):
        prev_co2 = pools.CO2
        pools, swc = onemonth_step_rothc(
            pools, clay, depth,
            tmp[t], rain[t], evap[t], int(pc[t]),
            dpm_rpm, c_inp[t], fym[t], swc
        )
        soc_monthly[t] = pools.SOC
        co2_monthly[t] = pools.CO2 - prev_co2
    soc_annual = soc_monthly[11::12]
    co2_annual = co2_monthly.reshape(n_years, 12).sum(axis=1)
    
    return soc_annual, co2_annual
