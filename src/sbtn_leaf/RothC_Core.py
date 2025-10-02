## File: rothc_core.py
"""
Core RothC model routines: pool initialization, decomposition (with CO2 tracking), optional spin-up, and simulation
"""

# -----------------------------------------------------------------------------
# MODULES
# -----------------------------------------------------------------------------
import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple, Optional
import polars as pl

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
    """Calculates the Tillage Rate Modifiers (TRM) of applying conservation tillage according to Hyun & Yoo (2024). This uses sand and SOC stock content from the soil to calculate TRM

    Args:
        sand (_type_): Sand Soil Content as percentage 
        SOC (_type_): SOC stock content of soil as t C/ha

    Returns:
        _type_: _description_
    """
    # Defines TRMs
    TRMs = {'DPM': [1.54, 1.71, 1.54, 0.72],
            'RPM': [0.35, 0.35, 2.15, 0.97],
            'BIO': [1.42, 0.38, 2.38, 0.99],
            'HUM': [0.42, 0.87, 2.93, 0.94]
            }
    
    # pull out the per-node lists as NumPy arrays
    dpm_vals = np.array(TRMs['DPM'], dtype=float)   # shape (4,)
    rpm_vals = np.array(TRMs['RPM'], dtype=float)
    bio_vals = np.array(TRMs['BIO'], dtype=float)
    hum_vals = np.array(TRMs['HUM'], dtype=float)

    # Load variables as array
    sand = np.asarray(sand, dtype=float)
    SOC = np.asarray(SOC, dtype=float)

    # Initializes TRMs
    TRM_DPM = np.zeros_like(SOC, dtype=float)
    TRM_RPM = np.zeros_like(SOC, dtype=float)
    TRM_BIO = np.zeros_like(SOC, dtype=float)
    TRM_HUM = np.zeros_like(SOC, dtype=float)

    # Initialize node classification
    nodes = np.zeros_like(SOC, dtype=int)

    # Classifies each cell into 4 different nodes
    SN1 = sand > 37.6   # (SN1)

    # Goes through SN 2
    high_soc = (SOC > 75.7) & (SN1)
    nodes = np.where(high_soc, 1, 2)    # TN1 if SOC > 75.7, TN2 if <=
    
    # Goes through SN 3
    high_sand = (sand > 35.0) & (~SN1)
    nodes = np.where(high_sand, 3, 4)

    # if sand > 37.6:     # (SN1)
    #    # Goes through Separation Node 2
    #    high_soc = SOC > 75.7
    #    nodes = np.where(high_soc, 1, 2)    # TN1 if SOC > 75.7, TN2 if <=
    # else:
    #    # Goes through SN3
    #    high_sand = sand > 35.0
    #    nodes = np.where(high_sand, 3, 4)

    # Transaling nodes into 0-index array
    nodes_idx = (nodes.astype(int) - 1)

    # Now updating TRMs
    TRM_DPM = dpm_vals[nodes_idx]   # shape (ny, nx)
    TRM_RPM = rpm_vals[nodes_idx]
    TRM_BIO = bio_vals[nodes_idx]
    TRM_HUM = hum_vals[nodes_idx]

    # Why does the above work:
    # When you index a 1-D array with another integer array of shape (ny, nx), NumPy builds a new array of shape (ny, nx) where each element is TRM_DPM[i,j] = dpm_vals[ node_idx[i,j] ]. 
    # In other words, at each pixel location (i,j) you look up which “node” index lives there, then grab the corresponding DPM value from dpm_vals.

    # Finally returning this
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
    dpm_rpm: float,
    trm: Optional[Tuple[float, float, float, float]] = None,
) -> CarbonPools:
    """
    Update carbon pools and track CO2 emitted this timestep.
    """
    # unpack current
    DPM_0, RPM_0, BIO_0, HUM_0, IOM_v, CO2_0 = (
        pools.DPM, pools.RPM, pools.BIO, pools.HUM, pools.IOM, pools.CO2
    )
    k = dict(DPM=10.0, RPM=0.3, BIO=0.66, HUM=0.02)
    if trm is None:
        trm_vals = dict(DPM=1.0, RPM=1.0, BIO=1.0, HUM=1.0)
    else:
        trm_vals = dict(DPM=trm[0], RPM=trm[1], BIO=trm[2], HUM=trm[3])
    tstep = 1.0 / time_fact

    # decay existing pools
    DPM_1 = DPM_0 * np.exp(-rate_m * trm_vals['DPM'] * k['DPM'] * tstep)
    RPM_1 = RPM_0 * np.exp(-rate_m * trm_vals['RPM'] * k['RPM'] * tstep)
    BIO_1 = BIO_0 * np.exp(-rate_m * trm_vals['BIO'] * k['BIO'] * tstep)
    HUM_1 = HUM_0 * np.exp(-rate_m * trm_vals['HUM'] * k['HUM'] * tstep)

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
    swc: float,
    trm: Optional[Tuple[float, float, float, float]] = None
) -> Tuple[CarbonPools, float]:
    """Advance pools one month"""
    rm_tmp = RMF_Tmp(temp)
    rm_moist, new_swc = RMF_Moist(rain, evap, clay, depth, cover, swc)
    rm_pc = RMF_PC(cover)
    rate_m = rm_tmp * rm_moist * rm_pc
    new_pools = decomp(12.0, pools, rate_m, clay, c_inp, fym, dpm_rpm, trm=trm)
    return new_pools, new_swc

# -----------------------------------------------------------------------------
# Shared monthly loop utilities
# -----------------------------------------------------------------------------

TRMProvider = Callable[[int, CarbonPools], Optional[Tuple[float, float, float, float]]]


def _run_monthly_sequence(
    pools: CarbonPools,
    clay: float,
    depth: float,
    tmp: np.ndarray,
    rain: np.ndarray,
    evap: np.ndarray,
    pc: np.ndarray,
    dpm_rpm: float,
    c_inp: Optional[np.ndarray],
    fym: Optional[np.ndarray],
    swc: float,
    *,
    irr: Optional[np.ndarray] = None,
    trm_provider: Optional[TRMProvider] = None,
) -> Tuple[CarbonPools, float, np.ndarray, np.ndarray]:
    """Advance the RothC state through a sequence of monthly timesteps.

    Parameters
    ----------
    pools : CarbonPools
        Current pool state that will be updated in-place.
    tmp, rain, evap, pc : array-like
        Monthly time-series driving the model. All arrays must share the same
        first dimension length.
    dpm_rpm : float
        Decomposable-to-resistant plant material ratio.
    c_inp, fym : array-like or None
        Monthly external carbon and FYM inputs. When ``None`` they are treated
        as zero inputs.
    swc : float
        Soil water content carried between months.
    irr : array-like or None, optional
        Additional monthly irrigation amounts. When provided they are added to
        ``rain`` before computing soil moisture.
    trm_provider : callable, optional
        Function returning per-pool tillage rate modifiers for the current
        month. It receives the month index (0-based) and the current
        ``CarbonPools`` instance and must return a ``tuple`` of four floats
        ``(DPM, RPM, BIO, HUM)`` or ``None`` when no modifier is applied.

    Returns
    -------
    CarbonPools
        Updated pool state after processing all months.
    float
        Updated soil water content to carry forward.
    np.ndarray
        Monthly SOC log for the processed months.
    np.ndarray
        Monthly CO2 emissions for the processed months.
    """

    tmp_arr = np.asarray(tmp, dtype=float)
    rain_arr = np.asarray(rain, dtype=float)
    evap_arr = np.asarray(evap, dtype=float)
    pc_arr = np.asarray(pc, dtype=int)
    months = tmp_arr.shape[0]

    c_arr = np.zeros(months, dtype=float) if c_inp is None else np.asarray(c_inp, dtype=float)
    f_arr = np.zeros(months, dtype=float) if fym is None else np.asarray(fym, dtype=float)
    irr_arr = None if irr is None else np.asarray(irr, dtype=float)

    soc_monthly = np.empty(months, dtype=float)
    co2_monthly = np.empty(months, dtype=float)
    swc_curr = float(swc)

    for month_idx in range(months):
        prev_co2 = pools.CO2
        water = rain_arr[month_idx]
        if irr_arr is not None:
            water = water + irr_arr[month_idx]
        trm_val = trm_provider(month_idx, pools) if trm_provider is not None else None
        pools, swc_curr = onemonth_step_rothc(
            pools,
            clay,
            depth,
            tmp_arr[month_idx],
            water,
            evap_arr[month_idx],
            int(pc_arr[month_idx]),
            dpm_rpm,
            c_arr[month_idx],
            f_arr[month_idx],
            swc_curr,
            trm=trm_val,
        )
        soc_monthly[month_idx] = pools.SOC
        co2_monthly[month_idx] = pools.CO2 - prev_co2

    return pools, swc_curr, soc_monthly, co2_monthly

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
    irr12: Optional[np.ndarray] = None,
    trm_provider: Optional[TRMProvider] = None,
    tol: float = 1e-6,
    max_cycles: int = 10000
) -> Tuple[
    CarbonPools, float,
    np.ndarray, np.ndarray,
    np.ndarray, np.ndarray
]:
    """Iterate yearly cycles until SOC change ``< tol``.

    Parameters
    ----------
    pools : CarbonPools
        Initial carbon pools, typically produced by :func:`initialize_pools`.
    clay, depth : float
        Soil properties used by the moisture modifier.
    tmp12, rain12, evap12, pc12 : array-like
        Monthly climate forcing for a single year.
    dpm_rpm : float
        Decomposable-to-resistant plant material ratio.
    c12, f12 : array-like or None, optional
        Monthly external carbon and FYM inputs. ``None`` defaults to zeros.
    irr12 : array-like or None, optional
        Optional monthly irrigation totals (same length as ``tmp12``). These
        are unique to reduced tillage scenarios and are added to precipitation
        before soil moisture is evaluated.
    trm_provider : callable, optional
        Hook that injects per-pool tillage rate modifiers. Supply a callable
        returning ``(DPM, RPM, BIO, HUM)`` multipliers to enable reduced
        tillage behaviour; pass ``None`` for conventional runs.
    tol : float, optional
        SOC convergence threshold between yearly cycles.
    max_cycles : int, optional
        Maximum number of annual cycles to evaluate.

    Returns
    -------
    tuple
        Updated pools, final soil water content, and monthly/annual SOC and
        CO2 logs for each cycle.
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
        pools, swc_cycle, monthly_soc, monthly_co2 = _run_monthly_sequence(
            pools,
            clay,
            depth,
            tmp12,
            rain12,
            evap12,
            pc12,
            dpm_rpm,
            c12,
            f12,
            swc,
            irr=irr12,
            trm_provider=trm_provider,
        )

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
    fym: Optional[np.ndarray] = None,
    irr: Optional[np.ndarray] = None,
    trm_provider: Optional[TRMProvider] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run RothC for ``n_years`` and return annual SOC and CO2 series.

    Parameters
    ----------
    clay, depth, iom, soc0 : float
        Site-specific soil characteristics.
    tmp, rain, evap, pc : array-like
        Monthly climate and cover time-series of length ``n_years * 12``.
    dpm_rpm : float
        Decomposable-to-resistant plant material ratio.
    n_years : int
        Simulation horizon in years.
    do_equilibrium : bool, optional
        When ``True`` perform an equilibrium spin-up using the first 12 months.
    c_inp, fym : array-like or None, optional
        Monthly plant carbon and manure inputs. Defaults to zeros when ``None``.
    irr : array-like or None, optional
        Optional monthly irrigation totals. Provide these, together with a
        ``trm_provider``, to model reduced tillage scenarios.
    trm_provider : callable, optional
        Callable returning tillage rate modifiers ``(DPM, RPM, BIO, HUM)`` for
        each month. Leave ``None`` for conventional management.

    Returns
    -------
    tuple of np.ndarray
        Annual SOC and CO2 totals (length ``n_years``).
    """
    pools = initialize_pools(soc0, clay)
    pools.IOM = iom
    swc = 0.0

    months = n_years * 12
    tmp = np.asarray(tmp, dtype=float)
    rain = np.asarray(rain, dtype=float)
    evap = np.asarray(evap, dtype=float)
    pc = np.asarray(pc, dtype=int)
    c_inp = np.zeros(months, dtype=float) if c_inp is None else np.asarray(c_inp, dtype=float)
    fym = np.zeros(months, dtype=float) if fym is None else np.asarray(fym, dtype=float)
    irr = None if irr is None else np.asarray(irr, dtype=float)

    # optional spin-up
    if do_equilibrium:
        _, swc, _, _, _, _ = run_equilibrium(
            pools, clay, depth,
            tmp[:12], rain[:12], evap[:12], pc[:12],
            dpm_rpm, c_inp[:12], fym[:12],
            irr12=None if irr is None else irr[:12],
            trm_provider=trm_provider,
        )

    # run full simulation
    pools, swc, soc_monthly, co2_monthly = _run_monthly_sequence(
        pools,
        clay,
        depth,
        tmp,
        rain,
        evap,
        pc,
        dpm_rpm,
        c_inp,
        fym,
        swc,
        irr=irr,
        trm_provider=trm_provider,
    )
    soc_annual = soc_monthly[11::12]
    co2_annual = co2_monthly.reshape(n_years, 12).sum(axis=1)

    return soc_annual, co2_annual
