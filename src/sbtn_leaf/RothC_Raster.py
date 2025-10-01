"""
Raster I/O and processing utilities for RothC modeling
"""

# -----------------------------------------------------------------------------
# MODULES
# -----------------------------------------------------------------------------
import xarray as xr
import rioxarray as rxr
import rasterio
from rasterio.enums import Resampling
import numpy as np
from typing import List, Optional, Tuple
import polars as pl
from tqdm import trange
from tqdm import tqdm

from sbtn_leaf.RothC_Core import RMF_Tmp, RMF_Moist, RMF_PC, RMF_TRM

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------
def load_single_band(path: str) -> xr.DataArray:
    """
    Load a single-band raster as an xarray DataArray with spatial coords.
    """
    da = rxr.open_rasterio(path, masked=True).squeeze()
    
    return da

def load_multiband(path: str) -> xr.DataArray:
    """
    Load a multi-band raster (e.g., 12 bands) as xarray DataArray with 'band' dimension.
    """
    da = rxr.open_rasterio(path, masked=True)
    da = da.rename({'band': 'time'})
    
    return da

def align_and_resample(
    reference: xr.DataArray,
    others: List[xr.DataArray],
    resampling: Resampling = Resampling.nearest
) -> List[xr.DataArray]:
    """
    Align and resample each DataArray in 'others' to match the grid of 'reference'.
    """
    aligned = []
    
    for da in others:
        da2 = da.rio.reproject_match(reference, resampling=resampling)  # reprojects other array into the reference one
        aligned.append(da2)
    
    return aligned


def mask_by_landuse(data: xr.DataArray, mask: xr.DataArray) -> xr.DataArray:
    """
    Apply a binary land-use mask (0/1) to a DataArray, setting values to NaN where mask==0.
    """
    return data.where(mask == 1)


def stack_time_series(baseline: xr.DataArray, years: int) -> xr.DataArray:
    """
    Repeat a 12-step baseline DataArray for multiple years along a new 'time' dimension.
    """
    arr = np.tile(baseline.values, (years, 1, 1))
    times = np.arange(len(arr))
    return xr.DataArray(arr, 
                        dims=('time', 'y', 'x'), 
                        coords={'time': times, 'y': baseline.y, 'x': baseline.x}
                        )

def build_pc_mask(
    cover_df: pl.DataFrame,
    template: xr.DataArray
) -> xr.DataArray:
    """
    Build a monthly plant-cover mask (0/1) using a Polars DataFrame of monthly cover values and a template grid.

    Parameters:
      - cover_df: Polars DataFrame with columns ['Month','Plant_Cover'] for 12 months.
      - template: xr.DataArray with dims ('time','y','x'), time length multiple of 12.

    Returns:
      - xr.DataArray mask with dims ('time','y','x') matching template coords.
    """
    # Extract sorted 12-month cover values from Polars DataFrame
    vals = np.array(cover_df.sort('Month')['Plant_Cover'].to_list())

    n_time = template.sizes.get('time')
    if n_time is None or n_time % 12 != 0:
        raise ValueError("Template 'time' dimension must exist and its length be a multiple of 12")
    years = n_time // 12

    # Tile 1D cover curve across all months
    cover_1d = np.tile(vals, years)

    # Build mask array and broadcast to spatial dims
    mask_arr = cover_1d.reshape(n_time, 1, 1)
    mask_arr = np.broadcast_to(
        mask_arr,
        (n_time, template.sizes['y'], template.sizes['x'])
    )

    return xr.DataArray(
        mask_arr,
        dims=('time', 'y', 'x'),
        coords={'time': template.time, 'y': template.y, 'x': template.x},
        name='plant_cover'
    )


def write_single_band_tif(
    data: xr.DataArray,
    out_path: str,
    template: xr.DataArray
) -> None:
    """
    Write a single-band DataArray (dims 'y','x') to a GeoTIFF using template metadata.
    """
    meta = template.rio.profile.copy()
    meta.update(count=1, dtype='float32')
    with rasterio.open(out_path, 'w', **meta) as dst:
        dst.write(
            data.values.astype('float32'),
            1
        )


def write_multiband_tif(
    data: xr.DataArray,
    out_path: str,
    template: xr.DataArray
) -> None:
    """
    Write a multi-band DataArray (time as band) to a GeoTIFF using template for metadata.
    """
    meta = template.rio.profile.copy()
    meta.update(count=data.sizes['time'], dtype='float32')
    with rasterio.open(out_path, 'w', **meta) as dst:
        for i in range(data.sizes['time']):
            dst.write(data.isel(time=i).astype('float32').values, i+1)


# -----------------------------------------------------------------------------
# Raster RothC simulation
# -----------------------------------------------------------------------------
def raster_rothc_annual_only(
    clay: np.ndarray,
    depth: float,
    soc0: np.ndarray,
    tmp: np.ndarray,
    rain: np.ndarray,
    evap: np.ndarray,
    pc: np.ndarray,
    c_inp: Optional[np.ndarray] = None,
    fym:   Optional[np.ndarray] = None,
    dpm_rpm: float = 1.44,
    n_years: int   = 100
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized RothC that returns only annual SOC and CO2.
    
    Parameters
    ----------
    clay, soc0 : 2D (y, x)
    tmp, rain, evap, pc, c_inp, fym : 3D (time, y, x)
      time can be 12 (annual cycle) or n_years*12.
    depth : float (cm)
    dpm_rpm : float
    n_years : int
    
    Returns
    -------
    soc_annual : ndarray (years, y, x)
    co2_annual : ndarray (years, y, x)
    """
    def _expand(arr: np.ndarray, months: int) -> np.ndarray:
        t, y, x = arr.shape
        if t == months:
            return arr
        if t == 12:
            return np.tile(arr, (n_years, 1, 1))
        raise ValueError(f"time dim {t}, expected 12 or {months}")
    
    months = n_years * 12
    tmp  = _expand(tmp,  months)
    rain = _expand(rain, months)
    evap = _expand(evap, months)
    pc   = _expand(pc,   months)
    c_inp = _expand(c_inp, months) if c_inp is not None else np.zeros_like(tmp)
    fym   = _expand(fym,   months) if fym   is not None else np.zeros_like(tmp)
    
    # Initialize pools
    with np.errstate(invalid='ignore'):
        IOM = 0.049 * soc0**1.139
        RPM = (0.1847*soc0 + 0.1555)*(clay + 1.275)**(-0.1158)
        HUM = (0.7148*soc0)      *(clay + 0.3421)**(0.0184)
        BIO = (0.014*soc0 + 0.0075)*(clay + 8.8473)**(0.0567)
    DPM = soc0 - (IOM + RPM + HUM + BIO)
    SOC = soc0.copy()
    swc = np.zeros_like(soc0)
    
    # Prepare annual outputs
    soc_annual = np.zeros((n_years,)+soc0.shape, dtype=np.float32)
    co2_annual = np.zeros_like(soc_annual)
    annual_co2_acc = np.zeros_like(soc0, dtype=np.float32)
    
    dt = 1.0 / 12.0
    for t in trange(months, desc="RothC months (annual only)"):
        # Rate-modifying factors
        rm_tmp = RMF_Tmp(tmp[t])
        rm_moist, swc = RMF_Moist(rain[t], evap[t], clay, depth, pc[t], swc)
        rm_pc = RMF_PC(pc[t])
        rate_m = rm_tmp * rm_moist * rm_pc
        
        # Decomposition
        D1 = DPM * np.exp(-rate_m * 10.0 * dt)
        R1 = RPM * np.exp(-rate_m *  0.3 * dt)
        B1 = BIO * np.exp(-rate_m *  0.66 * dt)
        H1 = HUM * np.exp(-rate_m *  0.02 * dt)
        
        lossD, lossR, lossB, lossH = DPM - D1, RPM - R1, BIO - B1, HUM - H1
        x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay))
        resp_frac = x / (x + 1.0)
        total_co2 = (lossD + lossR + lossB + lossH) * resp_frac
        
        # Pool partition
        def part(arr):
            return arr * (0.46/(x+1.0)), arr * (0.54/(x+1.0))
        D2B, D2H = part(lossD); R2B, R2H = part(lossR)
        B2B, B2H = part(lossB); H2B, H2H = part(lossH)
        
        # Update pools
        DPM = D1 + (dpm_rpm/(dpm_rpm+1.0))*c_inp[t] + 0.49*fym[t]
        RPM = R1 + (1.0/(dpm_rpm+1.0))*c_inp[t] + 0.49*fym[t]
        BIO = B1 + D2B + R2B + B2B + H2B
        HUM = H1 + D2H + R2H + B2H + H2H
        SOC = DPM + RPM + BIO + HUM + IOM
        
        # Accumulate CO2 this month
        annual_co2_acc += total_co2.astype(np.float32)
        
        # End-of-year: record and reset CO2 accumulator
        if (t + 1) % 12 == 0:
            yi = (t + 1)//12 - 1
            soc_annual[yi] = SOC.astype(np.float32)
            co2_annual[yi] = annual_co2_acc
            annual_co2_acc[:] = 0
    
    return soc_annual, co2_annual



def raster_rothc_annual_results_1yrloop(
    n_years: int,
    clay: np.ndarray,
    soc0: np.ndarray,
    tmp: np.ndarray,
    rain: np.ndarray,
    evap: np.ndarray,
    pc: np.ndarray,
    irr: Optional[np.ndarray]   = None,
    c_inp: Optional[np.ndarray] = None,
    fym:   Optional[np.ndarray] = None,
    depth: float = 23,
    dpm_rpm: float = 1.44,
    soc0_nodatavalue = -32768.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized RothC that returns only annual SOC and CO2.
    
    Parameters
    ----------
    clay, soc0 : 2D (y, x)
    tmp, rain, evap, pc, c_inp, fym : 3D (time, y, x)
      time can be 12 (annual cycle) or n_years*12.
    depth : float (cm)
    dpm_rpm : float
    n_years : int
    
    Returns
    -------
    soc_annual : ndarray (years, y, x)
    co2_annual : ndarray (years, y, x)
    """
    t_dim, y, x = tmp.shape
    months = n_years * 12

    # Initialize c_inp and fym if no input given
    c_inp = c_inp if c_inp is not None else np.zeros_like(tmp)
    fym   = fym   if fym   is not None else np.zeros_like(tmp)
    
    # Initialize pools
    with np.errstate(invalid='ignore'):
        IOM = 0.049 * soc0**1.139
        RPM = (0.1847*soc0 + 0.1555)*(clay + 1.275)**(-0.1158)
        HUM = (0.7148*soc0)      *(clay + 0.3421)**(0.0184)
        BIO = (0.014*soc0 + 0.0075)*(clay + 8.8473)**(0.0567)
    DPM = soc0 - (IOM + RPM + HUM + BIO)
    SOC = soc0.copy()
    swc = np.zeros_like(soc0)
    
    # Prepare annual outputs
    soc_annual = np.zeros((n_years + 1,y,x), dtype=np.float32)
    co2_annual = np.zeros_like(soc_annual)
    annual_co2_acc = np.zeros((y, x), dtype=np.float32)

    # Year 0 state
    arr0 = soc0.astype(np.float32)
    arr0[arr0 == soc0_nodatavalue] = np.nan
    soc_annual[0] = arr0
    co2_annual[0] = 0
    
    dt = 1.0 / 12.0
    for t_abs in trange(months, desc="RothC months", position=1):
        # Re-assigning t_abs into t
        t = t_abs % 12

        # Includes irrigation if provided
        if irr is not None:
            wat = rain[t] + irr[t]
        else:
            wat = rain[t]

        # Rate-modifying factors
        rm_tmp = RMF_Tmp(tmp[t])
        rm_moist, swc = RMF_Moist(wat, evap[t], clay, depth, pc[t], swc)
        rm_pc = RMF_PC(pc[t])
        rate_m = rm_tmp * rm_moist * rm_pc
        
        # Decomposition
        D1 = DPM * np.exp(-rate_m * 10.0 * dt)
        R1 = RPM * np.exp(-rate_m *  0.3 * dt)
        B1 = BIO * np.exp(-rate_m *  0.66 * dt)
        H1 = HUM * np.exp(-rate_m *  0.02 * dt)
        
        lossD, lossR, lossB, lossH = DPM - D1, RPM - R1, BIO - B1, HUM - H1
        x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay))
        resp_frac = x / (x + 1.0)
        total_co2 = (lossD + lossR + lossB + lossH) * resp_frac
        
        # Pool partition
        def part(arr):
            return arr * (0.46/(x+1.0)), arr * (0.54/(x+1.0))
        D2B, D2H = part(lossD); R2B, R2H = part(lossR)
        B2B, B2H = part(lossB); H2B, H2H = part(lossH)
        
        # Update pools
        DPM = D1 + (dpm_rpm/(dpm_rpm+1.0))*c_inp[t] + 0.49*fym[t]
        RPM = R1 + (1.0/(dpm_rpm+1.0))*c_inp[t] + 0.49*fym[t]
        BIO = B1 + D2B + R2B + B2B + H2B
        HUM = H1 + D2H + R2H + B2H + H2H
        SOC = DPM + RPM + BIO + HUM + IOM
        
        # Accumulate CO2 this month
        annual_co2_acc += total_co2.astype(np.float32)
        
        # End-of-year: record and reset CO2 accumulator
        if (t_abs + 1) % 12 == 0:
            yi = (t_abs + 1)//12
            soc_annual[yi] = SOC.astype(np.float32)
            co2_annual[yi] = annual_co2_acc
            annual_co2_acc[:] = 0
    
    return soc_annual, co2_annual


def raster_rothc_ReducedTillage_annual_results_1yrloop(
    n_years: int,
    clay: np.ndarray,
    soc0: np.ndarray,
    tmp: np.ndarray,
    rain: np.ndarray,
    evap: np.ndarray,
    pc: np.ndarray,
    sand: np.ndarray,
    irr: Optional[np.ndarray]   = None,
    c_inp: Optional[np.ndarray] = None,
    fym:   Optional[np.ndarray] = None,
    depth: float = 23,
    dpm_rpm: float = 1.44,
    soc0_nodatavalue: float = -32768
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized RothC for reduced tillage that returns only annual SOC and CO2.
    
    Parameters
    ----------
    clay, soc0 : 2D (y, x)
    tmp, rain, evap, pc, c_inp, fym : 3D (time, y, x)
      time can be 12 (annual cycle) or n_years*12.
    depth : float (cm)
    dpm_rpm : float
    n_years : int
    
    Returns
    -------
    soc_annual : ndarray (years, y, x)
    co2_annual : ndarray (years, y, x)
    """
    t_dim, y, x = tmp.shape
    months = n_years * 12

    # Initialize c_inp and fym if no input given
    c_inp = c_inp if c_inp is not None else np.zeros_like(tmp)
    fym   = fym   if fym   is not None else np.zeros_like(tmp)
    
    # Initialize pools
    with np.errstate(invalid='ignore'):
        IOM = 0.049 * soc0**1.139
        RPM = (0.1847*soc0 + 0.1555)*(clay + 1.275)**(-0.1158)
        HUM = (0.7148*soc0)      *(clay + 0.3421)**(0.0184)
        BIO = (0.014*soc0 + 0.0075)*(clay + 8.8473)**(0.0567)
    DPM = soc0 - (IOM + RPM + HUM + BIO)
    SOC = soc0.copy()
    swc = np.zeros_like(soc0)
    
    # Prepare annual outputs
    soc_annual = np.zeros((n_years + 1,y,x), dtype=np.float32)
    co2_annual = np.zeros_like(soc_annual)
    annual_co2_acc = np.zeros((y, x), dtype=np.float32)

    # Year 0 state
    arr0 = soc0.astype(np.float32)
    arr0[arr0 == soc0_nodatavalue] = np.nan
    soc_annual[0] = arr0
    co2_annual[0] = 0
    
    dt = 1.0 / 12.0
    for t_abs in trange(months, desc="RothC months", position=1):
        # Re-assigning t_abs into t
        t = t_abs % 12 

        # Includes irrigation if provided
        if irr is not None:
            wat = rain[t] + irr[t]
        else:
            wat = rain[t]

        # Rate-modifying factors
        rm_tmp = RMF_Tmp(tmp[t])
        rm_moist, swc = RMF_Moist(wat, evap[t], clay, depth, pc[t], swc)
        rm_pc = RMF_PC(pc[t])
        rate_m = rm_tmp * rm_moist * rm_pc

        # Tillage Rate Modifiers (TRMs)
        TRM_DPM, TRM_RPM, TRM_BIO, TRM_HUM = RMF_TRM(sand[t], SOC[t])
        
        # Decomposition
        D1 = DPM * np.exp(-rate_m * TRM_DPM * 10.0 * dt)
        R1 = RPM * np.exp(-rate_m * TRM_RPM * 0.3  * dt)
        B1 = BIO * np.exp(-rate_m * TRM_BIO * 0.66 * dt)
        H1 = HUM * np.exp(-rate_m * TRM_HUM * 0.02 * dt)
        
        lossD, lossR, lossB, lossH = DPM - D1, RPM - R1, BIO - B1, HUM - H1
        x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay))
        resp_frac = x / (x + 1.0)
        total_co2 = (lossD + lossR + lossB + lossH) * resp_frac
        
        # Pool partition
        def part(arr):
            return arr * (0.46/(x+1.0)), arr * (0.54/(x+1.0))
        D2B, D2H = part(lossD); R2B, R2H = part(lossR)
        B2B, B2H = part(lossB); H2B, H2H = part(lossH)
        
        # Update pools
        DPM = D1 + (dpm_rpm/(dpm_rpm+1.0))*c_inp[t] + 0.49*fym[t]
        RPM = R1 + (1.0/(dpm_rpm+1.0))*c_inp[t] + 0.49*fym[t]
        BIO = B1 + D2B + R2B + B2B + H2B
        HUM = H1 + D2H + R2H + B2H + H2H
        SOC = DPM + RPM + BIO + HUM + IOM
        
        # Accumulate CO2 this month
        annual_co2_acc += total_co2.astype(np.float32)
        
        # End-of-year: record and reset CO2 accumulator
        if (t_abs + 1) % 12 == 0:
            yi = (t_abs + 1)//12
            soc_annual[yi] = SOC.astype(np.float32)
            co2_annual[yi] = annual_co2_acc
            annual_co2_acc[:] = 0
    
    return soc_annual, co2_annual


# -----------------------------------------------------------------------------
# Other Related RothC Functions
# -----------------------------------------------------------------------------

# Function to save annual results
def save_annual_results(results_array, reference_raster, n_years, var_name, save_path, data_description: str, units: str='t C/ha', long_name: str = "Soil Organic Carbon", model_description: str = "RothC rasterized vectorized"):
    
    years = np.arange(1, n_years+1+1) # To include year 0
    
    # Construct the results array
    data_array = xr.DataArray(results_array, 
                              dims=('year','y','x'), 
                              coords={'year': years, 
                                      'y': reference_raster.y, 
                                      'x': reference_raster.x},
                              name=var_name)
    
    # Saving CRS and Transform
    data_array = (
        data_array
        .rio.write_crs(reference_raster.rio.crs)
        .rio.write_transform(reference_raster.rio.transform())
)

    # Define which dims correspond to spatial dims (only needed if not auto-detected)
    # da = da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)

    # 3. Write nodata value and additional tags
    data_array = data_array.rio.write_nodata(np.nan, inplace=False)

    # Write data description
    if data_description is None:
        data_description = f"RothC model results for {results_array} after {n_years}"

    # Update text metadata
    data_array.attrs.update({
        "units": units,
        "long_name": long_name,
        "model": model_description,
        "description": data_description 
    })

    # 4. Export to GeoTIFF with tags        
    data_array.rio.to_raster(save_path)


# Function to prepare all data
def load_environmental_data(lu_rp: str):
    # Loads data
    tmp = rxr.open_rasterio("../data/soil_weather/uhth_monthly_avg_temp_celsius.tif", masked=True)  # in °C
    rain = rxr.open_rasterio("../data/soil_weather/uhth_monthly_avg_precip.tif", masked=True)
    clay = rxr.open_rasterio("../data/soil_weather/uhth_clay_15-30cm_mean_perc.tif", masked=False).squeeze()
    soc0 = rxr.open_rasterio("../data/soil_weather/uhth_soc_0-30cm_mean.tif", masked=False).squeeze()
    sand = rxr.open_rasterio("../data/soil_weather/uhth_sand_15-30cm_mean_perc.tif", masked=False).squeeze()
    lu_raster = rxr.open_rasterio(lu_rp, masked=False).squeeze()

    # Creates IOM
    iom = 0.049 * soc0**1.139
    iom.attrs["units"]       = "t C/ha"
    iom.attrs["description"] = "IOM derived from SOC_initial"

    # Rename bands 
    tmp   = tmp[0].rename({'band': 'time'}) if isinstance(tmp, list) else tmp.rename({'band': 'time'})
    rain  = rain[0].rename({'band': 'time'}) if isinstance(rain, list) else rain.rename({'band': 'time'})

    # Mask data to land use requirementes
    lu_maks = (lu_raster==1)

    # Single‐band rasters
    clay   = clay.where(lu_maks)
    soc0   = soc0.where(lu_maks)
    iom    = iom.where(lu_maks)
    sand   = sand.where(lu_maks)

    # Multiband rasters (‘time’ × y × x)
    tmp    = tmp.where(lu_maks)
    rain   = rain.where(lu_maks)

    return tmp, rain, soc0, iom, clay, sand

def load_crop_data(lu_fp: str, evap_fp: str,  pc_fp: str, irr_fp: Optional[str], pr_fp: Optional[str], fym_fp: Optional[str]):
    # Opens land use data
    lu_raster = rxr.open_rasterio(lu_fp, masked=False).squeeze()
    lu_maks = (lu_raster==1)
        
    
    # Opens evap and pc, and process it
    evap = rxr.open_rasterio(evap_fp, masked=True)  # (12-band: Jan–Dec)
    evap  = evap.rename({"band": "time"})
    evap = evap.where(lu_maks).fillna(0)

    pc = rxr.open_rasterio(pc_fp, masked=True)
    pc    = pc.rename({"band": "time"})
    pc = (pc).where(lu_maks)
    pc = pc.where(lu_maks).fillna(0)
    
    # Optional inputs
    if irr_fp:
        irr = rxr.open_rasterio(irr_fp, masked=True)  # (12-band: Jan–Dec)
        irr = irr.rename({'band': 'time'})
        irr = (irr).where(lu_maks)
        irr = irr.where(lu_maks).fillna(0)
    else:
        irr = (xr.zeros_like(pc)).where(lu_maks)

    if pr_fp:
        pr = rxr.open_rasterio(pr_fp, masked=True)  # (12-band: Jan–Dec)
        pr = pr.rename({'band': 'time'})
        pr = (pr).where(lu_maks)
        pr = pr.where(lu_maks).fillna(0)
    else:
        pr    = (xr.zeros_like(pc)).where(lu_maks)
    
    if fym_fp:
        fym = rxr.open_rasterio(fym_fp, masked=True) # No farmyard manure in this case
        fym   = fym.rename({'band': 'time'})
        fym = (fym).where(lu_maks)
        fym = fym.where(lu_maks).fillna(0)
    else:
        fym    = (xr.zeros_like(pc)).where(lu_maks)


    return lu_raster, evap, pc, irr, pr, fym

def run_RothC(crop_name: str, practices_string_id: str, n_years: int, save_folder: str, data_description: str, lu_fp: str, evap_fp: str,  pc_fp: str, irr_fp: Optional[str] = None, pr_fp: Optional[str] = None, fym_fp: Optional[str] = None, red_till = False, save_CO2 = False):
    # Loads environmental data:
    print("Loading environmental data...")
    tmp, rain, soc0, iom, clay, sand = load_environmental_data(lu_fp)

    # Prepares crop data
    print("Loading crop data...")
    lu_raster, evap, pc, irr, pr, fym = load_crop_data(lu_fp, evap_fp,  pc_fp, irr_fp, pr_fp, fym_fp)
    
    # Convert to values
    clay_a, soc0_a, iom_a, sand_a = np.asarray(clay.values), np.asarray(soc0.values), np.asarray(iom.values), np.asarray(sand.values)
    tmp_a, rain_a, evap_a = np.asarray(tmp.values), np.asarray(rain.values), np.asarray(evap.values)
    pc_a, c_a, fym_a, irr_a = np.asarray(pc.values), np.asarray(pr.values), np.asarray(fym.values), np.asarray(irr.values)

    # Run model
    print("Running RothC...")
    if irr is not None:
        if red_till:
            SOC_results, CO2_results = raster_rothc_ReducedTillage_annual_results_1yrloop(
                n_years = n_years,
                clay    = clay_a,
                soc0    = soc0_a,
                tmp     = tmp_a,
                rain    = rain_a,
                evap    = evap_a,  
                pc      = pc_a,
                irr     = irr_a,
                c_inp   = c_a,
                fym     = fym_a,
                sand    = sand_a
            )
        else:
            SOC_results, CO2_results = raster_rothc_annual_results_1yrloop(
                n_years = n_years,
                clay    = clay_a,
                soc0    = soc0_a,
                tmp     = tmp_a,
                rain    = rain_a,
                evap    = evap_a,  
                pc      = pc_a,
                irr     = irr_a,
                c_inp   = c_a,
                fym     = fym_a
            )
    else:
        if red_till:
            SOC_results, CO2_results = raster_rothc_ReducedTillage_annual_results_1yrloop(
                n_years = n_years,
                clay    = clay_a,
                soc0    = soc0_a,
                tmp     = tmp_a,
                rain    = rain_a,
                evap    = evap_a,  
                pc      = pc_a,
                c_inp   = c_a,
                fym     = fym_a,
                sand    = sand_a
            )
        else:
            SOC_results, CO2_results = raster_rothc_annual_results_1yrloop(
                n_years = n_years,
                clay    = clay_a,
                soc0    = soc0_a,
                tmp     = tmp_a,
                rain    = rain_a,
                evap    = evap_a,  
                pc      = pc_a,
                c_inp   = c_a,
                fym     = fym_a
            )

    # Saving results
    string_save = f"{crop_name}_{practices_string_id}_{n_years}y_SOC.tif"
    save_path =f"{save_folder}/{string_save}"
    
    # SOC
    save_annual_results(SOC_results, lu_raster, n_years, "SOC", save_path, data_description, 't C/ha', long_name = "Soil Organic Carbon", model_description = "RothC rasterized vectorized")

    if save_CO2:
        save_annual_results(CO2_results, lu_raster, n_years, "CO2", save_path, data_description, 't CO2/ha', long_name = "CO2", model_description = "RothC rasterized vectorized")

    return SOC_results


def run_rothC_sceneraios_from_csv(csv_filepath):
    # 1) Read & cast your CSV exactly as before
    scenarios = (
        pl.read_csv(csv_filepath, null_values=["", "None"])
        .with_columns([
            pl.col("n_years").cast(pl.Int64),
            pl.col("red_till").cast(pl.Boolean),
            pl.col("save_CO2").cast(pl.Boolean),
        ])
    )

    # 2) Turn into a list of dicts once (so we know the total count)
    scenario_list = scenarios.to_dicts()

    # 3) Iterate with tqdm
    for scenario in scenario_list:
    #for scenario in tqdm(scenario_list, desc="Running RothC Scenarios", unit="scenario", position=0):
        # tqdm.write(f"Processing scenario: {scenario['practices_string_id']}")
        print(f"Running {scenario["crop_name"]} - {scenario["practices_string_id"]}")
        run_RothC(**scenario)
        print(f"\n\n")

##########################################
#### OTHER USEFUL FUNCTIONS FOR ROTHC ####
##########################################

def calcuate_annual_perc_changes(raster_path):
    # Open the raster
    da = rxr.open_rasterio(raster_path, masked=True)
    if isinstance(da, list):
        da = da[0]

    # Baseline at year 0
    baseline = da.isel(band = 0)

    # Calculate percentage changes
    pct = (da/baseline - 1) * 100

    # Eliminates infinites and replace them with NaNs
    pct = pct.where(np.isfinite(pct))
    
    return pct

def calcuate_practice_change_benefit(raster1_fp, raster2_fp, band_r1, band_r2):
    # Open the raster
    da1 = rxr.open_rasterio(raster1_fp, masked=True)
    da2 = rxr.open_rasterio(raster2_fp, masked=True)

    # Soc at year of band #
    da1_data = da1.isel(band = band_r1)
    da2_data = da2.isel(band = band_r2)

    # Calculate percentage changes
    pct = (da2_data/da1_data - 1) * 100

    # Eliminates infinites and replace them with NaNs
    pct = pct.where(np.isfinite(pct))
    
    return pct