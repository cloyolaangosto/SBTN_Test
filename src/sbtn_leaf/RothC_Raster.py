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
from typing import Any, Callable, Dict, List, Optional, Tuple, Union
import polars as pl
from tqdm.auto import trange, tqdm
from pathlib import Path
import json
import os

from sbtn_leaf.RothC_Core import RMF_Tmp, RMF_Moist, RMF_PC, RMF_TRM, _partition_to_bio_hum
import sbtn_leaf.cropcalcs as cropcalcs
from sbtn_leaf.paths import data_path

PathLike = Union[str, Path]


def _as_path(value: PathLike) -> Path:
    """Return ``value`` as a :class:`~pathlib.Path` instance."""

    return value if isinstance(value, Path) else Path(value)


def _resolve_optional_path(
    candidate: Optional[PathLike], *default: PathLike
) -> Path:
    """Return ``candidate`` coerced to :class:`Path` or fall back to ``data_path``.

    Parameters
    ----------
    candidate:
        Optional path override supplied by the caller.
    default:
        Path segments pointing to the default dataset beneath the repository
        ``data`` directory.  When ``candidate`` is :data:`None`,
        :func:`sbtn_leaf.paths.data_path` is invoked with these segments to
        produce the resolved location.
    """

    if candidate is None:
        return data_path(*default)
    return _as_path(candidate)


def _resolve_data_path(path: PathLike) -> Path:
    """Resolve a data path, falling back to :func:`data_path` when needed."""

    candidate = _as_path(path)
    if candidate.exists():
        return candidate

    data_candidate = data_path(candidate)
    if data_candidate.exists():
        return data_candidate

    return candidate

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------
def load_single_band(path: PathLike) -> xr.DataArray:
    """
    Load a single-band raster as an xarray DataArray with spatial coords.
    """
    da = rxr.open_rasterio(_as_path(path), masked=True).squeeze()
    
    return da

def load_multiband(path: PathLike) -> xr.DataArray:
    """
    Load a multi-band raster (e.g., 12 bands) as xarray DataArray with 'band' dimension.
    """
    da = rxr.open_rasterio(_as_path(path), masked=True)
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
    out_path: PathLike,
    template: xr.DataArray
) -> None:
    """
    Write a single-band DataArray (dims 'y','x') to a GeoTIFF using template metadata.
    """
    meta = template.rio.profile.copy()
    meta.update(count=1, dtype='float32')
    with rasterio.open(_as_path(out_path), 'w', **meta) as dst:
        dst.write(
            data.values.astype('float32'),
            1
        )


def write_multiband_tif(
    data: xr.DataArray,
    out_path: PathLike,
    template: xr.DataArray
) -> None:
    """
    Write a multi-band DataArray (time as band) to a GeoTIFF using template for metadata.
    """
    meta = template.rio.profile.copy()
    meta.update(count=data.sizes['time'], dtype='float32')
    with rasterio.open(_as_path(out_path), 'w', **meta) as dst:
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



TRMHandler = Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]


def _raster_rothc_annual_results(
    *,
    n_years: int,
    clay: np.ndarray,
    soc0: np.ndarray,
    tmp: np.ndarray,
    rain: np.ndarray,
    evap: np.ndarray,
    pc: np.ndarray,
    irr: Optional[np.ndarray],
    c_inp: Optional[np.ndarray],
    fym: Optional[np.ndarray],
    depth: float,
    commodity_type: str,
    soc0_nodatavalue: float,
    crop_name: Optional[str] = None,
    spam_crop_raster: Optional[str] = None,
    practices_string_id: Optional[str] = None,
    irr_yield_scaling: Optional[str] = None,
    spam_all_fp: Optional[str] = None,
    spam_irr_fp: Optional[str] = None,
    spam_rf_fp: Optional[str] = None,
    commodity_lu_fp: Optional[PathLike] = None,    # For grassland, crops residues calculations
    sand: Optional[np.ndarray] = None,
    forest_age:  Optional[np.ndarray] = None,
    forest_type: Optional[str] = None,
    grassland_type: Optional[str] = None,
    residue_runs: int = 100,
    weather_type: Optional[str] = None,
    TP_IPCC_bool: bool = False,
    trm_handler: Optional[TRMHandler],
    progress_desc: str = "RothC months",
    progress_position: Optional[int] = None,
    spam_outlier_strategy = str | None,
    spam_outlier_percentile = Tuple[float, float] | None,
    spam_outlier_k = float | None
) -> Tuple[np.ndarray, np.ndarray]:
    """Shared implementation for baseline and reduced tillage raster RothC runs."""

    months = n_years * 12
    _, y, x = tmp.shape
    crop_type: Optional[str] = None

    # Assigns dpm_rpm according to commodity type
    # checks if commodity type is correct
    if commodity_type not in ["annual_crop", "permanent_crop", "forest", "grassland"]:
        raise ValueError("Commodity type not valid. Valid types: annual_crop, permanent_crop, forest, grassland")
    
    c_inp_provided = c_inp is not None

    # Assigns dpm_rpm factor and other values if needed
    if commodity_type in ["annual_crop", "grassland"]:
        dpm_rpm = 1.44

        # Checks if grassland type is valid
        if grassland_type is not None and grassland_type not in ["natural", "managed"]:
            raise ValueError("Grassland type not valid. Valid types: natural, managed")

        if commodity_type == "annual_crop":
            crop_type = "annual"

            # initialize c_inp
            if practices_string_id is not None and "roff" in practices_string_id:
                print(f"        Residues removed. C_inputs are 0's")
                c_inp = np.zeros_like(rain)
            else:
                print(f"        Calculating baseline residue inputs for {crop_type} {crop_name} using {residue_runs} stochastic runs")

                c_inp = cropcalcs.calculate_monthly_residues_array(
                    lu_fp=commodity_lu_fp,
                    crop_name=crop_name,
                    crop_type=crop_type,
                    spam_crop_raster = spam_crop_raster,
                    irr_yield_scaling = irr_yield_scaling,
                    spam_all_fp = spam_all_fp,
                    spam_irr_fp = spam_irr_fp,
                    spam_rf_fp = spam_rf_fp,
                    random_runs=residue_runs,
                    print_outputs= False,
                    spam_outlier_strategy = spam_outlier_strategy,
                    spam_outlier_percentile = spam_outlier_percentile,
                    spam_outlier_k = spam_outlier_k
                )
            if c_inp is None:
                if practices_string_id is not None and "roff" in practices_string_id:
                    print(f"        C_inputs are 0's")
                    c_inp = np.zeros_like(rain)
                else:
                    print(f"        Calculating baseline residue inputs for {crop_type} {crop_name}")
                    c_inp = cropcalcs.calculate_monthly_residues_array(
                        lu_fp=commodity_lu_fp,
                        crop_name=crop_name,
                        crop_type=crop_type,
                        spam_crop_raster = spam_crop_raster,
                        irr_yield_scaling = irr_yield_scaling,
                        spam_all_fp = spam_all_fp,
                        spam_irr_fp = spam_irr_fp,
                        spam_rf_fp = spam_rf_fp,
                        random_runs=residue_runs,
                        spam_outlier_strategy = spam_outlier_strategy,
                        spam_outlier_percentile = spam_outlier_percentile,
                        spam_outlier_k = spam_outlier_k
                    )
            c_inp = np.squeeze(np.asarray(c_inp))
        
    elif commodity_type == "permanent_crop":
        dpm_rpm = 1
        crop_type = "permanent"

        # initialize c_inp
        if c_inp is None:
            print(f"        Calculating baseline residue inputs for {crop_type} {crop_name}")
            c_inp = cropcalcs.calculate_monthly_residues_array(
                lu_fp=commodity_lu_fp,
                crop_name=crop_name,
                crop_type=crop_type,
                spam_crop_raster=spam_crop_raster,
                irr_yield_scaling=irr_yield_scaling,
                spam_all_fp=spam_all_fp,
                spam_irr_fp=spam_irr_fp,
                spam_rf_fp=spam_rf_fp,
                random_runs=residue_runs,
                spam_outlier_strategy = spam_outlier_strategy,
                spam_outlier_percentile = spam_outlier_percentile,
                spam_outlier_k = spam_outlier_k
            )
        c_inp = np.squeeze(np.asarray(c_inp))

    else: # forest type
        dpm_rpm = 0.25
        
        # Checks that all forest inputs are there
        if forest_age is None or forest_type is None or weather_type is None or TP_IPCC_bool is None:
            raise ValueError("Missing forest inputs. Specify forest_age, forest_type, weather_type and TP_IPCC_bool")

        # initialize c_inp
        c_inp = cropcalcs.get_forest_litter_monthlyrate_fromda(
            forest_age,
            forest_type,
            weather_type,
            TP_IPCC_bool,
            residue_runs=residue_runs
        )

    # Initialize c_inp and fym if no input given
    if commodity_type != "forest":
        c_inp = c_inp if c_inp is not None else np.zeros_like(tmp)
    if c_inp is not None:
        c_inp = np.asarray(c_inp)
        if c_inp.ndim > 3:
            c_inp = np.squeeze(c_inp)
        if c_inp.ndim == 2:
            c_inp = np.broadcast_to(c_inp, tmp.shape)
        if c_inp.ndim != 3:
            raise ValueError("C input must be 3-D after squeezing/broadcasting")
    
    fym = fym if fym is not None else np.zeros_like(tmp)
    fym = np.asarray(fym)
    if fym.ndim > 3:
        fym = np.squeeze(fym)
    if fym.ndim not in (2, 3):
        raise ValueError("FYM input must be 2-D or 3-D after squeezing")

    # Initialize pools
    with np.errstate(invalid="ignore"):
        IOM = 0.049 * soc0 ** 1.139
        RPM = (0.1847 * soc0 + 0.1555) * (clay + 1.275) ** (-0.1158)
        HUM = (0.7148 * soc0) * (clay + 0.3421) ** (0.0184)
        BIO = (0.014 * soc0 + 0.0075) * (clay + 8.8473) ** (0.0567)
    DPM = soc0 - (IOM + RPM + HUM + BIO)
    SOC = soc0.copy()
    swc = np.zeros_like(soc0)

    # Prepare annual outputs
    soc_annual = np.zeros((n_years + 1, y, x), dtype=np.float32)
    co2_annual = np.zeros_like(soc_annual)
    annual_co2_acc = np.zeros((y, x), dtype=np.float32)

    # Year 0 state
    arr0 = soc0.astype(np.float32)
    arr0[arr0 == soc0_nodatavalue] = np.nan
    soc_annual[0] = arr0
    co2_annual[0] = 0

    dt = 1.0 / 12.0
    sand_has_time_dim = sand is not None and sand.ndim == 3

    position = 0 if progress_position is None else progress_position

    for t_abs in trange(months, desc=progress_desc, position=position):
        t = t_abs % 12
        year = t_abs // 12

        wat = rain[t] + irr[t] if irr is not None else rain[t]

        # Rate-modifying factors
        rm_tmp = RMF_Tmp(tmp[t])
        rm_moist, swc = RMF_Moist(wat, evap[t], clay, depth, pc[t], swc)
        rm_pc = RMF_PC(pc[t])
        rate_m = rm_tmp * rm_moist * rm_pc

        #TRM when sand is given for reduced tillage
        trm_dpm = trm_rpm = trm_bio = trm_hum = 1.0
        if trm_handler is not None:  # Only applicable for reduced tillage option
            if sand is None:
                raise ValueError("sand must be provided when a TRM handler is supplied")
            sand_current = sand[t] if sand_has_time_dim else sand
            trm_dpm, trm_rpm, trm_bio, trm_hum = trm_handler(sand_current, SOC)

        # Decomposition
        D1 = DPM * np.exp(-rate_m * trm_dpm * 10.0 * dt)
        R1 = RPM * np.exp(-rate_m * trm_rpm * 0.3 * dt)
        B1 = BIO * np.exp(-rate_m * trm_bio * 0.66 * dt)
        H1 = HUM * np.exp(-rate_m * trm_hum * 0.02 * dt)

        lossD, lossR, lossB, lossH = DPM - D1, RPM - R1, BIO - B1, HUM - H1
        x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay))
        resp_frac = x / (x + 1.0)
        total_co2 = (lossD + lossR + lossB + lossH) * resp_frac

        # Pool partition
        D2B, D2H = _partition_to_bio_hum(x, lossD)
        R2B, R2H = _partition_to_bio_hum(x, lossR)
        B2B, B2H = _partition_to_bio_hum(x, lossB)
        H2B, H2H = _partition_to_bio_hum(x, lossH)

        # Calculates carbon residue input if it's forest or grassland
        if commodity_type == "forest":
            c_inp_month = cropcalcs.get_forest_litter_monthlyrate_fromda(
                forest_age,
                forest_type,
                weather_type,
                TP_IPCC_bool,
                year_offset=year,
                residue_runs=residue_runs
            )
            c_inp_month = np.squeeze(np.asarray(c_inp_month))
            if c_inp_month.ndim != 2:
                raise ValueError("Forest litter input must be 2-D after squeezing")
        elif commodity_type == "grassland":
            c_annual = cropcalcs.generate_grassland_residue_map(grass_lu_fp=commodity_lu_fp, random_runs=residue_runs)  # Returns raster for 1 year
            pr_monthly = c_annual/12
            c_inp_month = np.squeeze(np.asarray(pr_monthly))
        else:
            c_inp_month = c_inp[t]

        # Harmonize fym
        if commodity_type == "grassland":
            fym = fym.squeeze()  # Forces fym to have only spatial dimensions
        fym_slice = fym if fym.ndim == 2 else fym[t]

        # Update pools
        DPM = D1 + (dpm_rpm / (dpm_rpm + 1.0)) * c_inp_month + 0.49 * fym_slice
        RPM = R1 + (1.0 / (dpm_rpm + 1.0)) * c_inp_month + 0.49 * fym_slice
        BIO = B1 + D2B + R2B + B2B + H2B
        HUM = H1 + D2H + R2H + B2H + H2H
        SOC = DPM + RPM + BIO + HUM + IOM

        annual_co2_acc += total_co2.astype(np.float32)

        if (t_abs + 1) % 12 == 0:
            yi = (t_abs + 1) // 12
            soc_annual[yi] = SOC.astype(np.float32)
            co2_annual[yi] = annual_co2_acc
            annual_co2_acc[:] = 0

            # Updates plant residue inputs for crops, skips final iteration:
            if (
                (t_abs + 1) < months
                and commodity_type in ("permanent_crop", "annual_crop")
                and not c_inp_provided
            ):
                # initialize c_inp
                if practices_string_id is not None and "roff" in practices_string_id:
                    print(f"        C_inputs are 0's")
                    c_inp = np.zeros_like(rain)
                else:
                    print(f"        Calculating baseline residue inputs for {crop_type} {crop_name}")
                    c_inp = cropcalcs.calculate_monthly_residues_array(
                        lu_fp=commodity_lu_fp,
                        crop_name=crop_name,
                        crop_type=crop_type,
                        spam_crop_raster = spam_crop_raster,
                        irr_yield_scaling = irr_yield_scaling,
                        spam_all_fp = spam_all_fp,
                        spam_irr_fp = spam_irr_fp,
                        spam_rf_fp = spam_rf_fp,
                        random_runs=residue_runs,
                        print_outputs= True,
                        spam_outlier_strategy = spam_outlier_strategy,
                        spam_outlier_percentile = spam_outlier_percentile,
                        spam_outlier_k = spam_outlier_k
                    )
                c_inp = np.squeeze(np.asarray(c_inp))

    return soc_annual, co2_annual


def raster_rothc_annual_results(
    n_years: int,
    clay: np.ndarray,
    soc0: np.ndarray,
    tmp: np.ndarray,
    rain: np.ndarray,
    evap: np.ndarray,
    pc: np.ndarray,
    commodity_type: str,
    irr: Optional[np.ndarray] = None,
    c_inp: Optional[np.ndarray] = None,
    fym: Optional[np.ndarray] = None,
    crop_name: Optional[str] = None,
    spam_crop_raster: Optional[str] = None,
    practices_string_id: Optional[str] = None,
    irr_yield_scaling: Optional[str] = None,
    spam_all_fp: Optional[str] = None,
    spam_irr_fp: Optional[str] = None,
    spam_rf_fp: Optional[str] = None,
    forest_age: Optional[np.ndarray] = None,
    forest_type: Optional[str] = None,
    commodity_lu_fp: Optional[PathLike] = None,
    grassland_type: Optional[str] = None,
    residue_runs: int = 100,
    weather_type: Optional[str] = None,
    TP_IPCC_bool: bool = False,
    depth: float = 15,
    soc0_nodatavalue: float = -32768.0,
    red_till: bool = False,
    sand: Optional[np.ndarray] = None,
    spam_outlier_strategy = str | None,
    spam_outlier_percentile = Tuple[float, float] | None,
    spam_outlier_k = float | None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized RothC that returns annual SOC and CO2.

    Parameters
    ----------
    clay, soc0 : 2D (y, x)
    tmp, rain, evap, pc, c_inp, fym : 3D (time, y, x)
      time can be 12 (annual cycle) or n_years*12.
    red_till : bool
        When ``True``, apply reduced-tillage rate modifiers using ``sand``.
    depth : float (cm)
    dpm_rpm : float
    n_years : int

    Returns
    -------
    soc_annual : ndarray (years, y, x)
    co2_annual : ndarray (years, y, x)
    """

    if red_till and sand is None:
        raise ValueError("sand must be provided when red_till is True")

    trm_handler = RMF_TRM if red_till else None

    return _raster_rothc_annual_results(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        irr=irr,
        c_inp=c_inp,
        fym=fym,
        sand=sand,
        depth=depth,
        commodity_type=commodity_type,
        soc0_nodatavalue=soc0_nodatavalue,
        trm_handler=trm_handler,
        forest_type=forest_type,
        commodity_lu_fp=commodity_lu_fp,
        grassland_type=grassland_type,
        residue_runs=residue_runs,
        weather_type=weather_type,
        TP_IPCC_bool=TP_IPCC_bool,
        forest_age=forest_age,
        crop_name=crop_name,
        spam_crop_raster=spam_crop_raster,
        practices_string_id=practices_string_id,
        irr_yield_scaling=irr_yield_scaling,
        spam_all_fp=spam_all_fp,
        spam_irr_fp=spam_irr_fp,
        spam_rf_fp=spam_rf_fp,
        spam_outlier_strategy = spam_outlier_strategy,
        spam_outlier_percentile = spam_outlier_percentile,
        spam_outlier_k = spam_outlier_k
    )


# -----------------------------------------------------------------------------
# Other Related RothC Functions
# -----------------------------------------------------------------------------

# Function to save annual results
def save_annual_results(
    results_array,
    reference_raster,
    n_years,
    var_name,
    save_path: PathLike,
    data_description: str,
    units: str = 't C/ha',
    long_name: str = "Soil Organic Carbon",
    model_description: str = "RothC rasterized vectorized",
):
    
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
    data_array.rio.to_raster(_as_path(save_path))


# Function to prepare all data
def _load_environmental_data(
    lu_rp: PathLike,
    *,
    tmp_fp: Optional[PathLike] = None,
    rain_fp: Optional[PathLike] = None,
    clay_fp: Optional[PathLike] = None,
    soc0_fp: Optional[PathLike] = None,
    sand_fp: Optional[PathLike] = None,
):
    # Loads data
    tmp = rxr.open_rasterio(
        _resolve_optional_path(
            tmp_fp,
            "soil_weather",
            "uhth_monthly_avg_temp_celsius.tif",
        ),
        masked=True,
    )  # in °C
    rain = rxr.open_rasterio(
        _resolve_optional_path(
            rain_fp,
            "soil_weather",
            "uhth_monthly_avg_precip.tif",
        ),
        masked=True,
    )
    clay = rxr.open_rasterio(
        _resolve_optional_path(
            clay_fp,
            "soil_weather",
            "uhth_clay_15-30cm_mean_perc.tif",
        ),
        masked=False,
    ).squeeze()
    soc0 = rxr.open_rasterio(
        _resolve_optional_path(
            soc0_fp,
            "soil_weather",
            "uhth_soc_0-30cm_mean.tif",
        ),
        masked=False,
    ).squeeze()
    sand = rxr.open_rasterio(
        _resolve_optional_path(
            sand_fp,
            "soil_weather",
            "uhth_sand_15-30cm_mean_perc.tif",
        ),
        masked=False,
    ).squeeze()
    lu_raster = rxr.open_rasterio(_as_path(lu_rp), masked=False).squeeze()

    # Creates IOM
    iom = 0.049 * soc0**1.139
    iom.attrs["units"]       = "t C/ha"
    iom.attrs["description"] = "IOM derived from SOC_initial"

    # Rename bands 
    tmp   = tmp[0].rename({'band': 'time'}) if isinstance(tmp, list) else tmp.rename({'band': 'time'})
    rain  = rain[0].rename({'band': 'time'}) if isinstance(rain, list) else rain.rename({'band': 'time'})

    # Mask data to land use requirementes
    lu_mask = (lu_raster==1)

    # Single‐band rasters
    clay   = clay.where(lu_mask)
    soc0   = soc0.where(lu_mask)
    iom    = iom.where(lu_mask)
    sand   = sand.where(lu_mask)

    # Multiband rasters (‘time’ × y × x)
    tmp    = tmp.where(lu_mask)
    rain   = rain.where(lu_mask)

    return tmp, rain, soc0, iom, clay, sand

def _load_crop_data(
    lu_fp: PathLike,
    evap_fp: PathLike,
    pc_fp: PathLike,
    irr_fp: Optional[PathLike],
    pr_fp: Optional[PathLike],
    fym_fp: Optional[PathLike],
):
    # Opens land use data
    lu_raster = rxr.open_rasterio(_as_path(lu_fp), masked=False).squeeze()
    lu_mask = (lu_raster==1)
        
    
    # Opens evap and pc, and process it
    evap = rxr.open_rasterio(_as_path(evap_fp), masked=True)  # (12-band: Jan–Dec)
    evap  = evap.rename({"band": "time"})
    evap = evap.where(lu_mask).fillna(0)

    pc = rxr.open_rasterio(_as_path(pc_fp), masked=True)
    pc    = pc.rename({"band": "time"})
    pc = (pc).where(lu_mask)
    pc = pc.where(lu_mask).fillna(0)
    
    # Optional inputs
    if irr_fp:
        irr = rxr.open_rasterio(_as_path(irr_fp), masked=True)  # (12-band: Jan–Dec)
        irr = irr.rename({'band': 'time'})
        irr = (irr).where(lu_mask)
        irr = irr.where(lu_mask).fillna(0)
    else:
        irr = (xr.zeros_like(pc)).where(lu_mask)

    if pr_fp:
        pr = rxr.open_rasterio(_as_path(pr_fp), masked=True)  # (12-band: Jan–Dec)
        pr = pr.rename({'band': 'time'})
        pr = (pr).where(lu_mask)
        pr = pr.where(lu_mask).fillna(0)
    else:
        pr    = (xr.zeros_like(pc)).where(lu_mask)
    
    if fym_fp:
        fym = rxr.open_rasterio(_as_path(fym_fp), masked=True) # No farmyard manure in this case
        fym   = fym.rename({'band': 'time'})
        fym = (fym).where(lu_mask)
        fym = fym.where(lu_mask).fillna(0)
    else:
        fym    = (xr.zeros_like(pc)).where(lu_mask)


    return lu_raster, evap, pc, irr, pr, fym


def _load_forest_data(lu_fp: PathLike, evap_fp: PathLike, age_fp: PathLike):
    # 1) Land-use (mask where class == 1)
    lu_raster = rxr.open_rasterio(_as_path(lu_fp), masked=True).squeeze()   # (y, x)
    lu_mask = (lu_raster == 1)  
    
    # Opens evap and and process it
    evap = rxr.open_rasterio(_as_path(evap_fp), masked=True)  # (12-band: Jan–Dec)
    evap  = evap.rename({"band": "time"})
    evap = evap.where(lu_mask).fillna(0)

    # Creates 12-month plant cover: 1 where lu_maks is True, 0 elsewhere 
    base = lu_mask.astype('float32')
    pc = xr.concat([base] * 12, dim='time')
    pc = pc.assign_coords(time=np.arange(1, 13))
    pc.name = "plant_cover"
    # Write spatial metadata so pc lines up with lu/evap
    pc = pc.rio.write_crs(lu_raster.rio.crs)
    pc = pc.rio.write_transform(lu_raster.rio.transform())

    # open age
    age = rxr.open_rasterio(_as_path(age_fp), masked=True).squeeze()
    age = age.where(lu_mask)

    return lu_raster, evap, pc, age
    
def _load_grassland_data(
    lu_fp: PathLike,
    evap_fp: PathLike,
    grassland_type: str,
    fym_fp: list[PathLike],
    pc_fp: Optional[PathLike] = None,
    irr_fp: Optional[PathLike] = None,
    pr_fp: Optional[PathLike] = None,
    residue_runs: int = 100,
):
    # Opens land use data
    lu_raster = rxr.open_rasterio(_as_path(lu_fp), masked=False).squeeze()
    lu_mask = (lu_raster==1)

    # Opens evap and pc, and process it
    evap = rxr.open_rasterio(_as_path(evap_fp), masked=True)  # (12-band: Jan–Dec)
    evap  = evap.rename({"band": "time"})
    evap = evap.where(lu_mask).fillna(0)

    # Plant residues
    if pr_fp is None:  # Calculates plant_residues based on a different raster if given
        pr_fp = lu_fp
    pr_annual = cropcalcs.generate_grassland_residue_map(
        grass_lu_fp=_as_path(pr_fp),
        random_runs=residue_runs,
    )  # Returns raster for 1 year
    pr_monthly = pr_annual/12
    pr = np.where(lu_mask, pr_monthly, 0)
    
    # Plant Cover
    if grassland_type == "natural":
        # Creates 12-month plant cover: 1 where lu_maks is True, 0 elsewhere 
        base = lu_mask.astype('int16')
        pc = xr.concat([base] * 12, dim='time')
        pc = pc.assign_coords(time=np.arange(1, 13))
        pc.name = "plant_cover"
        # Write spatial metadata so pc lines up with lu/evap
        pc = pc.rio.write_crs(lu_raster.rio.crs)
        pc = pc.rio.write_transform(lu_raster.rio.transform())
    else:   # TODO: managed grassland
        raise ValueError("Managed grasslands not yet implemented")
        # pc = rxr.open_rasterio(_as_path(pc_fp), masked=True)
        # pc = pc.rename({"band": "time"})
        # pc = (pc).where(lu_mask)
        # pc = pc.where(lu_mask).fillna(0)
    
    # Opens and returns each fym_fp for each animal
    fym_all = np.zeros_like(lu_mask)
    for fp in fym_fp:
        fym_annual = rxr.open_rasterio(_as_path(fp), masked=True) # No farmyard manure in this case
        fym_monthly = fym_annual/12
        fym_all = fym_all + fym_monthly
        fym = (fym_all).where(lu_mask)
        fym = fym.where(lu_mask).fillna(0)

    # Optional irrigation input - POTENTIALLY USED ON MANAGED GRASSLANDS
    if irr_fp:
        irr = rxr.open_rasterio(_as_path(irr_fp), masked=True)  # (12-band: Jan–Dec)
        irr = irr.rename({'band': 'time'})
        irr = (irr).where(lu_mask)
        irr = irr.where(lu_mask).fillna(0)
    else:
        irr = (xr.zeros_like(pc)).where(lu_mask)

    return lu_raster, evap, pr, pc, fym, irr


def _run_rothc_scenario(
    *,
    lu_fp: PathLike,
    n_years: int,
    save_folder: PathLike,
    data_description: str,
    result_basename: str,
    loader: Callable[..., Tuple[xr.DataArray, Dict[str, Any]]],
    loader_kwargs: Optional[Dict[str, Any]] = None,
    runner: Callable[..., Any],
    runner_kwargs: Optional[Dict[str, Any]] = None,
    loader_message: Optional[str] = None,
    save_CO2: bool = False,
    env_overrides: Optional[Dict[str, PathLike]] = None,
):
    """Shared workflow for RothC scenario execution and persistence."""

    loader_kwargs = loader_kwargs or {}
    runner_kwargs = runner_kwargs or {}

    print("    Loading environmental data...")
    tmp, rain, soc0, iom, clay, sand = _load_environmental_data(
        lu_fp, **(env_overrides or {})
    )

    env_arrays = {
        "tmp": np.asarray(tmp.values),
        "rain": np.asarray(rain.values),
        "soc0": np.asarray(soc0.values),
        "iom": np.asarray(iom.values),
        "clay": np.asarray(clay.values),
        "sand": np.asarray(sand.values),
    }

    if loader_message:
        print(loader_message)
    lu_raster, scenario_inputs = loader(lu_fp=_as_path(lu_fp), **loader_kwargs)

    print("    Running RothC...")

    results = runner(
        n_years=n_years,
        env=env_arrays,
        scenario=scenario_inputs,
        **runner_kwargs,
    )

    if isinstance(results, tuple):
        SOC_results, CO2_results = results
    else:
        SOC_results, CO2_results = results, None

    save_path = _as_path(save_folder) / result_basename

    save_annual_results(
        SOC_results,
        lu_raster,
        n_years,
        "SOC",
        save_path,
        data_description,
        't C/ha',
        long_name="Soil Organic Carbon",
        model_description="RothC rasterized vectorized",
    )

    print(f"    RothC results saved into {save_path}")

    if save_CO2 and CO2_results is not None:
        save_annual_results(
            CO2_results,
            lu_raster,
            n_years,
            "CO2",
            save_path,
            data_description,
            't CO2/ha',
            long_name="CO2",
            model_description="RothC rasterized vectorized",
        )

    return SOC_results


def run_RothC_crops(
    crop_name: str,
    commodity_type: str,
    n_years: int,
    save_folder: PathLike,
    data_description: str,
    lu_fp: PathLike,
    evap_fp: PathLike,
    pc_fp: PathLike,
    practices_string_id: Optional[str] = None,
    irr_fp: Optional[PathLike] = None,
    pr_fp: Optional[PathLike] = None,
    fym_fp: Optional[PathLike] = None,
    residue_runs: int = 1,
    spam_crop_raster:  Optional[PathLike] = None,
    irr_yield_scaling: Optional[str] = None,
    spam_all_fp: Optional[PathLike] = None,
    spam_irr_fp: Optional[PathLike] = None,
    spam_rf_fp: Optional[PathLike] = None,
    red_till: bool = False,
    save_CO2: bool = False,
    env_path_overrides: Optional[Dict[str, PathLike]] = None,
    spam_outlier_strategy: str = "spam_sd",
    spam_outlier_percentile: Tuple[float, float] = (5.0, 95.0),
    spam_outlier_k: float = 2.0
):
    def _crop_loader(
        *,
        lu_fp: str,
        evap_fp: str,
        pc_fp: str,
        irr_fp: Optional[str],
        pr_fp: Optional[str],
        fym_fp: Optional[str],
    ) -> Tuple[xr.DataArray, Dict[str, Any]]:
        lu_raster, evap, pc, irr, pr, fym = _load_crop_data(
            _as_path(lu_fp),
            _as_path(evap_fp),
            _as_path(pc_fp),
            None if irr_fp is None else _as_path(irr_fp),
            None if pr_fp  is None else _as_path(pr_fp),
            None if fym_fp is None else _as_path(fym_fp),
        )
        return lu_raster, {
            "evap": evap,
            "pc": pc,
            "irr": irr,
            "c_inp": pr,
            "fym": fym,
        }

    def _crop_runner(
        *,
        n_years: int,
        env: Dict[str, np.ndarray],
        scenario: Dict[str, Any],
        commodity_type: str,
        red_till: bool,
        practices_string_id: Optional[str],
        irr_yield_scaling: Optional[str],
        spam_crop_raster: Optional[str],
        spam_all_fp: Optional[str],
        spam_irr_fp: Optional[str],
        spam_rf_fp: Optional[str],
        spam_outlier_strategy: str | None,
        spam_outlier_percentile: Tuple[float, float] | None,
        spam_outlier_k: float | None
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        evap_a = np.asarray(scenario["evap"].values)
        pc_a = np.asarray(scenario["pc"].values)
        irr_val = scenario.get("irr")
        irr_a = np.asarray(irr_val.values) if hasattr(irr_val, "values") else np.asarray(irr_val) if irr_val is not None else None
        c_inp_a = np.asarray(scenario["c_inp"].values)
        fym_a = np.asarray(scenario["fym"].values)

        base_kwargs = dict(
            n_years=n_years,
            clay=env["clay"],
            soc0=env["soc0"],
            tmp=env["tmp"],
            rain=env["rain"],
            evap=evap_a,
            pc=pc_a,
            c_inp=c_inp_a,
            fym=fym_a,
            commodity_type=commodity_type,
            residue_runs = residue_runs,
            commodity_lu_fp=lu_fp,
            crop_name=crop_name,
            spam_crop_raster = spam_crop_raster,
            practices_string_id = practices_string_id,
            irr_yield_scaling = irr_yield_scaling,
            spam_all_fp = spam_all_fp,
            spam_irr_fp = spam_irr_fp,
            spam_rf_fp = spam_rf_fp,
            spam_outlier_strategy = spam_outlier_strategy,
            spam_outlier_percentile = spam_outlier_percentile,
            spam_outlier_k = spam_outlier_k
        )

        if irr_a is not None:
            base_kwargs["irr"] = irr_a

        base_kwargs["red_till"] = red_till
        if red_till:
            base_kwargs["sand"] = env["sand"]

        return raster_rothc_annual_results(**base_kwargs)
    
    # Creating results_basename
    if commodity_type == "permanent_crop":
        result_basename = f"{crop_name}_{irr_yield_scaling}_{2016+n_years}y_SOC.tif"
    else:
        result_basename = f"{crop_name}_{practices_string_id}_{2016+n_years}y_SOC.tif"

    return _run_rothc_scenario(
        lu_fp=lu_fp,
        n_years=n_years,
        save_folder=save_folder,
        data_description=data_description,
        result_basename=result_basename,
        loader=_crop_loader,
        loader_kwargs={
            "evap_fp": evap_fp,
            "pc_fp": pc_fp,
            "irr_fp": irr_fp,
            "pr_fp": pr_fp,
            "fym_fp": fym_fp,
        },
        runner=_crop_runner,
        runner_kwargs={
            "commodity_type": commodity_type,
            "red_till": red_till,
            "irr_yield_scaling": irr_yield_scaling,
            "practices_string_id": practices_string_id,
            "spam_crop_raster": spam_crop_raster,
            "spam_all_fp": spam_all_fp,
            "spam_irr_fp": spam_irr_fp,
            "spam_rf_fp": spam_rf_fp,
            "spam_outlier_strategy": spam_outlier_strategy,
            "spam_outlier_percentile": spam_outlier_percentile,
            "spam_outlier_k": spam_outlier_k
        },
        loader_message="    Loading crop data...",
        save_CO2=save_CO2,
        env_overrides=env_path_overrides,
    )

# Forest version
def run_RothC_forest(
    forest_type: str,
    weather_type: str,
    n_years: int,
    save_folder: PathLike,
    data_description: str,
    lu_fp: PathLike,
    evap_fp: PathLike,
    age_fp: PathLike,
    practices_string_id: Optional[str] = None,
    TP_IPCC_bool: bool = False,
    save_CO2: bool = False,
    residue_runs = 100,
    env_path_overrides: Optional[Dict[str, PathLike]] = None,
):
    def _forest_loader(
        *,
        lu_fp: str,
        evap_fp: str,
        age_fp: str,
    ) -> Tuple[xr.DataArray, Dict[str, Any]]:
        lu_raster, evap, pc, age = _load_forest_data(
            _as_path(lu_fp),
            _as_path(evap_fp),
            _as_path(age_fp),
        )
        return lu_raster, {
            "evap": evap,
            "pc": pc,
            "age": age,
        }

    def _forest_runner(
        *,
        n_years: int,
        env: Dict[str, np.ndarray],
        scenario: Dict[str, Any],
        forest_type: str,
        weather_type: str,
        TP_IPCC_bool: bool,
        residue_runs: int
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        evap_a = np.asarray(scenario["evap"].values)
        pc_a = np.asarray(scenario["pc"].values)
        age_vals = scenario["age"]
        age_a = np.asarray(age_vals.values) if hasattr(age_vals, "values") else np.asarray(age_vals)
        if age_a.ndim != 2:
            age_a = np.squeeze(age_a)
        if age_a.ndim != 2:
            raise ValueError("Forest age raster must be 2-D after squeezing")

        return raster_rothc_annual_results(
            n_years=n_years,
            clay=env["clay"],
            soc0=env["soc0"],
            tmp=env["tmp"],
            rain=env["rain"],
            evap=evap_a,
            pc=pc_a,
            forest_age=age_a,
            commodity_type="forest",
            forest_type=forest_type,
            weather_type=weather_type,
            TP_IPCC_bool=TP_IPCC_bool,
            residue_runs=residue_runs
        )

    if practices_string_id is not None:
        result_basename = f"{forest_type}_{weather_type}_{practices_string_id}_{2016+n_years}y_SOC.tif"
    else:
        result_basename=f"{forest_type}_{weather_type}_{2016+n_years}y_SOC.tif"

    return _run_rothc_scenario(
        lu_fp=lu_fp,
        n_years=n_years,
        save_folder=save_folder,
        data_description=data_description,
        result_basename=result_basename,
        loader=_forest_loader,
        loader_kwargs={
            "evap_fp": evap_fp,
            "age_fp": age_fp,
        },
        runner=_forest_runner,
        runner_kwargs={
            "forest_type": forest_type,
            "weather_type": weather_type,
            "TP_IPCC_bool": TP_IPCC_bool,
            "residue_runs": residue_runs
        },
        loader_message="    Loading forest data...",
        save_CO2=save_CO2,
        env_overrides=env_path_overrides,
    )


def run_RothC_grassland(
    grassland_type: str,
    n_years: int,
    save_folder: PathLike,
    data_description: str,
    lu_fp: PathLike,
    evap_fp: PathLike,
    fym_fp_list: list[PathLike],
    pc_fp: Optional[PathLike] = None,   # Left for future development of commercial grasslands
    irr_fp: Optional[PathLike] = None,  # Left for future development of commercial grasslands
    pr_fp: Optional[PathLike] = None,   # Left for future development of commercial grasslands
    string_id: Optional[str] = None,
    residue_runs = 100,
    save_CO2: bool = False,
    env_path_overrides: Optional[Dict[str, PathLike]] = None,
):
    def _grassland_loader(
        *,
        lu_fp: Union[str, Path],
        evap_fp: Union[str, Path],
        grassland_type: str,
        fym_fp_list: list[str | Path],
        pc_fp: Optional[Union[str, Path]],
        irr_fp: Optional[Union[str, Path]],
        pr_fp: Optional[Union[str, Path]],
        residue_runs: int,
    ) -> Tuple[xr.DataArray, Dict[str, Any]]:
        lu_raster, evap, pr, pc, fym, irr = _load_grassland_data(
            _as_path(lu_fp),
            _as_path(evap_fp),
            grassland_type,
            [_as_path(fp) for fp in fym_fp_list],
            None if pc_fp is None else _as_path(pc_fp),
            None if irr_fp is None else _as_path(irr_fp),
            None if pr_fp is None else _as_path(pr_fp),
            residue_runs,
        )
        return lu_raster, {
            "evap": evap,
            "pc": pc,
            "irr": irr,
            "c_inp": pr,
            "fym": fym,
        }

    def _grassland_runner(
        *,
        n_years: int,
        env: Dict[str, np.ndarray],
        scenario: Dict[str, Any],
        grassland_type: str,
        grassland_lu_fp: PathLike,
        residue_runs: int,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        
        evap_a = np.asarray(scenario["evap"].values)
        pc_a = np.asarray(scenario["pc"].values)
        irr_val = scenario.get("irr")
        irr_a = np.asarray(irr_val.values) if hasattr(irr_val, "values") else np.asarray(irr_val) if irr_val is not None else None
        c_inp_a = np.asarray(scenario["c_inp"])
        fym_val = scenario.get("fym")
        fym_a = np.asarray(fym_val.values) if hasattr(fym_val, "values") else np.asarray(fym_val)

        base_kwargs = dict(
            n_years=n_years,
            clay=env["clay"],
            soc0=env["soc0"],
            tmp=env["tmp"],
            rain=env["rain"],
            evap=evap_a,
            pc=pc_a,
            c_inp=c_inp_a,
            fym=fym_a,
            commodity_type="grassland",
            commodity_lu_fp=_as_path(grassland_lu_fp),
            grassland_type=grassland_type,
            residue_runs=residue_runs,
        )

        if irr_a is not None:
            base_kwargs["irr"] = irr_a

        return raster_rothc_annual_results(**base_kwargs)

    return _run_rothc_scenario(
        lu_fp=lu_fp,
        n_years=n_years,
        save_folder=save_folder,
        data_description=data_description,
        result_basename=f"{grassland_type}_grassland_{string_id}_{2016+n_years}y_SOC.tif",
        loader=_grassland_loader,
        loader_kwargs={
            "evap_fp": evap_fp,
            "grassland_type": grassland_type,
            "fym_fp_list": fym_fp_list,
            "pc_fp": pc_fp,
            "irr_fp": irr_fp,
            "pr_fp": pr_fp,
            "residue_runs": residue_runs,
        },
        runner=_grassland_runner,
        runner_kwargs={
            "grassland_type": grassland_type,
            "grassland_lu_fp": lu_fp,
            "residue_runs": residue_runs,
        },
        loader_message=f"    Loading {grassland_type} grassland data...",
        save_CO2=save_CO2,
        env_overrides=env_path_overrides,
    )


def DEPRECATED_run_rothC_crop_scenarios_from_csv(csv_filepath: PathLike):
    # 1) Read & cast your CSV exactly as before
    scenarios = (
        pl.read_csv(_resolve_data_path(csv_filepath), null_values=["", "None"])
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
        print(f"Running {scenario['crop_name']} - {scenario['practices_string_id']}")
        run_RothC_crops(**scenario)
        print("\n\n")

def run_rothc_crops_scenarios_from_excel(excel_filepath: PathLike, all_new_files: bool = False, run_test: bool = False, scenario_sheet_name = "scenarios"):
    # 1) Read & cast your CSV exactly as before
    scenarios = (
        pl.read_excel(_resolve_data_path(excel_filepath), has_header=True, sheet_name=scenario_sheet_name)
        .with_columns([
            pl.col("n_years").cast(pl.Int64)
        ])
    )

    if run_test:
        print("Running test. Only top 2 scenarios are run. Residue runs forced to 2, years to 5.")
        scenarios = scenarios[0:2]

    # 2) Turn into a list of dicts once (so we know the total count)
    scenario_list = scenarios.to_dicts()

    # 3) Iterate with tqdm
    for scenario in scenario_list:
        if "force_new_file" in scenario:
            file_fnw = scenario["force_new_file"]

        if scenario["commodity_type"] == "permanent_crop":
            crop_type_string = "Permanent"
        else:
            crop_type_string = "Annual"
        
        scenario_description = (scenario["practices_string_id"] if "practices_string_id" in scenario else scenario    ["irr_yield_scaling"])
        scn_string_text = f"{crop_type_string} crop - {scenario['crop_name']} - {scenario_description}"

        if run_test:
            scenario["residue_runs"] = 2
            scenario["n_years"] = 5

        # Checks if output filepath exist
        output_folder = scenario["save_folder"]
        output_string = f"{scenario['crop_name']}_{scenario_description}_{2016 + scenario['n_years']}y_SOC.tif"
        output_path = f"{output_folder}/{output_string}"

        # Remove 'force_new_file' so it's not forwarded to run_RothC_crops
        scenario.pop("force_new_file", None)

        if all_new_files:
            print(f"Running {scn_string_text}")
            run_RothC_crops(**scenario)
        elif file_fnw is not None and file_fnw is True:
            print(f"Running {scn_string_text}")
            run_RothC_crops(**scenario)
        else:
            if os.path.exists(output_path):
                print(f"{scn_string_text} already exists. Skipping...")
                continue
            else:
                print(f"Running {scn_string_text}")
                run_RothC_crops(**scenario)

        print(f"{scn_string_text} calculated. Continuing...\n\n")


def run_rothc_grassland_scenarios_from_excel(excel_filepath: PathLike, force_new_files: bool = False, run_test: bool = False):
    # 1) Read & cast your CSV exactly as before
    scenarios = (
        pl.read_excel(_resolve_data_path(excel_filepath), has_header=True)
        .with_columns([
            pl.col("n_years").cast(pl.Int64)
        ])
    )

    if run_test:
        print("Running test. Only top 2 scenarios are run.")
        scenarios = scenarios[0:2]

    # 2) Turn into a list of dicts once (so we know the total count)
    scenario_list = scenarios.to_dicts()

    # 3) Iterate with tqdm
    for scenario in scenario_list:
        scn_string_text = f"Grassland - {scenario['grassland_type']} - {scenario['string_id']}"

        # Checks if output filepath exist
        output_folder = scenario["save_folder"]
        output_string = f"{scenario['grassland_type']}_grassland_{scenario['string_id']}_{2016 + scenario['n_years']}y_SOC.tif"
        output_path = f"{output_folder}/{output_string}"

        # Loads fym_fp
        scenario['fym_fp_list'] = json.loads(scenario["fym_fp_list"])

        if force_new_files:
            print(f"Running {scn_string_text}")
            run_RothC_grassland(**scenario)
        else:
            if os.path.exists(output_path):
                print(f"{scn_string_text} already exists. Skipping...")
                continue
            else:
                print(f"Running {scn_string_text}")
                run_RothC_grassland(**scenario)

        print(f"{scn_string_text} calculated. Continuing...\n\n")

def run_rothC_forest_scenarios_from_excel(excel_filepath: PathLike, force_new_files: bool = False, run_test: bool = False):
    # 1) Read & cast your CSV exactly as before
    scenarios = (
        pl.read_excel(_resolve_data_path(excel_filepath), has_header=True)
        .with_columns([
            pl.col("n_years").cast(pl.Int64)
        ])
    )

    if run_test:
        print("Running test. Only top 2 scenarios are run.")
        scenarios = scenarios[0:2]

    # 2) Turn into a list of dicts once (so we know the total count)
    scenario_list = scenarios.to_dicts()

    # 3) Iterate with tqdm
    for scenario in scenario_list:
        scn_string_text = f"Forest - {scenario['forest_type']} - {scenario['weather_type']}"

        # Checks if output filepath exist
        output_folder = scenario["save_folder"]
        output_string = f"{scenario['forest_type']}_{scenario['weather_type']}_{2016 + scenario['n_years']}y_SOC.tif"
        output_path = f"{output_folder}/{output_string}"

        if force_new_files:
            print(f"Running {scn_string_text}")
            run_RothC_forest(**scenario)
        else:
            if os.path.exists(output_path):
                print(f"{scn_string_text} already exists. Skipping...")
                continue
            else:
                print(f"Running {scn_string_text}")
                run_RothC_forest(**scenario)

        print(f"{scn_string_text} calculated. Continuing...\n\n")

##########################################
#### OTHER USEFUL FUNCTIONS FOR ROTHC ####
##########################################

def calcuate_annual_perc_changes(raster_path: PathLike):
    # Open the raster
    da = rxr.open_rasterio(_as_path(raster_path), masked=True)
    if isinstance(da, list):
        da = da[0]

    # Baseline at year 0
    baseline = da.isel(band = 0)

    # Calculate percentage changes
    pct = (da/baseline - 1) * 100

    # Eliminates infinites and replace them with NaNs
    pct = pct.where(np.isfinite(pct))
    
    return pct

def calcuate_practice_change_benefit(
    raster1_fp: PathLike,
    raster2_fp: PathLike,
    band_r1,
    band_r2,
):
    # Open the raster
    da1 = rxr.open_rasterio(_as_path(raster1_fp), masked=True)
    da2 = rxr.open_rasterio(_as_path(raster2_fp), masked=True)

    # Soc at year of band #
    da1_data = da1.isel(band = band_r1)
    da2_data = da2.isel(band = band_r2)

    # Calculate percentage changes
    pct = (da2_data/da1_data - 1) * 100

    # Eliminates infinites and replace them with NaNs
    pct = pct.where(np.isfinite(pct))
    
    return pct
