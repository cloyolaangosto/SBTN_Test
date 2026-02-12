############## Crop Calcs ##############
# Script to calculate yields and residues for different crops

#### MODULES ####
from pathlib import Path
import logging
from functools import lru_cache
import hashlib
from collections import defaultdict
import numpy as np
import pandas as pd
import polars as pl
import xarray as xr
from dataclasses import dataclass
from typing import Mapping, Optional, Tuple, Dict, Union, List, NamedTuple, Set
import geopandas as gpd
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.warp import reproject
import rioxarray as rxr
from scipy import ndimage

from sbtn_leaf.PET import calculate_crop_based_PET_raster_vPipeline
from sbtn_leaf.data_loader import (
    get_crop_coefficients_table,
    get_crop_ag_residue_table,
    get_crop_naming_index_table,
    get_crop_residue_ratio_table,
    get_ecoregions_shapefile,
    get_fao_crop_yields_table,
    get_fao_statistics_table,
    get_country_boundaries,
    get_thermal_climate_tables,
    get_absolute_day_table
)
from sbtn_leaf.paths import data_path
from sbtn_leaf.map_calculations import resample_raster_to_match


# ``Resampling`` is re-exported for backwards compatibility in this module.
from rasterio.enums import Resampling
from rasterio.fill import fillnodata
from scipy import ndimage


@dataclass(frozen=True)
class CropYieldRasterConfig:
    """Configuration for crop yield raster generation.

    SPAM yields are clipped per FAO zone to avoid outliers using FAO statistics.
    """

    fao_avg_yield_name: str
    fao_yield_ratio_name: str
    fao_sd_yield_name: str
    irr_yield_scaling: Optional[str] = None
    all_fp: Optional[str] = None
    irr_fp: Optional[str] = None
    rf_fp: Optional[str] = None
    spam_band: int = 1
    resampling_method: Resampling = Resampling.bilinear
    apply_ecoregion_fill: bool = False
    random_runs: int = 1
    rng: Optional[np.random.Generator] = None
    write_output: bool = True
    return_array: bool = False
    print_outputs: bool = False
    spam_outlier_strategy: str = "spam_sd"
    spam_outlier_percentile: Tuple[float, float] = (1.0, 99.0)
    spam_outlier_k: float = 2.0
    ylds_src: str = "GAEZ"
    p1_window: int = 5     # 5 pixels
    p1_stat: str = "median"
    # pass 2
    p2_enable: bool = True
    p2_window: int = 3     # 3 pixels
    p2_stat: str = "median"
    p2_fallback: str = "clip"  # "clip" or "neighborhood_fill"


@dataclass(frozen=True)
class CropYieldRasterResult:
    """Result object for crop yield raster generation."""

    averaged_result: np.ndarray
    yield_result: np.ndarray
    fao_avg_yields_array: np.ndarray
    fao_sd_yields_array: np.ndarray
    spam_on_lu: np.ndarray
    zone_array: np.ndarray
    lu_meta: Dict[str, object]
    lu_mask: np.ndarray
    lu_transform: Affine
    lu_crs: CRS
    avg_wat_ratio: float
    scaling_mode: Optional[str]


@dataclass(frozen=True)
class IrrigationScalingResult:
    """Result object for irrigation scaling."""

    fao_avg_yields_array: np.ndarray
    avg_wat_ratio: float
    scaling_mode: Optional[str]


##### DATA ####
rain_monthly_fp = data_path("soil_weather", "uhth_monthly_avg_precip.tif")
uhth_climates_fp = data_path("soil_weather", "uhth_thermal_climates.tif")
crop_types      = ["annual", "permanent"]

#### FUNCTIONS ####
def _resolve_crop_coefficient_table(crop_table: Optional[pl.DataFrame] = None) -> pl.DataFrame:
    """Return the provided crop coefficient table or the shared cached copy."""

    if crop_table is not None:
        return crop_table
    return get_crop_coefficients_table()


def _resolve_climate_lookup(
    climate_zone_lookup: Optional[Mapping[int, str]] = None,
) -> Mapping[int, str]:
    """Return the provided climate lookup or the shared cached mapping."""

    if climate_zone_lookup is not None:
        return climate_zone_lookup
    _, lookup, _ = get_thermal_climate_tables()
    return lookup


def _get_crop_naming_index_table() -> pl.DataFrame:
    """Fetch the cached crop naming index table."""

    return get_crop_naming_index_table()


def _get_fao_statistics_table() -> pl.DataFrame:
    """Fetch the cached FAO production statistics table."""

    return get_fao_statistics_table()


def _get_fao_crop_yields_table() -> pl.DataFrame:
    """Fetch the cached FAO crop yields table."""

    return get_fao_crop_yields_table()


def _get_country_boundaries() -> gpd.GeoDataFrame:
    """Fetch the cached country boundary GeoDataFrame."""

    return get_country_boundaries()


def _get_ecoregions_shapefile() -> gpd.GeoDataFrame:
    """Fetch the cached ecoregions GeoDataFrame."""

    return get_ecoregions_shapefile()


def _get_crop_ag_residue_table() -> pl.DataFrame:
    """Fetch the cached crop above-ground residue table."""

    return get_crop_ag_residue_table()


def _get_crop_residue_ratio_table() -> pl.DataFrame:
    """Fetch the cached crop residue ratio table."""

    return get_crop_residue_ratio_table()


def _get_absolute_day_table()-> pl.DataFrame:
    """Fetch the cached cached copy of the absolute day lookup table."""

    return get_absolute_day_table()


# Backwards compatibility: expose lazy proxies for legacy imports expecting
# module-level tables.  The proxies load the underlying dataset on first use
# and then delegate attribute/item access to the cached object.


class _LazyDatasetProxy:
    """Proxy that exposes a lazily loaded dataset via ``__getattr__``/``__getitem__``."""

    def __init__(self, loader):
        self._loader = loader
        self._cached = None

    def _get(self):
        if self._cached is None:
            self._cached = self._loader()
        return self._cached

    def __call__(self):
        return self._get()

    def __getattr__(self, name):
        return getattr(self._get(), name)

    def __getitem__(self, item):
        return self._get()[item]

    def __iter__(self):
        return iter(self._get())


# Legacy attribute names for external callers
crops_name_table = _LazyDatasetProxy(_get_crop_naming_index_table)
fao_stats = _LazyDatasetProxy(_get_fao_statistics_table)
fao_crop_yields_1423 = _LazyDatasetProxy(_get_fao_crop_yields_table)
country_shp = _LazyDatasetProxy(_get_country_boundaries)
crop_ag_res_table = _LazyDatasetProxy(_get_crop_ag_residue_table)
crop_res_table = _LazyDatasetProxy(_get_crop_residue_ratio_table)
er_17 = _LazyDatasetProxy(_get_ecoregions_shapefile)


def index_files(folder_path: str, output_csv: str):
    """
    Walks through `folder_path`, indexes all files, and writes a CSV with:
      - file_name
      - file_path (absolute)
      - suffix (file extension)
    """

    base = Path(folder_path)
    rows = []
    
    for p in base.rglob('*'):   # rglob goes through all existing files
        if p.is_file():         # Checkes if it's actually a file
            rows.append({
                "file_name":     p.name,
                "file_path":     str(p.resolve()),
                "suffix":        p.suffix,
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)

def create_crop_yield_shapefile(
    fao_crop: str,
    *,
    crop_table: Optional[pl.DataFrame] = None,
    yields_table: Optional[pl.DataFrame] = None,
    country_shapes: Optional[gpd.GeoDataFrame] = None,
):
    """Create a country-level shapefile containing FAO yield statistics."""

    crop_table = _get_crop_naming_index_table() if crop_table is None else crop_table
    yields_table = _get_fao_crop_yields_table() if yields_table is None else yields_table
    country_shapes = _get_country_boundaries() if country_shapes is None else country_shapes

    # Checks if the crop is in the list
    if fao_crop not in crop_table['FAO_Crop'].unique():
        raise ValueError(f'{fao_crop} not found or has no data')

    # Extract needed data
    yields_df = (
        yields_table.filter(pl.col("Item") == fao_crop)
        .select(["Area", "Unit", "avg_yield_1423", "ratio_yield_20_toavg", "sd_yields_1423"])
        .to_pandas()
    )

    # rename to shorter names
    yields_df = yields_df.rename(columns={"avg_yield_1423": "avg_yield",
                                          "ratio_yield_20_toavg": "yld_ratio",
                                          "sd_yields_1423": "sd_yield"})

    # Merges with shapefile
    yield_shp = country_shapes.merge(yields_df, how='left', left_on='ADM0_NAME', right_on='Area').drop(columns='Area')

    return yield_shp


def _apply_uncertainty_to_yields(
    result: np.ndarray,
    fao_avg_yields_array: np.ndarray,
    fao_sd_yields_array: np.ndarray,
    lu_mask: np.ndarray,
    *,
    random_runs: int,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Return the average crop yield after sampling FAO uncertainty draws."""

    # Ensure downstream math happens on a predictable floating type.  The
    # baseline input is treated as the deterministic run that will always be
    # included in the ensemble average.
    baseline = np.asarray(result, dtype="float32")

    # Convert FAO reported standard deviations to coefficients of variation
    # (expressed as a fraction of the mean) so we can treat the stochastic
    # component as a percentage delta from the deterministic baseline.  Pixels
    # outside the land-use mask are set to NaN so they do not contribute to any
    # random draws.
    coefficient = np.divide(
        fao_sd_yields_array,
        fao_avg_yields_array,
        out=np.zeros_like(fao_sd_yields_array, dtype="float32"),
        where=fao_avg_yields_array != 0,
    ).astype("float32", copy=False)
    coefficient[~lu_mask] = np.nan

    # When only a deterministic result is requested, return it directly to keep
    # the output shape aligned with the baseline array.
    if random_runs <= 1:
        return baseline

    baseline_valid = baseline[lu_mask]
    coef_valid = coefficient[lu_mask]

    # Generate normally distributed perturbations with a standard deviation set
    # by the coefficient of variation.  The RNG is injectable to support
    # reproducible testing.
    rng = np.random.default_rng() if rng is None else rng
    n_valid = baseline_valid.size
    sum_draws = np.zeros(n_valid, dtype="float32")
    remaining_runs = random_runs - 1
    chunk_size = 1024
    while remaining_runs > 0:
        chunk_runs = min(chunk_size, remaining_runs)
        draws = rng.normal(
            loc=0.0,
            scale=coef_valid,
            size=(chunk_runs, n_valid),
        ).astype("float32", copy=False)
        sum_draws += draws.sum(axis=0, dtype="float32")
        remaining_runs -= chunk_runs

    mean_draws = sum_draws / random_runs
    averaged_valid = baseline_valid * (1.0 + mean_draws)

    # Preallocate result full of NaNs
    averaged = np.full_like(baseline, np.nan, dtype="float32")
    averaged[lu_mask] = averaged_valid

    return averaged.astype("float32", copy=False)


def _read_cropland_raster(
    croplu_grid_raster: str,
) -> tuple[dict, np.ndarray, Affine, CRS, int, int]:
    with rasterio.open(croplu_grid_raster) as crop_lu:
        lu_meta = crop_lu.meta.copy()
        lu_crs = crop_lu.crs
        lu_transform = crop_lu.transform
        lu_height = crop_lu.height
        lu_width = crop_lu.width
        lu_data = crop_lu.read(1)
        lu_nodata = crop_lu.nodata

    lu_mask = (lu_data == 1) & (lu_data != lu_nodata) & (~np.isnan(lu_data))
    return lu_meta, lu_mask, lu_transform, lu_crs, lu_height, lu_width


def _reproject_spam_to_lu(
    spam_crop_raster: str,
    *,
    spam_band: int,
    lu_height: int,
    lu_width: int,
    lu_transform: Affine,
    lu_crs: CRS,
    resampling_method: Resampling,
) -> np.ndarray:
    with rasterio.open(spam_crop_raster) as spam:
        spam_data = spam.read(spam_band)
        spam_on_lu = np.full((lu_height, lu_width), np.nan, dtype="float32")
        reproject(
            source=spam_data,
            destination=spam_on_lu,
            src_transform=spam.transform,
            src_crs=spam.crs,
            src_nodata=spam.nodata,
            dst_transform=lu_transform,
            dst_crs=lu_crs,
            dst_nodata=np.nan,
            resampling=resampling_method,
        )
    return spam_on_lu


def _rasterize_fao_yields(
    fao_crop_shp: gpd.GeoDataFrame,
    lu_crs: CRS,
    lu_height: int,
    lu_width: int,
    lu_transform: Affine,
    fao_avg_yield_name: str,
    fao_yield_ratio_name: str,
    fao_sd_yield_name: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, gpd.GeoDataFrame]:
    """Transform FAO yields from ton/ha to kg/ha and rasterize the shapefile into a raster aligned with the LU raster. 
    Returns fao_average_yields_array, fao_sd_yields_array, zone_array, global_fao_ratio, and fao_gdf_on_lu reprojected to LU raster.

    Args:
        fao_crop_shp (gpd.GeoDataFrame): _description_
        lu_crs (CRS): _description_
        lu_height (int): _description_
        lu_width (int): _description_
        lu_transform (Affine): _description_
        fao_avg_yield_name (str): _description_
        fao_yield_ratio_name (str): _description_
        fao_sd_yield_name (str): _description_

    Raises:
        KeyError: _description_

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, float, gpd.GeoDataFrame]: fao_average_yields_array, fao_sd_yields_array, zone_array, global_fao_ratio, and fao_gdf_on_lu reprojected
    """
    
    
    fao_gdf_on_lu = fao_crop_shp.to_crs(lu_crs).reset_index(drop=True)
    for field in (fao_avg_yield_name, fao_yield_ratio_name, fao_sd_yield_name):
        if field not in fao_gdf_on_lu.columns:
            raise KeyError(f"Missing '{field}' in FAO shapefile")

    fao_gdf_on_lu[fao_avg_yield_name] = fao_gdf_on_lu[fao_avg_yield_name] / 1000.0
    fao_gdf_on_lu[fao_sd_yield_name] = fao_gdf_on_lu[fao_sd_yield_name] / 1000.0
    global_fao_ratio = fao_gdf_on_lu[fao_yield_ratio_name].dropna().mean()

    fao_gdf_on_lu["zone_id"] = fao_gdf_on_lu.index.astype("int32")
    shapes = ((geom, zid) for geom, zid in zip(fao_gdf_on_lu.geometry, fao_gdf_on_lu.zone_id))
    zone_array = rasterize(
        shapes=shapes,
        out_shape=(lu_height, lu_width),
        transform=lu_transform,
        fill=-1,
        dtype="int32",
    )

    fao_avg_yields_array = np.full((lu_height, lu_width), np.nan, dtype="float32")
    fao_sd_yields_array = np.full((lu_height, lu_width), np.nan, dtype="float32")
    for _, row in fao_gdf_on_lu.iterrows():
        zid = int(row["zone_id"])
        zid_mask = zone_array == zid
        fao_avg_yields_array[zid_mask] = row[fao_avg_yield_name]
        fao_sd_yields_array[zid_mask] = row[fao_sd_yield_name]

    return (
        fao_avg_yields_array,
        fao_sd_yields_array,
        zone_array,
        global_fao_ratio,
        fao_gdf_on_lu,
    )


def _apply_irrigation_scaling_toFAO_yields(
    fao_avg_yields_array: np.ndarray,
    valid_fao: np.ndarray,
    irr_yield_scaling: str,
    *,
    spam_all: np.ndarray,
    spam_irr: np.ndarray,
    spam_rf: np.ndarray,
    print_outputs: bool,
) -> IrrigationScalingResult:
    scaling_mode = irr_yield_scaling.lower()
    if scaling_mode not in {"irr", "rf"}:
        raise ValueError("irr_yield_scaling must be either 'irr' or 'rf'")
    if any(path is None for path in (spam_all, spam_irr, spam_rf)):
        raise ValueError("Need all_fp, irr_fp and rf_fp for irrigation scaling")

    irr_ratios, rf_ratios = _calculate_SPAM_yield_modifiers(
        all_yields=spam_all,
        irr_yields=spam_irr,
        rf_yields=spam_rf,
        print_outputs=print_outputs,
    )

    watering_ratio = irr_ratios if scaling_mode == "irr" else rf_ratios

    valid_wat = ~np.isnan(watering_ratio)
    avg_wat_ratio = np.nanmean(watering_ratio)

    # Creates an irrigations scale yield raster applying watering ratios 
    fao_scaled = np.where(valid_wat, fao_avg_yields_array * watering_ratio, np.nan)
    
    # Fill close pixels that are nan
    fao_scaled = fillnodata(
        fao_scaled, mask=np.isnan(fao_scaled), max_search_distance=1, smoothing_iterations=2
    )

    # Applies a global watering ratio to pixels where fao_scaled is nan but fao yields are valid
    fao_scaled = np.where(
        np.isnan(fao_scaled) & valid_fao,
        fao_avg_yields_array * avg_wat_ratio,
        fao_scaled,
    )

    return IrrigationScalingResult(
        fao_avg_yields_array=fao_scaled,
        avg_wat_ratio=avg_wat_ratio,
        scaling_mode=scaling_mode,
    )


def _compose_yield_result(
    spam_on_lu: np.ndarray,
    fao_gdf: gpd.GeoDataFrame,
    zone_array: np.ndarray,
    fao_avg_yields_array: np.ndarray,
    lu_mask: np.ndarray,
    valid_fao: np.ndarray,
    global_fao_ratio: float,
    fao_yield_ratio_name: str,
    avg_wat_ratio: float,
    *,
    all_fp_on_lu: Optional[np.ndarray],
    scaling_mode: Optional[str],
    print_outputs: bool,
    ylds_src: str = "GAEZ"
) -> np.ndarray:
    """Compose the yield raster while clipping SPAM outliers per FAO zone."""
    # Creates an empty array of yields with the same size as fao yields raster
    result = np.full_like(fao_avg_yields_array, np.nan, dtype="float32")
    
    # Fill the array
    for _, row in fao_gdf.iterrows():
        zid = int(row["zone_id"])  # Gets zone id
        fao_yield_adjustment_ratio = 1 if (ylds_src == "GAEZ") else row[fao_yield_ratio_name]
        
        zid_mask = zone_array == zid  # Creates a mask to apply only for the given fao zone
        spam_scaled = spam_on_lu * fao_yield_adjustment_ratio  # Loads the spam yields
        
        # Fills results with spam scaled factors
        valid_mask = zid_mask & ~np.isnan(spam_scaled)
        result[valid_mask] = spam_scaled[valid_mask]

        # Checks where there are not results for spam values
        mask_need_avg = zid_mask & np.isnan(result)
        result[mask_need_avg] = fao_avg_yields_array[mask_need_avg] * avg_wat_ratio

    # Fill results where there are still missing pixels
    if all_fp_on_lu is not None and scaling_mode is not None:
        label = "rainfed" if scaling_mode == "rf" else "irrigation"
        mask_missing = lu_mask & np.isnan(result) & ~np.isnan(all_fp_on_lu)
        
        if np.any(mask_missing):
            if print_outputs:
                print(f"Pixels still missing values...  → Applying {label} scaling to all‐SPAM yields…")
            result[mask_missing] = spam_on_lu[mask_missing] * global_fao_ratio

    # If there are still missing pixels, fills them with fao yields
    mask_fao = lu_mask & np.isnan(result) & valid_fao
    if np.any(mask_fao):
        if print_outputs:
            print(f"Pixels still missing values... Applying FAO scaled yields")
        result[mask_fao] = fao_avg_yields_array[mask_fao]

    result_mean = np.nanmean(result)
    result_median = np.nanmedian(result)
    result_max = np.nanmax(result)
    result_min = np.nanmin(result)
    
    if print_outputs:
        print(f"Current mean and median yields are {result_mean} and {result_median}. Max is {result_max} and min is {result_min}")

    return result


def _fill_with_ecoregions(
    result: np.ndarray,
    croplu_grid_raster: str,
    lu_mask: np.ndarray,
    global_fao_ratio: float,
) -> np.ndarray:
    ecoregion_avg, biome_avg, zone_array, biome_name_map = calculate_average_yield_by_ecoregion_and_biome(
        result, croplu_grid_raster
    )

    remaining = lu_mask & np.isnan(result)
    if np.any(remaining):
        zone_max = int(zone_array.max())
        if zone_max >= 0:
            zone_lookup = np.full(zone_max + 1, np.nan, dtype=float)
            for zid, avg in ecoregion_avg.items():
                if 0 <= zid <= zone_max:
                    zone_lookup[zid] = avg

            unique_zones = np.unique(zone_array)
            unique_zones = unique_zones[unique_zones >= 0]
            if unique_zones.size:
                missing = np.isnan(zone_lookup[unique_zones])
                if np.any(missing):
                    for zid in unique_zones[missing]:
                        biome = biome_name_map.get(int(zid))
                        if isinstance(biome, str):
                            zone_lookup[zid] = biome_avg.get(biome, global_fao_ratio)
                        else:
                            zone_lookup[zid] = global_fao_ratio
            zone_lookup = np.where(np.isnan(zone_lookup), global_fao_ratio, zone_lookup)

            remaining_zones = zone_array[remaining]
            fill_vals = np.full(remaining_zones.shape, global_fao_ratio, dtype=float)
            valid_zones = remaining_zones >= 0
            if np.any(valid_zones):
                fill_vals[valid_zones] = zone_lookup[remaining_zones[valid_zones]]
            result[remaining] = fill_vals
        else:
            result[remaining] = global_fao_ratio

    remaining = lu_mask & np.isnan(result)
    if np.any(remaining):
        valid = ~np.isnan(result)
        dist, (iy, ix) = ndimage.distance_transform_edt(
            ~valid, return_distances=True, return_indices=True
        )
        filled = result[iy, ix]
        result[remaining] = filled[remaining]

    return result

def _pre_filter_yields_rasters(
    yield_arrays: tuple[np.ndarray, ...],
    fao_gdf: gpd.GeoDataFrame,
    zone_array: np.ndarray,
    fao_avg_yield_name: str,
    *,
    filter_outlier_strategy: str,
    percentile_bounds: Tuple[float, float] = (1.0, 99.0),
    k_sd: float,
    # pass 1
    p1_window: int = 5,     # 5 pixels
    p1_stat: str = "median",
    # pass 2
    p2_enable: bool = True,
    p2_window: int = 3,     # 3 pixels
    p2_stat: str = "median",
    p2_fallback: str = "clip",  # "clip" or "neighborhood_fill"
):
    #  Initialize out
    out_arrays = [np.full_like(a, np.nan, dtype="float32") for a in yield_arrays]

    # Goes through each FAO Zone (Country)
    for _, row in fao_gdf.iterrows():
        zid = int(row["zone_id"])                       #  Gets the zone id
        zid_mask = zone_array == zid                    #  Mask of the country
        fao_zone_avg = row[fao_avg_yield_name]          #  Country yield average (10-year)

        # For each array:
        for i, array in enumerate(yield_arrays):
            # --- start zone subset ---
            valid_zone = zid_mask & np.isfinite(array)
            yld_vals = array[valid_zone]

            # Check if there are any values
            if yld_vals.size == 0:
                continue

            # ========== PASS 1 ==========
            min1, max1 = _compute_zone_bounds(
                yld_vals,
                strategy=filter_outlier_strategy,
                fao_zone_avg=fao_zone_avg,
                pct=percentile_bounds,
                k_sd=k_sd,
            )
            if not (np.isfinite(min1) and np.isfinite(max1)):
                continue

            out1 = array.astype("float32", copy=True)
            outlier1 = valid_zone & ((out1 < min1) | (out1 > max1))

            if np.any(outlier1):
                out1 = _replace_masked_with_neighborhood(
                    out1,
                    replace_mask=outlier1,
                    valid_mask=valid_zone,
                    window=p1_window,
                    replace_with=p1_stat,
                )
            # guarantee bounds after fill (prevents a few corner cases)
            out1 = np.where(valid_zone, np.clip(out1, min1, max1), out1)

            # ========== PASS 2 (optional) ==========
            if p2_enable:
                # recompute on pass1 output (still per-zone)
                vals2 = out1[valid_zone]
                min2, max2 = _compute_zone_bounds(
                    vals2,
                    strategy=filter_outlier_strategy,
                    fao_zone_avg=fao_zone_avg,
                    pct=percentile_bounds,
                    k_sd=k_sd,
                )

                if np.isfinite(min2) and np.isfinite(max2):
                    outlier2 = valid_zone & ((out1 < min2) | (out1 > max2))
                    if np.any(outlier2):
                        if p2_fallback == "neighborhood_fill":
                            out1 = _replace_masked_with_neighborhood(
                                out1,
                                replace_mask=outlier2,
                                valid_mask=valid_zone,
                                window=p2_window,
                                replace_with=p2_stat,
                            )
                            out1 = np.where(valid_zone, np.clip(out1, min2, max2), out1)
                        elif p2_fallback == "clip":
                            out1 = np.where(outlier2, np.clip(out1, min2, max2), out1)
                        else:
                            raise ValueError("p2_fallback must be 'clip' or 'neighborhood_fill'")

            # write only inside valid_zone for this FAO zone
            out_arrays[i][valid_zone] = out1[valid_zone]

    return out_arrays


def _compute_zone_bounds(
    yld_vals: np.ndarray,
    *,
    strategy: str,
    fao_zone_avg: float,
    pct: tuple[float, float],
    k_sd: float,
) -> tuple[float, float]:
    if strategy == "ratio_percentile":
        if not np.isfinite(fao_zone_avg) or fao_zone_avg <= 0:
            return (np.nan, np.nan)
        ratios = yld_vals / fao_zone_avg
        lo_r, hi_r = np.nanpercentile(ratios, [pct[0], pct[1]])
        return (lo_r * fao_zone_avg, hi_r * fao_zone_avg)

    if strategy == "spam_sd":
        mu = np.nanmean(yld_vals)
        sd = np.nanstd(yld_vals)
        return (max(0, mu - k_sd * sd), (mu + k_sd * sd))

    raise ValueError(f"Unknown strategy: {strategy}")


def _fill_remaining_with_nearest_in_zone(arr, replace_mask, valid_mask):
    """Fill replace_mask pixels using nearest valid pixel *within valid_mask*."""
    out = arr.astype("float32", copy=True)

    # valid donors are pixels that are valid_mask and NOT being replaced and finite
    donors = valid_mask & ~replace_mask & np.isfinite(out)

    if not np.any(donors):
        return out  # nothing we can do

    # distance transform: for pixels where donors==False, returns indices of nearest True in donors
    dist, (iy, ix) = ndimage.distance_transform_edt(~donors, return_indices=True)

    # fill only where requested and donor exists
    out[replace_mask] = out[iy[replace_mask], ix[replace_mask]]
    return out

def _replace_masked_with_neighborhood(
    arr: np.ndarray,
    replace_mask: np.ndarray,
    valid_mask: np.ndarray,
    *,
    window: int = 5,
    replace_with: str = "median",
    mode: str = "nearest",
    fallback_nearest: bool = True,
) -> np.ndarray:
    if window < 3 or window % 2 == 0:
        raise ValueError("window must be odd and >= 3")

    work = arr.astype("float32", copy=True)
    work[~valid_mask] = np.nan
    work[replace_mask] = np.nan

    footprint = np.ones((window, window), dtype=bool)

    if replace_with == "median":
        neigh = ndimage.generic_filter(work, np.nanmedian, footprint=footprint, mode=mode)
    elif replace_with == "mean":
        neigh = ndimage.generic_filter(work, np.nanmean, footprint=footprint, mode=mode)
    else:
        raise ValueError("replace_with must be 'median' or 'mean'")

    out = arr.astype("float32", copy=True)

    can_replace = replace_mask & np.isfinite(neigh)
    out[can_replace] = neigh[can_replace]

    # fallback: if neighborhood stat couldn't find any donors (NaN), use nearest in-country donor
    remaining = replace_mask & ~np.isfinite(neigh)
    if fallback_nearest and np.any(remaining):
        out = _fill_remaining_with_nearest_in_zone(out, remaining, valid_mask)

    return out

def _create_crop_yield_raster_core(
    croplu_grid_raster: str,
    fao_crop_shp: gpd.GeoDataFrame,
    spam_crop_raster: str,
    output_rst_path: Optional[str],
    config: CropYieldRasterConfig,
) -> CropYieldRasterResult:
    """Shared implementation for the crop yield raster generators.

    When ``config.return_array`` is ``True``, raster output is skipped and the
    caller is expected to use the arrays in the returned result instead.
    """
    write_output = config.write_output
    if config.return_array:
        if write_output:
            logging.info(
                "return_array=True requested; skipping raster write for crop yield results."
            )
        write_output = False
    
    # Step 1 - Read cropland raster
    (
        lu_meta,
        lu_mask,
        lu_transform,
        lu_crs,
        lu_height,
        lu_width,
    ) = _read_cropland_raster(croplu_grid_raster)

    # Step 2 - Reprojects all spam bands
    spam_on_lu = _reproject_spam_to_lu(
        spam_crop_raster,
        spam_band=config.spam_band,
        lu_height=lu_height,
        lu_width=lu_width,
        lu_transform=lu_transform,
        lu_crs=lu_crs,
        resampling_method=config.resampling_method,
    )

    spam_all = _reproject_spam_to_lu(
        spam_crop_raster=config.all_fp,
        spam_band=config.spam_band,
        lu_height=lu_height,
        lu_width=lu_width,
        lu_transform=lu_transform,
        lu_crs=lu_crs,
        resampling_method=config.resampling_method,
    )

    spam_irr = _reproject_spam_to_lu(
        spam_crop_raster=config.irr_fp,
        spam_band=config.spam_band,
        lu_height=lu_height,
        lu_width=lu_width,
        lu_transform=lu_transform,
        lu_crs=lu_crs,
        resampling_method=config.resampling_method,
    )

    spam_rf = _reproject_spam_to_lu(
        spam_crop_raster=config.rf_fp,
        spam_band=config.spam_band,
        lu_height=lu_height,
        lu_width=lu_width,
        lu_transform=lu_transform,
        lu_crs=lu_crs,
        resampling_method=config.resampling_method,
    )

    # Step 3 - Rasterize fao yields
    (
        fao_avg_yields_array,
        fao_sd_yields_array,
        zone_array,
        global_fao_ratio,
        fao_gdf,
    ) = _rasterize_fao_yields(
        fao_crop_shp,
        lu_crs,
        lu_height,
        lu_width,
        lu_transform,
        config.fao_avg_yield_name,
        config.fao_yield_ratio_name,
        config.fao_sd_yield_name,
    )

    valid_fao = ~np.isnan(fao_avg_yields_array)
    irrigation_scaling = IrrigationScalingResult(
        fao_avg_yields_array=fao_avg_yields_array,
        avg_wat_ratio=np.nan,
        scaling_mode=None,
    )

    #  Pre process SPAM yields
    spam_arrays = (spam_on_lu, spam_all, spam_irr, spam_rf)
   
    (spam_on_lu_filt, spam_all_filt, spam_irr_filt, spam_rf_filt) = _pre_filter_yields_rasters(
        yield_arrays=spam_arrays,
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name=config.fao_avg_yield_name,
        filter_outlier_strategy=config.spam_outlier_strategy,
        percentile_bounds=config.spam_outlier_percentile,
        k_sd=config.spam_outlier_k,
        # pass 1
        p1_window = config.p1_window,     # 5 pixels
        p1_stat = config.p1_stat,
        # pass 2
        p2_enable = config.p2_enable,
        p2_window = config.p2_window,
        p2_stat= config.p2_stat,
        p2_fallback= config.p2_fallback,  # "clip" or "neighborhood_fill"
    )

    # Apply yields scaling based on irrigation technique and SPAM yields
    if config.irr_yield_scaling is not None:
        irrigation_scaling = _apply_irrigation_scaling_toFAO_yields(
            fao_avg_yields_array,
            valid_fao,
            config.irr_yield_scaling,
            spam_all=spam_all_filt,
            spam_irr=spam_irr_filt,
            spam_rf=spam_rf_filt,
            print_outputs=config.print_outputs,
        )

    # Prepare results
    result = _compose_yield_result(
        spam_on_lu_filt,
        fao_gdf,
        zone_array,
        irrigation_scaling.fao_avg_yields_array,
        lu_mask,
        valid_fao,
        global_fao_ratio,
        config.fao_yield_ratio_name,
        all_fp_on_lu=spam_on_lu_filt,
        avg_wat_ratio=irrigation_scaling.avg_wat_ratio,
        scaling_mode=irrigation_scaling.scaling_mode,
        print_outputs=config.print_outputs,
        ylds_src = config.ylds_src
    )

    if config.apply_ecoregion_fill:
        result = _fill_with_ecoregions(
            result,
            croplu_grid_raster,
            lu_mask,
            global_fao_ratio,
        )

    result[~lu_mask] = np.nan

    # Apply uncertainty to results
    averaged_result = _apply_uncertainty_to_yields(
        result,
        irrigation_scaling.fao_avg_yields_array,
        fao_sd_yields_array,
        lu_mask,
        random_runs=config.random_runs,
        rng=config.rng,
    )

    mean = np.nanmean(averaged_result)
    median = np.nanmedian(averaged_result)

    print(f"        Final mean is {mean:.1f} and median is {median:.1f}.")

    # Output block
    if write_output:
        _write_yield_raster(output_rst_path, lu_meta, averaged_result)

    return CropYieldRasterResult(
        averaged_result=averaged_result,
        yield_result=result,
        fao_avg_yields_array=irrigation_scaling.fao_avg_yields_array,
        fao_sd_yields_array=fao_sd_yields_array,
        spam_on_lu=spam_on_lu,
        zone_array=zone_array,
        lu_meta=lu_meta,
        lu_mask=lu_mask,
        lu_transform=lu_transform,
        lu_crs=lu_crs,
        avg_wat_ratio=irrigation_scaling.avg_wat_ratio,
        scaling_mode=irrigation_scaling.scaling_mode,
    )


def _write_yield_raster(
    output_rst_path: Optional[str],
    lu_meta: Dict[str, object],
    averaged_result: np.ndarray,
) -> None:
    """Write the averaged yield raster to disk."""

    if output_rst_path is None:
        raise ValueError("output_rst_path is required when write_output is True")

    lu_meta.update(dtype="float32", count=1, nodata=np.nan)
    with rasterio.open(output_rst_path, "w", **lu_meta) as dst:
        dst.write(averaged_result[np.newaxis, ...])

    print(f"Yield raster written to {output_rst_path}")


def create_crop_yield_raster(
    croplu_grid_raster: str,
    fao_crop_shp: gpd.GeoDataFrame,
    spam_crop_raster: str,
    output_rst_path: str,
    spam_band: int = 1,
    resampling_method: Resampling = Resampling.bilinear,
    ylds_src: str = "GAEZ"
) -> CropYieldRasterResult:
    """Create a crop yield raster without irrigation scaling."""

    config = CropYieldRasterConfig(
        fao_avg_yield_name="avg_yield_1423",
        fao_yield_ratio_name="ratio_yield_20_toavg",
        fao_sd_yield_name="sd_yields_1423",
        spam_band=spam_band,
        resampling_method=resampling_method,
        print_outputs=True,
        ylds_src = ylds_src
    )
    return _create_crop_yield_raster_core(
        croplu_grid_raster,
        fao_crop_shp,
        spam_crop_raster,
        output_rst_path,
        config,
    )


def create_crop_yield_raster_withIrrigationPracticeScaling(
    croplu_grid_raster: str,
    fao_crop_shp: "gpd.GeoDataFrame",
    spam_crop_raster: str,
    output_rst_path: str,
    spam_band: int = 1,
    resampling_method: Resampling = Resampling.bilinear,
    irr_yield_scaling: Optional[str] = None,
    all_fp: Optional[str] = None,
    irr_fp: Optional[str] = None,
    rf_fp: Optional[str] = None,
    fao_avg_yield_name: str = "avg_yield",
    fao_yield_ratio_name: str = "yld_ratio",
    fao_sd_yield_name: str = "sd_yield",
    apply_ecoregion_fill: bool = True,
    spam_outlier_strategy: str = "spam_sd",
    spam_outlier_percentile: Tuple[float, float] = (0.5, 0.95),
    spam_outlier_k: float = 2,
    p1_window: int = 5,     # 5 pixels
    p1_stat: str = "median",
    # pass 2
    p2_enable: bool = True,
    p2_window: int = 3,     # 3 pixels
    p2_stat: str = "median",
    p2_fallback: str = "clip"  # "clip" or "neighborhood_fill"
) -> CropYieldRasterResult:
    """Create a crop yield raster with optional irrigation/rainfed scaling.

    Parameters
    ----------
    irr_yield_scaling:
        Either ``"irr"`` or ``"rf"`` to select irrigated or rainfed scaling. ``None``
        preserves the FAO averages.
    apply_ecoregion_fill:
        When ``True`` (the default), use ecoregion and biome averages to fill any
        remaining nodata pixels, matching the historical pipeline behaviour.
    """

    config = CropYieldRasterConfig(
        fao_avg_yield_name=fao_avg_yield_name,
        fao_yield_ratio_name=fao_yield_ratio_name,
        fao_sd_yield_name=fao_sd_yield_name,
        irr_yield_scaling=irr_yield_scaling,
        all_fp=all_fp,
        irr_fp=irr_fp,
        rf_fp=rf_fp,
        spam_band=spam_band,
        resampling_method=resampling_method,
        apply_ecoregion_fill=apply_ecoregion_fill,
        print_outputs=True,
        spam_outlier_strategy = spam_outlier_strategy,
        spam_outlier_percentile = spam_outlier_percentile,
        spam_outlier_k = spam_outlier_k,
        p1_window =p1_window,
        p1_stat = p1_stat,
        # pass 2
        p2_enable = p2_enable,
        p2_window = p2_window,
        p2_stat = p2_stat,
        p2_fallback = p2_fallback,
    )
    return _create_crop_yield_raster_core(
        croplu_grid_raster,
        fao_crop_shp,
        spam_crop_raster,
        output_rst_path,
        config,
    )


def _calculate_SPAM_yield_modifiers(
    all_yields: np.ndarray,
    irr_yields: np.ndarray,
    rf_yields: np.ndarray,
    save_ratios: bool = False,
    all_rasters_fp: Optional[str] = None,
    irr_ratios_fp: Optional[str] = None,
    rf_ratios_fp: Optional[str] = None,
    print_outputs: bool = False
):
    '''
    Calculates the the ratios yields between irrigated:all and rainfed:all
    '''
    # Step 1 - Open all files
    if save_ratios:
        with rasterio.open(all_rasters_fp) as all_src:
            all_profile   = all_src.profile.copy()

    # Step 2 - Creates a new array for irrigation and rainfed
    irr_ratios = np.full_like(all_yields, fill_value=np.nan, dtype=float)
    rf_ratios = np.full_like(all_yields, fill_value=np.nan, dtype=float)

    # Step 3 - Creates mask for where is data
    all_mask = ~np.isnan(all_yields) & (all_yields != 0)
    irr_mask = ~np.isnan(irr_yields) & (irr_yields != 0)
    rf_mask  = ~np.isnan(rf_yields)  & (rf_yields  != 0)

    # Step 4 - Fills ratios array where both mask are true

    # only divide where both the overall and irrigation/rainfed masks are True
    np.divide(
        irr_yields,
        all_yields,
        out=irr_ratios,
        where=(all_mask & irr_mask)
    )

    np.divide(
        rf_yields,
        all_yields,
        out=rf_ratios,
        where=(all_mask & rf_mask)
    )

    # compute means safely and print using str.format to avoid f-string parsing issues
    avg_irr = np.nanmean(irr_ratios)
    avg_rf = np.nanmean(rf_ratios)
    
    if print_outputs:
        print("Average irrigated ratio: {:.2f}".format(avg_irr))
        print("Average rainfed ratio: {:.2f}".format(avg_rf))

    # Step 5 - Optional, Save as GeoTiff
    if save_ratios:
        # Checks if save path have been provided
        if (irr_ratios_fp is None) or (rf_ratios_fp is None) or (all_rasters_fp is None):
            print('Source path for all yields or Saving path not provided for irrigation or rainfed ratios. Skipping save...')

        else:
            with rasterio.open(all_rasters_fp) as all_src:
                all_profile   = all_src.profile.copy()
            
            # Updates profiles
            irr_profile = all_profile.copy()
            irr_profile.update(
                dtype='float32',
                count=1,
                nodata=np.nan,
                description='Irrigated to Irrigated+Rainfed SPAM yield ratio'
            )

            rf_profile = all_profile.copy()
            rf_profile.update(
                dtype='float32',
                count=1,
                nodata=np.nan,
                description='Rainfed to Irrigated+Rainfed SPAM yield ratio'
            )

            # Writing the GeoTiffs
            with rasterio.open(irr_ratios_fp, "w", **irr_profile) as dst_irr:
                dst_irr.write(irr_ratios.astype("float32"), 1)  # Must write data first
                dst_irr.update_tags(
                    model="SPAM",
                    scenario="irrigated",
                    units="ratio",
                    description="Irrigated yields ratios compared to all yields"
                )

            with rasterio.open(rf_ratios_fp, "w", **rf_profile) as dst_rf:
                dst_rf.write(rf_ratios.astype("float32"), 1)  # Must write data first
                dst_rf.update_tags(
                    model="SPAM",
                    scenario="irrigated",
                    units="ratio",
                    description="Rainfed yield ratios compared to all yields"
                )

    # Returning ratios
    return irr_ratios, rf_ratios


def calculate_average_yield_by_ecoregion_and_biome(
    result_arr: np.ndarray,
    croplu_grid_raster: str,
    er_shapefile: Optional[gpd.GeoDataFrame] = None,
) -> Tuple[Dict[int, float], Dict[str, float], np.ndarray, Dict[int, str]]:
    """
    Calculate average yields per ecoregion and per biome.
    Returns:
      - ecoregion_avg: zone_id -> average yield
      - biome_avg: biome_name -> average yield
      - zone_array: rasterized zone_id array
      - biome_name_map: zone_id -> biome_name mapping
    """
    er_shapefile = _get_ecoregions_shapefile() if er_shapefile is None else er_shapefile

    # load LU raster for transform & CRS
    with rasterio.open(croplu_grid_raster) as src:
        transform, crs = src.transform, src.crs
        height, width = src.height, src.width

    # load and reproject ecoregions
    er_gdf = er_shapefile.to_crs(crs).reset_index(drop=True)
    er_gdf["zone_id"] = er_gdf.index.astype("int32")

    # rasterize zone_id
    shapes = ((geom, zid) for geom, zid in zip(er_gdf.geometry, er_gdf.zone_id))
    zone_array = rasterize(
        shapes,
        out_shape=(height, width),
        transform=transform,
        fill=-1,
        dtype="int32"
    )

    # compute per-zone averages
    ecoregion_avg: Dict[int, float] = {}
    biome_pixels: Dict[str, list] = {}
    biome_name_map: Dict[int, str] = {}

    for _, row in er_gdf.iterrows():
        zid = row["zone_id"]
        mask = (zone_array == zid)
        vals = result_arr[mask]
        vals = vals[~np.isnan(vals)]
        if vals.size:
            ecoregion_avg[zid] = float(vals.mean())
        biome = row["BIOME_NAME"]
        biome_name_map[zid] = biome
        biome_pixels.setdefault(biome, []).append(mask)

    # biome averages
    biome_avg: Dict[str, float] = {}
    for biome, masks in biome_pixels.items():
        combined = np.zeros_like(result_arr, dtype=bool)
        for m in masks:
            combined |= m
        vals = result_arr[combined]
        vals = vals[~np.isnan(vals)]
        if vals.size:
            biome_avg[biome] = float(vals.mean())

    return ecoregion_avg, biome_avg, zone_array, biome_name_map


def calculate_crop_residues(crop: str, crop_yield: float, C_Content: float = 0.5):
    """
    Compute above- and below-ground residues (in C-content dry matter) for a given crop.

    Returns a dict with keys:
      - 'ABG': above-ground residue C-mass
      - 'BG':  below-ground residue C-mass
      - 'Total': sum of ABG + BG
    """
    
    crop_table = _get_crop_naming_index_table()
    res_table = _get_crop_residue_ratio_table()
    ag_table = _get_crop_ag_residue_table()

    if crop not in (crop_table["Crops"].to_list() + crop_table["IPCC_Crop"].to_list()):
        raise ValueError(f"{crop} not found in data table")
    
    # Initialize crop amounts
    ABG = 0.0
    BG = 0.0
    Res = ABG + BG

    # Getting IPCC name and calculating belowground
    ipcc_crop = crop_table.filter(pl.col('Crops') == crop).select('IPCC_Crop').item()
    res_crop_data = res_table.filter(pl.col('Crop') == ipcc_crop)
    dry = res_crop_data.select('DRY').item()
    RS = res_crop_data.select("RS").item()

    # Checks if ABG can be calculated with line equation
    if crop in ag_table["Crop"].to_list():
        AG_crop_data = ag_table.filter(pl.col("Crop") == crop)
        slope = AG_crop_data.select("Slope").item()
        intercept = AG_crop_data.select("Intercept").item()

        # Calcuating the plant residues
        ABG = slope * crop_yield + intercept
        BG = float(RS) * ABG 
    # If not, checks if it has above ground to yield ratios
    elif RS > 0:
        R_AG = res_crop_data.select("R_AG").item()
        ABG = crop_yield * R_AG
        BG = float(RS) * ABG 
    # If not, goes through total yield to total residues
    else:
        Res = res_crop_data.select("R_T").item() * crop_yield * dry * C_Content

    # Now translating into dry matter and carbon content
    if ABG > 0:
        ABG = ABG * dry * C_Content
        BG = BG  * dry * C_Content
        Res = ABG + BG
    else:
        ABG = np.nan
        BG = np.nan

    # Returning results
    return {
        'Res': Res, 
        'ABG': ABG, 
        'BG': BG
        }


def apply_residues_to_raster_flexible(
    crop: str,
    yield_raster_path: str,
    output_path: str,
    C_Content: float = 0.5,
    band: int = 1,
):
    """
    Reads a single‐band yield raster, applies calculate_crop_residues() to each pixel,
    and then writes out either:
      - a 1‐band TIFF with only 'Res' if ABG/BG are always NaN, or
      - a 3‐band TIFF with 'Res','ABG','BG' if all three are valid.
    """

    # 1) Load input raster
    with rasterio.open(yield_raster_path) as src:
        meta      = src.meta.copy()
        yields    = src.read(1).astype("float32")
        nodata    = src.nodata

    # mask out input nodata → NaN
    mask       = (yields == nodata)
    yields[mask]  = np.nan

    # 2) Sample the function at yield=1 to see which outputs are real
    sample = calculate_crop_residues(crop, 1.0, C_Content)
    has_abg = not np.isnan(sample["ABG"])

    # 3) Build vectorized functions for each needed output
    vec_res = np.vectorize(
        lambda y: calculate_crop_residues(crop, float(y), C_Content)["Res"],
        otypes=["float32"]
    )

    if has_abg:
        vec_abg = np.vectorize(
            lambda y: calculate_crop_residues(crop, float(y), C_Content)["ABG"],
            otypes=["float32"]
        )
        vec_bg = np.vectorize(
            lambda y: calculate_crop_residues(crop, float(y), C_Content)["BG"],
            otypes=["float32"]
        )

    # 4) Apply them to the full array
    res_arr = vec_res(yields)
    res_arr[mask] = np.nan

    if has_abg:
        abg_arr = vec_abg(yields)
        bg_arr  = vec_bg(yields)
        abg_arr[mask] = np.nan
        bg_arr[mask]  = np.nan

    # 5) Write output(s)
    if has_abg:
        # 3‐band output
        meta.update(count=3, dtype="float32", nodata=np.nan)
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(res_arr,  1)  # band 1 = Res
            dst.write(abg_arr,  2)  # band 2 = ABG
            dst.write(bg_arr,   3)  # band 3 = BG
        print(f"Wrote 3‐band raster to {output_path} (Res, ABG, BG)")
    else:
        # single‐band output (Res only)
        meta.update(count=1, dtype="float32", nodata=np.nan)
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(res_arr, 1)
        print(f"Wrote 1‐band raster to {output_path} (only Res)")


def create_residue_raster_rasterops(
    crop: str,
    yield_raster_path: str,
    output_path: str,
    C_Content: float = 0.50
):
    """
    Read a single‐band yield raster, choose one residue‐calculation path globally:
      1) Regression (Slope/Intercept) if available
      2) Else Ratio (R_AG, RS) if RS > 0
      3) Else Total‐residues (R_T, dry, C_Content)
    and write either:
      • a 3-band TIFF (Res, ABG, BG) for paths 1 & 2, or
      • a 1-band TIFF (Res only) for path 3.
    """

    # 1) Load the yield raster
    with rasterio.open(yield_raster_path) as src:
        meta   = src.meta.copy()
        yields    = src.read(1).astype("float32")
        nodata = src.nodata

    # mask nodata → NaN
    valid = (yields != nodata)
    yld_arr = np.where(valid, yields, np.nan)

    # 2) Map user crop → IPCC crop key
    crop_table = _get_crop_naming_index_table()
    res_table = _get_crop_residue_ratio_table()
    ag_table = _get_crop_ag_residue_table()

    ipcc_crop = (
        crop_table
        .filter(pl.col("Crops") == crop)
        .select("IPCC_Crop")
        .to_series()
        .item()
    )

    # 3) Pull core residue parameters
    res_row = res_table.filter(pl.col("Crop") == ipcc_crop)
    dry      = float(res_row["DRY"].to_list()[0])
    dry_C_content = dry * C_Content
    
    # Seeing if there's an RS value
    try:
        RS       = float(res_row["RS"].to_list()[0])
    except (ValueError, TypeError):
        RS = 0
    
    # Looks for 
    try:
        R_AG     = float(res_row["R_AG"].to_list()[0])
    except (ValueError, TypeError):
        R_T      = float(res_row["R_T"].to_list()[0])

    # 4) See if regression parameters exist
    ag_row = ag_table.filter(pl.col("Crop") == crop)
    if ag_row.height > 0:
        slope     = float(ag_row["Slope"].to_list()[0])
        intercept = float(ag_row["Intercept"].to_list()[0])
        branch = "regression"
    elif RS > 0:
        branch = "ratio"
    else:
        branch = "total"

    # 5) Compute according to the chosen branch
    if branch == "regression":
        ABG = slope * yld_arr + intercept
        BG  = RS    * ABG
        Res = (ABG + BG) * dry_C_content

    elif branch == "ratio":
        ABG = R_AG * yld_arr
        BG  = RS   * ABG
        Res = (ABG + BG) * dry_C_content

    else:  # total‐residues branch
        Res = R_T * yld_arr * dry_C_content
        ABG = np.full_like(Res, np.nan, dtype="float32")
        BG  = np.full_like(Res, np.nan, dtype="float32")

    # restore nodata
    Res[~valid] = np.nan  # Assigns nan where there's no valid data
    if branch in ("regression", "ratio"):
        ABG[~valid] = np.nan
        BG [~valid] = np.nan

    # 6) Write out
    if branch in ("regression", "ratio"):
        # 3-band: Res, ABG, BG
        meta.update(count=3, dtype="float32", nodata=np.nan)
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(Res,  1)
            dst.write(ABG,  2)
            dst.write(BG,   3)
        print(f"[{branch}] → wrote 3-band raster: Res, ABG, BG")
    else:
        # 1-band: Res only
        meta.update(count=1, dtype="float32", nodata=np.nan)
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(Res, 1)
        print(f"[{branch}] → wrote 1-band raster: Res only")


def create_plant_cover_monthly_curve(
    crop: str,
    climate: str,
    *,
    crop_table: Optional[pl.DataFrame] = None,
):
    # Check the crops exists
    crop_table = _resolve_crop_coefficient_table(crop_table)

    if crop not in crop_table['Crop'].unique():
        raise ValueError(f"Crop '{crop}' not found in K_Crops table.")

    if climate not in crop_table["Climate_Zone"]:
        raise ValueError(f"Climate zone '{climate}' not found in K_Crops table.")

    # Retrieve plant cover data for the specified crop
    pc_starts = crop_table.filter((pl.col('Crop') == crop) & (pl.col('Climate_Zone') == climate)).select('SCP_Starts').item()
    pc_ends = crop_table.filter((pl.col('Crop') == crop) & (pl.col('Climate_Zone') == climate)).select('SCP_End').item()

    # Create a DataFrame for the plant cover curve
    plant_cover_array = pl.DataFrame(
        {
            "Month": list(range(1, 13)),
            "Plant_Cover": [0] * 12
        }
    )

    # Fill the plant cover curve based on start and end dates
    plant_cover_array = plant_cover_array.with_columns(
        pl.when((pl.col('Month')>=pc_starts), (pl.col('Month')<=pc_ends)).then(1).otherwise(0).alias('Plant_Cover')
    )

    return plant_cover_array


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


def convert_K2C_raster(kelvin_raster: str, output_path):
    # 1. Load without automatic masking so we can see the raw nodata tag
    t_K = rxr.open_rasterio(kelvin_raster, masked=False)

    # 2. Extract the nodata value from the file metadata
    nodata_val = t_K.rio.nodata

    # 3. Build a mask of valid (non-nodata) pixels
    valid_mask = t_K != nodata_val

    # 4. Subtract 273.15 only on valid data, leave others as nodata_val
    t_C = (t_K.where(valid_mask) - 273.15).astype("float32")

    # 5. Reapply the nodata tag so NaNs get written as your original nodata
    t_C = t_C.rio.write_nodata(nodata_val)

    # 6. Update metadata
    t_C.attrs["units"]       = "°C"
    t_C.attrs["description"] = "Monthly mean temperature in Celsius"

    # 7. (Optional) rename the band dimension for clarity
    if "band" in t_C.dims:
        t_C = t_C.rename({"band": "month"})

    # 8. Write out the new GeoTIFF, preserving CRS & transform
    t_C.rio.to_raster(output_path)


def create_plant_cover_monthly_raster(
    crop: str,
    save_path: str,
    climate_raster_path: str = uhth_climates_fp,
    output_nodata: int = 255,
    *,
    crop_table: Optional[pl.DataFrame] = None,
    climate_zone_lookup: Optional[Mapping[int, str]] = None,
):
    """
    Build a monthly (12×y×x) plant-cover mask from a climate-ID GeoTIFF.
    - crop: crop name for phenology lookup
    - climate_raster_path: path to a 1-band climate-ID TIFF (IDs 1–12)
    - save_path: if provided, writes a 12-band GeoTIFF
    - output_nodata: integer nodata code for the output mask
    """
    # 1. Load the climate raster (raw values, no masking)
    da = rxr.open_rasterio(climate_raster_path, masked=False)
    # If band dim exists, drop it
    if "band" in da.dims and da.sizes["band"] == 1:
        da = da.squeeze("band", drop=True)

    # 2. Get the raw ID grid and its spatial coords
    ids = da.values        # 2D array (y, x)
    y = da.coords["y"]
    x = da.coords["x"]

    # 3. Prepare an output array filled with nodata
    n_months = 12
    mask = np.full((n_months, y.size, x.size),
                   fill_value=output_nodata,
                   dtype=np.uint8)

    crop_table = _resolve_crop_coefficient_table(crop_table)
    climate_lookup = _resolve_climate_lookup(climate_zone_lookup)

    # 4. Loop over each unique climate ID
    # Pull nodata from rioxarray metadata, falling back to legacy attrs.
    nodata = da.rio.nodata
    if nodata is None:
        nodata = da.attrs.get("nodata")
    # Start with all pixels valid, then progressively filter out invalid ones.
    valid_mask = np.ones(ids.shape, dtype=bool)
    # Drop NaN values for floating-point rasters.
    if np.issubdtype(ids.dtype, np.floating):
        valid_mask &= ~np.isnan(ids)
    if nodata is not None:
        # Skip nodata comparison if nodata itself is NaN or non-numeric.
        try:
            nodata_is_nan = np.isnan(nodata)
        except TypeError:
            nodata_is_nan = False
        # Exclude explicit nodata values for integer or float rasters.
        if not nodata_is_nan:
            valid_mask &= ids != nodata
    unique_ids = np.unique(ids[valid_mask]).astype(int)
    for cid in unique_ids:
        group = climate_lookup.get(cid)
        if group is None:
            continue
        
        # Get the 12-month vector (0/1) from your existing function
        pc_df = create_plant_cover_monthly_curve(crop, group, crop_table=crop_table)
        pc_vec = np.array(pc_df.select("Plant_Cover").to_series())  # shape (12,)

        # Assign that vector to all pixels where ids==cid
        rows, cols = np.where(ids == cid)
        mask[:, rows, cols] = pc_vec[:, None]

    # 5. Wrap into an xarray.DataArray with spatial metadata
    da_mask = xr.DataArray(
        mask,
        dims=("month", "y", "x"),
        coords={
            "month": np.arange(1, n_months+1),
            "y": y,
            "x": x
        },
        name=f"{crop}_pc_mask"
    )
    # 6. Write CRS, transform, and nodata
    da_mask = da_mask.rio.write_crs(da.rio.crs)
    da_mask = da_mask.rio.write_transform(da.rio.transform())
    da_mask = da_mask.rio.write_nodata(output_nodata)
    da_mask = da_mask.rio.set_spatial_dims(x_dim="x", y_dim="y")

    # 7. Save
    da_mask.rio.to_raster(save_path)

    print(f"Plant cover raster saved to {save_path}")

# -----------------------------------------------------------------------------
# Residues
# -----------------------------------------------------------------------------
def compute_residue_raster(
    k_curve_df: pl.DataFrame,
    plant_residue: xr.DataArray,
    save_path: Optional[str] = None
) -> xr.DataArray:
    """
    Allocate annual plant residue raster across months based on a K-curve,
    and optionally save as a 12-band GeoTIFF.

    Parameters:
      - k_curve_df: Polars DataFrame with columns ['Month','K'] for 12 months.
      - plant_residue: xarray DataArray of total annual plant residue (dims 'y','x').
      - save_path: optional file path to write out the 12-band raster.

    Returns:
      - xarray DataArray of monthly residue (dims 'time','y','x').
    """
    # Sort and extract K values
    sorted_df = k_curve_df.sort('Month')
    k_vals = np.array(sorted_df['K'].to_list(), dtype=float)
    total_k = k_vals.sum()
    if total_k <= 0:
        raise ValueError("Sum of K values must be positive")
    
    # Monthly fractions
    fractions = k_vals / total_k
    
    # Prepare output array
    n_months = 12
    y_size = plant_residue.sizes['y']
    x_size = plant_residue.sizes['x']
    arr = np.empty((n_months, y_size, x_size), dtype=float)
    for i, frac in enumerate(fractions):
        arr[i, :, :] = frac * plant_residue.values
    
    # Build DataArray, inheriting spatial metadata
    residue_da = xr.DataArray(
        arr,
        dims=('time','y','x'),
        coords={
            'time': sorted_df['Month'].to_list(),
            'y': plant_residue.y,
            'x': plant_residue.x
        },
        name='residue'
    )
    # Copy Geo metadata from plant_residue
    residue_da = residue_da.rio.write_crs(plant_residue.rio.crs)
    residue_da = residue_da.rio.write_transform(plant_residue.rio.transform())

    # Optionally save
    if save_path:
        write_multiband_tif(residue_da, save_path, plant_residue)

    return residue_da

def _distribute_residue_monthly(
    crop: str,
    crop_type: str,
    climate_ids: np.ndarray,
    residue: np.ndarray,
    *,
    output_nodata: float,
    climate_zone_lookup: Optional[Mapping[int, str]],
    crop_coeff_table: Optional[pl.DataFrame] = None,
    climate_nodata: Optional[float],
) -> np.ndarray:
    """Return a (12, y, x) monthly residue cube using climate-specific Kc curves."""

    if climate_ids.shape != residue.shape:
        raise ValueError(
            "Climate-ID raster and residue raster must share the same shape; "
            f"got {climate_ids.shape} and {residue.shape}."
        )

    climate_lookup = _resolve_climate_lookup(climate_zone_lookup)
    crop_coeff_table = _resolve_crop_coefficient_table(crop_coeff_table)
    abs_day_table = _get_absolute_day_table()

    out = np.full((12, *residue.shape), fill_value=output_nodata, dtype="float32")

    if climate_nodata is None or np.isnan(climate_nodata):
        valid_mask = ~np.isnan(climate_ids)
    else:
        valid_mask = climate_ids != climate_nodata

    if not np.any(valid_mask):
        return out

    valid_ids = climate_ids[valid_mask].astype(int)
    unique_clim_ids = np.unique(valid_ids)

    for clim_id in unique_clim_ids:
        clim_group = climate_lookup.get(clim_id)
        if clim_group is None:
            continue

        # Get the needed data row
        crop_clim_data = crop_coeff_table.filter((pl.col("Crop")==crop) & (pl.col("Climate_Zone")==clim_group))
        plant_date = crop_clim_data.select("Planting_Greenup_Date").item()
        pd_abs_day = (
            abs_day_table
            .filter(pl.col("Date") == plant_date)
            .select("Day_Num")
            .to_series()
            .item()
        )

        # Calculates harvest date
        cycle_days = crop_clim_data.select(pl.sum_horizontal("Initial_days", "Dev_days", "Mid_days", "Late_days").alias("total_cycle_days")).item()
        end_abs_day = ((pd_abs_day + cycle_days - 1) % 365) + 1  # back to 1..365

        # Gets harvest month
        harvest_month = abs_day_table.filter(pl.col("Day_Num") == end_abs_day).select(pl.col("Month")).item()
        hm_0index = harvest_month - 1

        # Create a fraction output
        month_frac = np.zeros(shape=(12,), dtype="float32")
        # Now 2 routes
        if crop_type == "annual":
            res_months = np.arange(hm_0index - 3, hm_0index, 1) % 12  # %12 is to standardize into 0-11 months
            
            # now assigining fractions into months:
            month_frac[hm_0index] = 0.5
            month_frac[res_months] = 0.5/3
        else:  # Permanent crops
            res_months = np.arange(hm_0index - 4, hm_0index, 1) % 12  # %12 is to standardize into 0-11 months

            # now assigining fractions into months:
            month_frac[hm_0index] = 0.7
            month_frac[res_months] = 0.3/4

        rows, cols = np.where((climate_ids == clim_id) & valid_mask)
        if rows.size == 0:
            continue

        pr_vals = residue[rows, cols]
        out[:, rows, cols] = month_frac[:, None] * pr_vals[None, :]

    return out


def compute_monthly_residue_raster(
    crop: str,
    crop_type: str,
    climate_raster_path: str,
    plant_residue: xr.DataArray,
    save_path: str,
    output_nodata: float = np.nan,
    *,
    climate_zone_lookup: Optional[Mapping[int, str]] = None,
    crop_table: Optional[pl.DataFrame] = None,
):
    """
    Allocate annual plant_residue into monthly residues per pixel,
    based on crop-specific Kc curves per climate group.

    Parameters
    ----------
    crop : str
        Crop name, passed to create_KC_Curve().
    climate_raster_path : str
        Path to a 1-band climate-ID GeoTIFF (values 1–12).
    plant_residue : xr.DataArray
        2D DataArray (y, x) of annual residue (t C/ha).
    save_path : str, optional
        If provided, writes out a 12-band GeoTIFF of monthly residues.
    output_nodata : float, default np.nan
        Value to use for pixels with no valid climate ID or missing residue.

    Returns
    -------
    xr.DataArray
        3D DataArray dims=('month','y','x') of monthly residue (t C/ha).
    """
    clim = rxr.open_rasterio(climate_raster_path, masked=False)
    if "band" in clim.dims:
        clim = clim.isel(band=0)
    ids = clim.values
    climate_nodata = clim.rio.nodata

    if "band" in plant_residue.dims:
        plant_residue = plant_residue.isel(band=0)
    if plant_residue.ndim != 2:
        raise ValueError(
            "Expected plant_residue to be 2D after squeezing, "
            f"got {plant_residue.shape}"
        )

    monthly = _distribute_residue_monthly(
        crop,
        crop_type,
        ids,
        plant_residue.values,
        output_nodata=output_nodata,
        climate_zone_lookup=climate_zone_lookup,
        climate_nodata=climate_nodata,
    )

    da = xr.DataArray(
        monthly,
        dims=("month", "y", "x"),
        coords={
            "month": np.arange(1, 13),
            "y": clim.coords["y"],
            "x": clim.coords["x"],
        },
        name=f"{crop}_residue_monthly",
    ).astype("float32")

    da = da.rio.write_crs(clim.rio.crs)
    da = da.rio.write_transform(clim.rio.transform())
    da = da.rio.write_nodata(output_nodata)
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    da.rio.to_raster(save_path)


def compute_monthly_residue_raster_fromAnnualRaster(
    crop: str,
    crop_type: str,
    plant_residue: str,
    save_path: str,
    climate_raster_path: str = uhth_climates_fp,
    output_nodata: float = np.nan,
    *,
    climate_zone_lookup: Optional[Mapping[int, str]] = None,
    crop_table: Optional[pl.DataFrame] = None,
):
    """
    Allocate annual plant_residue into monthly residues per pixel, based on crop-specific Kc curves per climate group.

    Parameters
    ----------
    crop : str
        Crop name, passed to create_KC_Curve().
    climate_raster_path : str
        Path to a 1-band climate-ID GeoTIFF (values 1–12).
    plant_residue : str
        Path to annual residue GeoTIFF raster
    save_path : str, optional
        If provided, writes out a 12-band GeoTIFF of monthly residues.
    output_nodata : float, default np.nan
        Value to use for pixels with no valid climate ID or missing residue.

    Returns
    -------
    xr.DataArray
        3D DataArray dims=('month','y','x') of monthly residue (t C/ha).
    """
    clim = rxr.open_rasterio(climate_raster_path, masked=False)
    if "band" in clim.dims:
        clim = clim.isel(band=0)
    ids = clim.values
    climate_nodata = clim.rio.nodata

    pr_da = rxr.open_rasterio(plant_residue, masked=True)
    if "band" in pr_da.dims:
        pr_da = pr_da.isel(band=0)
    if pr_da.ndim != 2:
        raise ValueError(
            "Expected plant_residue to be 2D after squeezing, "
            f"got {pr_da.shape}"
        )

    monthly = _distribute_residue_monthly(
        crop,
        crop_type,
        ids,
        pr_da.values,
        output_nodata=output_nodata,
        climate_zone_lookup=climate_zone_lookup,
        climate_nodata=climate_nodata,
    )

    da = xr.DataArray(
        monthly,
        dims=("month", "y", "x"),
        coords={
            "month": np.arange(1, 13),
            "y": pr_da.coords["y"],
            "x": pr_da.coords["x"],
        },
        name=f"{crop}_residue_monthly",
    ).astype("float32")

    da = da.rio.write_crs(clim.rio.crs)
    da = da.rio.write_transform(clim.rio.transform())
    da = da.rio.write_nodata(output_nodata)
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    da.rio.to_raster(save_path)


# -----------------------------------------------------------------------------
# Irrigation
# -----------------------------------------------------------------------------
def calculate_irrigation_fromArray(rain, evap):
    """Calculate monthly theoretical irrigation needs for a crop based on rain and evapotranspiration needs for said crop.

    Args:
        rain (array-like): Monthly or daily precipitation (mm/month).
        evap (array-like): Corresponding evapotranspiration demand (mm/month).

    Returns:
        irr (np.array): Irrigation required (mm/month), where ET > precipitation; 0 elsewhere.
    """

    # Transform data into array if needed
    rain = np.asarray(rain, dtype=float)
    evap = np.asarray(evap, dtype=float)

    # Creates a new empty array
    irr = np.zeros_like(rain, dtype=float)

    # See where it needs irrigation
    irr_needed = evap > rain

    # Fills the irrigation array
    irr = np.where(irr_needed, evap - rain, 0)

    return irr


def calculate_irrigation_fromTif(rain_fp, evap_fp, out_path: str):
    """Calculate monthly theoretical irrigation needs for a crop based on rain and evapotranspiration needs for said crop.

    Args:
        rain (array-like): Monthly or daily precipitation (mm/month).
        evap (array-like): Corresponding evapotranspiration demand (mm/month).

    Returns:
        irr (np.array): Irrigation required (mm/month), where ET > precipitation; 0 elsewhere.
    """
    with rasterio.open(rain_fp) as src_rain:
        rain = src_rain.read().astype("float32")
        rain_src = src_rain.crs

    with rasterio.open(evap_fp) as src_evap:
        evap = src_evap.read().astype("float32")
        evap_src = src_evap.crs
        evap_profile = src_evap.profile

    # Check if crs is the same
    if rain_src != evap_src:
        raise ValueError("Rasters have different crs. Please align before.")

    # Checks if size are the same
    if src_rain.shape != src_evap.shape:
        raise ValueError("Rasters different shape")

    # See where it needs irrigation
    irr_needed = evap > rain

    # Fills the irrigation array
    irr = np.where(irr_needed, evap - rain, np.nan)

    # Saves the result
    evap_profile.update(dtype='float32', count=12, nodata=np.nan)
    with rasterio.open(out_path, "w", **evap_profile) as dst:
        dst.write(irr.astype("float32"))
    print(f"Irrigation raster saved to {out_path}")


############################################
#### CROP DATA PREPARATION STREAMLINING ####
############################################
def prepare_crop_data(
    crop_name: str,
    crop_type: str,
    crop_practice_string: str,
    lu_data_path: str,
    spam_crop_raster: str,
    output_data_folder: str,
    irr_yield_scaling: str,
    spam_all_fp: str,
    spam_irr_fp: str,
    spam_rf_fp: str,
    all_new_files: bool = False,
):
    # Check if crop_type is valid
    if crop_type not in crop_types:
        raise ValueError(f"Crop type {crop_type} not valid. Choose between {crop_types}")

    # Output saving string bases
    output_base = Path(output_data_folder)
    output_crop_based = output_base / crop_name
    output_practice_based = output_base / f"{crop_name}_{crop_practice_string}"

    # Step 0 - rasterize input path
    lu_bin_output = output_practice_based.parent / f"{output_practice_based.name}_lu.tif"
    if all_new_files or not lu_bin_output.exists():
        print("Creating lu raster...")
        lu_array = _binarize_raster_pipeline(lu_data_path, str(lu_bin_output))
    else:
        print("Land use binary raster already exist. Skipping...")
        lu_array = rxr.open_rasterio(lu_bin_output, masked=False).squeeze()

    # Step 1 - Prepare PET and irrigation
    # Step 1.1 - PET
    pet_monthly_output_path = output_crop_based.parent / f"{output_crop_based.name}_pet_monthly.tif"
    if all_new_files or not pet_monthly_output_path.exists():
        print("Creating PET raster...")
        monthly_pet = calculate_crop_based_PET_raster_vPipeline(
            crop_name=crop_name,
            landuse_array=lu_array,
            output_monthly_path=str(pet_monthly_output_path)
        )
    else:
        print("PET raster already exists. Skipping...")
        monthly_pet = rxr.open_rasterio(pet_monthly_output_path, masked=True).values

    # Step 1.2 - Irrigation
    irr_monthly_output_path = output_crop_based.parent / f"{output_crop_based.name}_irr_monthly.tif"
    if all_new_files or not irr_monthly_output_path.exists():
        print("Creating irrigation raster...")
        irr = calculate_irrigation_vPipeline(
            evap=monthly_pet,
            output_path=str(irr_monthly_output_path)
        )
    else:
        print("Irrigation raster already exists — skipping computation.")

    # Step 2 - Calculate yields
    # preparing fao yield shapefile
    crop_names_table = _get_crop_naming_index_table()
    fao_crop_name = crop_names_table.filter(pl.col("Crops")== crop_name).select(pl.col("FAO_Crop")).item()
    print(f"Creating {fao_crop_name} helper shapefile...")
    fao_yield_shp = create_crop_yield_shapefile(fao_crop_name)

    # Create irrigation adjusted yields
    yield_output_path = output_practice_based.parent / f"{output_practice_based.name}_yield.tif"
    if all_new_files or not yield_output_path.exists():
        print("Creating yield raster...")
        _ = create_crop_yield_raster_with_irrigation_scaling_pipeline(
            croplu_grid_raster=str(lu_bin_output),
            fao_crop_shp=fao_yield_shp,
            spam_crop_raster=spam_crop_raster,
            output_rst_path=str(yield_output_path),
            irr_yield_scaling=irr_yield_scaling,
            all_fp = spam_all_fp,
            irr_fp = spam_irr_fp,
            rf_fp= spam_rf_fp
        )
    else:
        print("Yields raster already exists — skipping computation.")

    # Step 3 - Create plant cover raster
    plantcover_output_path = output_crop_based.parent / f"{output_crop_based.name}_pc_monthly.tif"
    if all_new_files or not plantcover_output_path.exists():
        print("Creating plant cover raster...")
        create_plant_cover_monthly_raster(crop_name, str(plantcover_output_path))
    else:
        print("Plant Cover raster already exists — skipping computation.")

    # Step 4 - Create plant residue raster
    plant_residue_output_path = output_practice_based.parent / f"{output_practice_based.name}_residues_monthly.tif"
    if all_new_files or not plant_residue_output_path.exists():
        print("Creating plant residue raster...")
        create_monthly_residue_vPipeline(
            crop_name,
            crop_type,
            yield_raster_path=str(yield_output_path),
            output_path=str(plant_residue_output_path),
            write_output=True
        )
    else:
        print("Plant Residues raster already exists — skipping computation.")

    print(f"All data created for {crop_name}, {crop_practice_string}!!!")

def prepare_crop_data_irrigation_plantcover(
    crop_name: str,
    crop_type: str,
    crop_practice_string: str,
    lu_data_path: str,
    output_data_folder: str,
    all_new_files: bool = False,
):
    # Check if crop_type is valid
    if crop_type not in crop_types:
        raise ValueError(f"Crop type {crop_type} not valid. Choose between {crop_types}")

    # Output saving string bases
    output_base = Path(output_data_folder)
    output_crop_based = output_base / crop_name
    output_practice_based = output_base / f"{crop_name}_{crop_practice_string}"

    # Step 0 - Opening input path
    # Load land-use raster so we can enforce single-band inputs before squeezing.
    lu_raster = rxr.open_rasterio(lu_data_path, masked=False)
    # Determine how many bands are present; default to 1 if the dimension is missing.
    band_count = lu_raster.sizes.get("band", 1)
    if band_count != 1:
        raise ValueError(
            f"Land-use raster '{lu_data_path}' has {band_count} bands; "
            "multiband land-use rasters are not supported. Please provide a single-band raster."
        )
    # Explicitly select the first band to guarantee a 2D (H, W) array for downstream logic.
    lu_array = lu_raster.isel(band=0).squeeze()
    if lu_array.ndim != 2:
        # Fail fast if the land-use raster still isn't 2D after band selection.
        raise ValueError(
            f"Land-use raster '{lu_data_path}' must be a 2D array shaped (H, W); "
            f"got shape {lu_array.shape}."
        )

    # Step 1 - Prepare PET and irrigation
    # Step 1.1 - PET
    pet_monthly_output_path = output_practice_based.parent / f"{output_practice_based.name}_pet_monthly.tif"
    if all_new_files or not pet_monthly_output_path.exists():
        print("Creating PET raster...")
        monthly_pet = calculate_crop_based_PET_raster_vPipeline(
            crop_name=crop_name,
            landuse_array=lu_array,
            output_monthly_path=str(pet_monthly_output_path)
        )
    else:
        print("PET raster already exists. Skipping...")
        monthly_pet = rxr.open_rasterio(pet_monthly_output_path, masked=True).values

    # Step 1.2 - Irrigation
    irr_monthly_output_path = None
    if "irr" in crop_practice_string:
        irr_monthly_output_path = output_practice_based.parent / f"{output_practice_based.name}_irr_monthly.tif"
        if all_new_files or not irr_monthly_output_path.exists():
            print("Creating irrigation raster...")
            irr = calculate_irrigation_vPipeline(
                evap=monthly_pet,
                output_path=str(irr_monthly_output_path)
            )
        else:
            print("Irrigation raster already exists — skipping computation.")
    else:
        print("Irrigation not needed. Skipping...")
    

    # Step 3 - Create plant cover raster
    plantcover_output_path = output_practice_based.parent / f"{output_practice_based.name}_pc_monthly.tif"
    if all_new_files or not plantcover_output_path.exists():
        print("Creating plant cover raster...")
        create_plant_cover_monthly_raster(crop_name, str(plantcover_output_path))
    else:
        print("Plant Cover raster already exists — skipping computation.")

    print(f"Irrigation, PET, and plant cover rasters created for {crop_name}, {crop_practice_string}!!!")

    # Return all lu path associated with the scenario run
    return lu_data_path, pet_monthly_output_path, irr_monthly_output_path, plantcover_output_path


def calculate_monthly_residues_array(
    lu_fp: str,
    crop_name: str,
    crop_type: str,
    spam_crop_raster: str,
    irr_yield_scaling: str,
    spam_all_fp: str,
    spam_irr_fp: str,
    spam_rf_fp: str,
    random_runs: int,
    print_outputs: bool = False,
    spam_outlier_strategy: str = "spam_sd",
    spam_outlier_percentile: Tuple[float, float] = (1.0, 99.0),
    spam_outlier_k: float = 2.0,
    ylds_src: str = "GAEZ"
):
    # print("    Calculating stochastic residue array...")

    # Step 1 - Prepare fao yield shapefile
    crop_names_table = _get_crop_naming_index_table()
    fao_crop_name = crop_names_table.filter(pl.col("Crops") == crop_name).select(pl.col("FAO_Crop")).item()
    
    # print(f"Creating {fao_crop_name} helper shapefile...")
    fao_yield_shp = create_crop_yield_shapefile(fao_crop_name)

    yield_result = calculate_crop_yield_array_with_irrigation_scaling(
        croplu_grid_raster_fp=   lu_fp,
        fao_crop_shp=fao_yield_shp,
        spam_crop_raster=spam_crop_raster,
        irr_yield_scaling=irr_yield_scaling,
        all_fp=spam_all_fp,
        irr_fp=spam_irr_fp,
        rf_fp=spam_rf_fp,
        random_runs=random_runs,
        print_outputs= print_outputs,
        spam_outlier_strategy=spam_outlier_strategy,
        spam_outlier_percentile=spam_outlier_percentile,
        spam_outlier_k=spam_outlier_k,
    )

    # Step 4 - Create plant residue raster
    plant_residues = create_monthly_residue_vPipeline(
        crop_name,
        crop_type,
        yield_array=yield_result.averaged_result,
        write_output=False,
        return_array=True
    )

    return plant_residues

def prepare_crop_scenarios(csv_filepath: str, override_params: dict | None = None):
    # Load scenarios
    csv = pl.read_csv(csv_filepath)
    scenarios = csv.to_dicts()

    # Run scenarions
    for scenario in scenarios:
        scenario = scenario.copy()

        # Apply overrides if given
        if override_params is not None:
            scenario.update(override_params)

        print(f"Preparing data for {scenario['crop_name']}, {scenario['crop_practice_string']}")
        prepare_crop_data(**scenario)
        print(f"Next!\n")


def prepare_crop_scenarios_PET_PlantCover_only(csv_filepath: str, override_params: dict | None = None):
    # Load scenarios
    scneraios_df = pl.read_csv(csv_filepath)
    all_scenarios = scneraios_df.to_dicts()

    # Store unique land use scenarios, to not duplicate unneeded, heavy files
    unique_scenarios: list[dict] = []

    # For each crop, remember which LU hashes we have already seen
    seen_hashes_by_crop: Dict[str, Set[str]] = defaultdict(set)

    # Cache raster hashes by path to avoid re-reading the same file
    hash_cache: Dict[str, str] = {}
    
    print("Creating a unique list of rasters with different land use maps...\n")

    #1.- Iterate row-by-row (as dicts)
    for row in scneraios_df.iter_rows(named=True):
        # Needed information
        crop_name = row["crop_name"]
        crop_lu_path = row["lu_data_path"]

        # Get or compute hash for this LU raster
        if crop_lu_path in hash_cache:
            lu_hash = hash_cache[crop_lu_path]
        else:
            lu_hash = _hash_raster(crop_lu_path)
            hash_cache[crop_lu_path] = lu_hash

        # Check if we've already seen this LU for this crop. If yes, skip all that follows
        if lu_hash in seen_hashes_by_crop[crop_name]:
            continue

        ### This only happens if it has not been seen ###
        # First time seeing this LU for this crop → keep the scenario
        seen_hashes_by_crop[crop_name].add(lu_hash)
        unique_scenarios.append(row)

    #2.- Add back irrigation scenarios, as the irrigation rasters always need to be computed
    irrigation_df = scneraios_df.filter(pl.col("crop_practice_string").str.contains("irr_"))

    # Turn unique_scenarios (list-of-dicts) into a Polars DF to compare
    unique_df = pl.DataFrame(unique_scenarios) if unique_scenarios else scneraios_df.clear()

    # check which irrigation scenarios are missing
    join_keys = ["crop_name", "crop_practice_string"]
    extra_irrigation = irrigation_df.join(unique_df, on=join_keys, how="anti")  # anti keeps only the results missing

    for row in extra_irrigation.iter_rows(named=True):
        crop_name = row["crop_name"]
        lu_path = row["lu_data_path"]

        # get or compute the LU hash
        if lu_path in hash_cache:
            lu_hash = hash_cache[lu_path]
        else:
            lu_hash = _hash_raster(lu_path)
            hash_cache[lu_path] = lu_hash

        # if this LU hash is already seen for this crop, skip it
        if lu_hash in seen_hashes_by_crop[crop_name]:
            continue

        # otherwise, this is a new LU for this crop → keep it
        seen_hashes_by_crop[crop_name].add(lu_hash)
        unique_scenarios.append(row)

    # Append these missing irrigation scenarios to the list
    unique_scenarios.extend(extra_irrigation.to_dicts())

    # 3.- Run scenarions
    # Creates an emtpy list
    all_results = []
    for scenario in unique_scenarios:
        scenario = scenario.copy()

        # Apply overrides if given
        if override_params is not None:
            scenario.update(override_params)

        print(f"Preparing irrigation and plant cover data for {scenario['crop_name']}, {scenario['crop_practice_string']}")
        
        # Prepare crop data scenario
        lu_path, pet_path, irr_path, plantcover_path = prepare_crop_data_irrigation_plantcover(**scenario)

        # Build one combined record
        result_entry = {
            "crop_name": scenario["crop_name"],
            "crop_type": scenario["crop_type"],
            "crop_practice_string": scenario["crop_practice_string"],
            "lu_path": lu_path,
            "pet_path": pet_path,
            "irr_path": irr_path,
            "plantcover_path": plantcover_path,
        }
        # Append it
        all_results.append(result_entry)
        
        # Continues
        print(f"Next!\n")

    all_results_df = pd.DataFrame(all_results)

    return all_results_df


#### PIPELINE SUPPORTING FUNCTIONS
def _hash_raster(path: str) -> str:
    """
    Compute a hash for the raster at 'path' based on its data. Used to detect identical LU rasters.
    """
    with rasterio.open(path) as src:
        # Read band 1 as masked array
        arr = src.read(1, masked=True)

        # pick correct nodata fill
        fill_val = src.nodata if src.nodata is not None else 255

        # Normalize the data: fill mask with a stable value
        data = np.ma.filled(arr, fill_value=fill_val)

        # Create hash object
        h = hashlib.sha1()
        # Include shape so different sized rasters never collide trivially
        h.update(str(data.shape).encode("utf-8"))
        # Include the transform as part of identity (optional but sensible)
        h.update(str(src.transform).encode("utf-8"))
        # Hash the actual pixel bytes
        h.update(data.tobytes())

        return h.hexdigest()
    

def _binarize_raster_pipeline(
    src_path: str,
    dst_path: str,
    nodata_value: int = 255,
    band: int = 1
):
    """
    Create a 1-band raster where:
      - pixels with any valid input value → 1
      - pixels that are src.nodata or NaN → nodata_value
    The output's nodata is set to src.nodata (if defined) or to nodata_value.
    """
    with rasterio.open(src_path) as src:
        data       = src.read(band)
        src_nodata = src.nodata

        # decide what our output nodata will be
        out_nodata = src_nodata if src_nodata is not None else nodata_value

        # start with everything set to nodata. Creates an array of shape given by data.shape and fills it with out_nodata values
        mask = np.full(shape=data.shape, fill_value=out_nodata, dtype="uint8")

        # build a boolean of "valid" pixels
        if src_nodata is not None:
            valid = (data != src_nodata)  # valid pixels are those that are not equal to the src.nodata value
        else:
            valid = ~np.isnan(data)  # valid pixels are those that are not NaN

        # assign 1 to all valid pixels
        mask[valid] = 1

        # update profile
        profile = src.profile.copy()
        profile.update(
            dtype="uint8",
            count=1,
            nodata=out_nodata
        )

    # write out
    with rasterio.open(dst_path, "w", **profile) as dst:
        dst.write(mask, 1)

    return mask

def calculate_irrigation_vPipeline(evap: np.ndarray, output_path: str, rain_fp = rain_monthly_fp):
    """Calculate monthly theoretical irrigation needs for a crop based on rain and evapotranspiration needs for said crop.

    Args:
        rain (array-like): Monthly or daily precipitation (mm/month).
        evap (array-like): Corresponding evapotranspiration demand (mm/month).

    Returns:
        irr (np.array): Irrigation required (mm/month), where ET > precipitation; 0 elsewhere.
    """

    def _get_rio_metadata(data: object) -> Tuple[Optional[CRS], Optional[Affine]]:
        if not isinstance(data, xr.DataArray):
            return None, None
        try:
            crs = data.rio.crs
        except Exception:
            return None, None
        if crs is None:
            return None, None
        try:
            transform = data.rio.transform()
        except Exception:
            transform = None
        return crs, transform

    # Transform data into array if needed
    with rasterio.open(rain_fp) as src:
        rain = src.read().astype("float32")
        rain_crs = src.crs
        rain_profile = src.profile

    rain = np.asarray(rain, dtype=float)

    evap_crs, evap_transform = _get_rio_metadata(evap)
    rain_transform = rain_profile.get("transform")
    rain_da = None

    shape_mismatch = np.asarray(evap).shape != rain.shape
    grid_mismatch = (
        evap_crs is not None
        and rain_crs is not None
        and (evap_crs != rain_crs or (evap_transform is not None and rain_transform is not None and evap_transform != rain_transform))
    )

    if shape_mismatch or grid_mismatch:
        if evap_crs is None or evap_transform is None:
            raise ValueError(
                "Evapotranspiration inputs do not align with the rainfall grid. "
                "Please provide evap inputs with matching shape and spatial alignment "
                "(CRS/transform) or supply a georeferenced xarray DataArray so it can be resampled."
            )
        if rain_da is None:
            rain_da = rxr.open_rasterio(rain_fp, masked=False)
        evap = evap.rio.reproject_match(rain_da)

    evap = np.asarray(evap, dtype=float)
    if evap.shape != rain.shape:
        raise ValueError(
            "Evapotranspiration inputs do not align with the rainfall grid after resampling. "
            "Please provide aligned inputs with matching shape and spatial metadata."
        )

    # Creates a new empty array
    irr = np.zeros_like(rain, dtype=float)

    # See where it needs irrigation
    irr_needed = evap > rain

    # Fills the irrigation array
    irr = np.where(irr_needed, evap - rain, np.nan)

    # Saves the result
    rain_profile.update(dtype='float32', count=12, nodata=np.nan)
    with rasterio.open(output_path, "w", **rain_profile) as dst:
        dst.write(irr.astype("float32"))
    print(f"Irrigation raster saved to {output_path}")

    return irr


def create_crop_yield_raster_with_irrigation_scaling_pipeline(
    croplu_grid_raster: str,
    fao_crop_shp: "gpd.GeoDataFrame",
    spam_crop_raster: str,
    output_rst_path: str,
    spam_band: int = 1,
    resampling_method: Resampling = Resampling.bilinear,
    irr_yield_scaling: Optional[str] = None,
    all_fp: Optional[str] = None,
    irr_fp: Optional[str] = None,
    rf_fp: Optional[str] = None,
    fao_avg_yield_name: str = "avg_yield",
    fao_yield_ratio_name: str = "yld_ratio",
    fao_sd_yield_name: str = "sd_yield",
    apply_ecoregion_fill: bool = True,
    random_runs: int = 1,
    print_outputs: bool = False,
    ylds_src = "GAEZ",
    spam_outlier_strategy: str = "spam_sd",
    spam_outlier_percentile: Tuple[float, float] = (1.0, 99.0),
    spam_outlier_k: float = 2.0,
    p1_window: int = 5,     # 5 pixels
    p1_stat: str = "median",
    # pass 2
    p2_enable: bool = True,
    p2_window: int = 3,     # 3 pixels
    p2_stat: str = "median",
    p2_fallback: str = "clip"  # "clip" or "neighborhood_fill"
) -> CropYieldRasterResult:
    """Pipeline wrapper around :func:`create_crop_yield_raster_withIrrigationPracticeScaling`."""

    config = CropYieldRasterConfig(
        fao_avg_yield_name=fao_avg_yield_name,
        fao_yield_ratio_name=fao_yield_ratio_name,
        fao_sd_yield_name=fao_sd_yield_name,
        irr_yield_scaling=irr_yield_scaling,
        all_fp=all_fp,
        irr_fp=irr_fp,
        rf_fp=rf_fp,
        spam_band=spam_band,
        resampling_method=resampling_method,
        apply_ecoregion_fill=apply_ecoregion_fill,
        random_runs=random_runs,
        print_outputs=print_outputs,
        ylds_src = ylds_src,
        spam_outlier_strategy = spam_outlier_strategy,
        spam_outlier_percentile = spam_outlier_percentile,
        spam_outlier_k = spam_outlier_k,
        p1_window =p1_window,
        p1_stat = p1_stat,
        # pass 2
        p2_enable = p2_enable,
        p2_window = p2_window,
        p2_stat = p2_stat,
        p2_fallback = p2_fallback,
    )
    return _create_crop_yield_raster_core(
        croplu_grid_raster,
        fao_crop_shp,
        spam_crop_raster,
        output_rst_path,
        config,
    )

def calculate_crop_yield_array_with_irrigation_scaling(
    croplu_grid_raster_fp: str,
    fao_crop_shp: "gpd.GeoDataFrame",
    spam_crop_raster: str,
    spam_band: int = 1,
    resampling_method: Resampling = Resampling.bilinear,
    irr_yield_scaling: Optional[str] = None,
    all_fp: Optional[str] = None,
    irr_fp: Optional[str] = None,
    rf_fp: Optional[str] = None,
    fao_avg_yield_name: str = "avg_yield",
    fao_yield_ratio_name: str = "yld_ratio",
    fao_sd_yield_name: str = "sd_yield",
    apply_ecoregion_fill: bool = True,
    random_runs: int = 1,
    print_outputs: bool = False,
    spam_outlier_strategy: str = "spam_sd",
    spam_outlier_percentile: Tuple[float, float] = (1.0, 99.0),
    spam_outlier_k: float = 2.0,
    ylds_src: str = "GAEZ",
    p1_window: int = 5,     # 5 pixels
    p1_stat: str = "median",
    # pass 2
    p2_enable: bool = True,
    p2_window: int = 3,     # 3 pixels
    p2_stat: str = "median",
    p2_fallback: str = "clip"  # "clip" or "neighborhood_fill"
) -> CropYieldRasterResult:
    """Pipeline wrapper around :func:`create_crop_yield_raster_withIrrigationPracticeScaling`."""

    config = CropYieldRasterConfig(
        fao_avg_yield_name=fao_avg_yield_name,
        fao_yield_ratio_name=fao_yield_ratio_name,
        fao_sd_yield_name=fao_sd_yield_name,
        irr_yield_scaling=irr_yield_scaling,
        all_fp=all_fp,
        irr_fp=irr_fp,
        rf_fp=rf_fp,
        spam_band=spam_band,
        resampling_method=resampling_method,
        apply_ecoregion_fill=apply_ecoregion_fill,
        random_runs=random_runs,
        write_output=False,
        return_array=True,
        print_outputs=print_outputs,
        spam_outlier_strategy = spam_outlier_strategy,
        spam_outlier_percentile = spam_outlier_percentile,
        spam_outlier_k = spam_outlier_k,
        ylds_src = ylds_src,
        p1_window =p1_window,
        p1_stat = p1_stat,
        # pass 2
        p2_enable = p2_enable,
        p2_window = p2_window,
        p2_stat = p2_stat,
        p2_fallback = p2_fallback,
    )
    return _create_crop_yield_raster_core(
        croplu_grid_raster= croplu_grid_raster_fp,
        fao_crop_shp = fao_crop_shp,
        spam_crop_raster = spam_crop_raster,
        output_rst_path=None,
        config=config,
    )

def create_monthly_residue_vPipeline(
    crop: str,
    crop_type: str,
    output_path: Optional[str] = None,
    yield_array: Optional[np.ndarray] = None,
    yield_raster_path: Optional[str] = None,
    output_nodata = np.nan,
    climate_raster_path: str = uhth_climates_fp,
    c_content: float = 0.40,
    *,
    climate_zone_lookup: Optional[Mapping[int, str]] = None,
    write_output: bool = False,
    return_array: bool = False,
):
    """
    Read a single‐band yield raster, choose one residue‐calculation path globally:
      1) Regression (Slope/Intercept) if available
      2) Else Ratio (R_AG, RS) if RS > 0
      3) Else Total‐residues (R_T, dry, C_Content)
    and write either:
      • a 3-band TIFF (Res, ABG, BG) for paths 1 & 2, or
      • a 1-band TIFF (Res only) for path 3.
    """

    # 1) Load the yield raster
    if yield_raster_path is not None:
        with rasterio.open(yield_raster_path) as src:
            shape   = src.shape
            yields  = src.read(1).astype("float32")
            nodata  = src.nodata
            src_crs = src.crs
            src_transform = src.transform
        
        valid = (yields != nodata)
    else:
        yields = yield_array
        valid = ~np.isnan(yields)

    # mask nodata → NaN
    
    yld_arr = np.where(valid, yields, np.nan)

    # 2) Map user crop → IPCC crop key
    crop_names_table = _get_crop_naming_index_table()
    res_table = _get_crop_residue_ratio_table()
    ag_table = _get_crop_ag_residue_table()

    ipcc_crop = (
        crop_names_table
        .filter(pl.col("Crops") == crop)
        .select("IPCC_Crop")
        .to_series()
        .item()
    )

    # 3) Pull core residue parameters
    res_row = res_table.filter(pl.col("Crop") == ipcc_crop)
    dry      = float(res_row["DRY"].to_list()[0])
    dry_C_content = dry * c_content
    
    # Seeing if there's an RS value
    try:
        RS       = float(res_row["RS"].to_list()[0])
    except (ValueError, TypeError):
        RS = 0
    
    # Looks for 
    try:
        R_AG     = float(res_row["R_AG"].to_list()[0])
    except (ValueError, TypeError):
        R_T      = float(res_row["R_T"].to_list()[0])

    # 4) See if regression parameters exist
    ag_row = ag_table.filter(pl.col("Crop") == crop)
    if ag_row.height > 0:
        slope     = float(ag_row["Slope"].to_list()[0])
        intercept = float(ag_row["Intercept"].to_list()[0])
        branch = "regression"
    elif RS > 0:
        branch = "ratio"
    else:
        branch = "total"

    # 5) Compute according to the chosen branch
    if branch == "regression":
        ABG = slope * yld_arr + intercept
        BG  = RS    * ABG
        Res = (ABG + BG) * dry_C_content

    elif branch == "ratio":
        ABG = R_AG * yld_arr
        BG  = RS   * ABG
        Res = (ABG + BG) * dry_C_content

    else:  # total‐residues branch
        Res = R_T * yld_arr * dry_C_content
        ABG = np.full_like(Res, np.nan, dtype="float32")
        BG  = np.full_like(Res, np.nan, dtype="float32")

    # restore nodata
    Res[~valid] = np.nan  # Assigns nan where there's no valid data
    if branch in ("regression", "ratio"):
        ABG[~valid] = np.nan
        BG [~valid] = np.nan

    ##############################
    ### Distributing per month ###
    ##############################

    clim = rxr.open_rasterio(climate_raster_path, masked=False)
    if "band" in clim.dims:
        clim = clim.isel(band=0)
    ids = clim.values
    climate_nodata = clim.rio.nodata

    monthly = _distribute_residue_monthly(
        crop,
        crop_type,
        ids,
        Res,
        output_nodata=output_nodata,
        climate_zone_lookup=climate_zone_lookup,
        climate_nodata=climate_nodata,
    )

    da = xr.DataArray(
        monthly,
        dims=("month", "y", "x"),
        coords={
            "month": np.arange(1, 13),
            "y":   clim.coords["y"],
            "x":   clim.coords["x"],
        },
        name=f"{crop}_residue_monthly",
    ).astype("float32")

    # 1️⃣ Tell rioxarray which dims are spatial
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    if write_output:
        # 2️⃣ Write CRS, transform, nodata in order
        da.rio.to_raster(
            output_path,
            driver="GTiff",
            crs=src_crs,
            transform=src_transform,
            dtype="float32",
            nodata=output_nodata,
        )

    if return_array:
        return da.values


##############################
#### FOREST CALCULATIONS #####
##############################
try:
    forest_litter_table = pl.read_excel(data_path("forest", "forest_residues_IPCC.xlsx"))
except FileNotFoundError:  # pragma: no cover - optional input tables
    forest_litter_table = pl.DataFrame(
        {
            "IPCC Climate": ["Temperate"],
            "BD_mean": [1.0],
            "NE_mean": [1.0],
            "BD_TP": [20.0],
            "NE_TP": [20.0],
        }
    )
def get_forest_litter_rate(da_fp: str, forest_type: str, weather_type: str, TP_IPCC_bool = False, year_offset: int = 0, base_year_offset = 6):
    # Opens the raster and loads the data 
    with rasterio.open(da_fp) as src:
        age = src.read(1)
        src_nd_value = src.nodata

    # Checks if weather type is in the list
    if weather_type not in forest_litter_table.select(pl.col("IPCC Climate")).to_series().to_list():
        raise ValueError("Weather type not valid")
    
    # Checks if forest type is valid
    if forest_type not in ["NEEV", "BRDC"]:
        raise ValueError("Forest type not valid")
    else:
        forest_id_mean = "BD_mean" if forest_type == "BRDC" else "NE_mean"
        forest_id_tp = "BD_TP" if forest_type == "BRDC" else "NE_TP"
    
    # Gets maturity litter
    res_rate = forest_litter_table.filter(pl.col("IPCC Climate") == weather_type)[forest_id_mean].item()

    # Sets transition period. If IPCC route, deafults to 20, if not, depends on weather and forest type
    if TP_IPCC_bool:
        TP = 20
    else:
        TP = forest_litter_table.filter(pl.col("IPCC Climate") == weather_type)[forest_id_tp].item()

    # Calculates litter
    if src_nd_value is None:
        age_mask = np.isfinite(age)
    else:
        age_mask = ~np.isnan(age) if np.isnan(src_nd_value) else (age != src_nd_value)

    litter = np.where(age_mask, np.minimum(res_rate, res_rate/TP * (age + base_year_offset + year_offset)), np.nan)

    return litter

def get_forest_litter_monthlyrate_fromda(da: np.ndarray, forest_type: str, weather_type: str, TP_IPCC_bool = False, year_offset: int = 0, base_year_offset = 6, residue_runs = 1)-> np.ndarray:
    # this assumes that the da has already been masked properly

    # Checks if weather type is in the list
    if weather_type not in forest_litter_table.select(pl.col("IPCC Climate")).to_series().to_list():
        raise ValueError("Weather type not valid")
    
    # Checks if forest type is valid
    if forest_type not in ["NEEV", "BRDC"]:
        raise ValueError("Forest type not valid")
    else:
        forest_id_mean = "BD_mean" if forest_type == "BRDC" else "NE_mean"
        forest_id_min = "BD_min" if forest_type == "BRDC" else "NE_min"
        forest_id_max = "BD_max" if forest_type == "BRDC" else "NE_max"
        forest_id_tp = "BD_TP" if forest_type == "BRDC" else "NE_TP"
    
    # Gets maturity litter
    res_rate = forest_litter_table.filter(pl.col("IPCC Climate") == weather_type)[forest_id_mean].item()

    # Sets transition period. If IPCC route, deafults to 20, if not, depends on weather and forest type
    if TP_IPCC_bool:
        TP = 20
    else:
        TP = forest_litter_table.filter(pl.col("IPCC Climate") == weather_type)[forest_id_tp].item()

    # Calculates litter
    valid_mask = ~np.isnan(da)

    if residue_runs > 1:
        # draw 100 samples from a triangular distribution (left=min, mode=mean, right=max)
        min_val = float(forest_litter_table.filter(pl.col("IPCC Climate") == weather_type)[forest_id_min].item())
        mode_val = float(forest_litter_table.filter(pl.col("IPCC Climate") == weather_type)[forest_id_mean].item())
        max_val = float(forest_litter_table.filter(pl.col("IPCC Climate") == weather_type)[forest_id_max].item())

        # ensure valid triangular parameters (left <= mode <= right)
        left, mode, right = sorted([min_val, mode_val, max_val])
        samples = np.random.triangular(left, mode, right, size=residue_runs)
        res_rate = float(samples.mean())

    litter = np.where(valid_mask, np.minimum(res_rate, res_rate/TP * (da + base_year_offset + year_offset)), np.nan)
    # build a 12-band array where each month is 1/12 of the annual litter
    monthly_litter = litter/12

    return monthly_litter

#################################
#### GRASSLAND CALCULATIONS #####
#################################
@lru_cache(maxsize=1)
def _load_grassland_residue_table() -> pl.DataFrame:
    """Load the IPCC grassland residue table from disk."""

    table_path = data_path("grasslands", "grassland_residues_IPCC.xlsx")
    try:
        return pl.read_excel(table_path)
    except FileNotFoundError as exc:  # pragma: no cover - optional input tables
        raise FileNotFoundError(
            "Grassland residue lookup table not found at"
            f" {table_path}. Provide an explicit file path or override the input."
        ) from exc


def _resolve_optional_path(path: Optional[Union[str, Path]], *default: str) -> Path:
    """Return a resolved path, falling back to the shared data directory."""

    if path is None:
        return data_path(*default)

    candidate = Path(path)
    if candidate.exists():
        return candidate

    default_path = data_path(*default)
    logging.getLogger(__name__).debug(
        "Falling back to default data path %s for missing input %s", default_path, path
    )
    return default_path


def generate_grassland_residue_map(
    grass_lu_fp: Optional[Union[str, Path]] = None,
    fao_climate_map: Optional[Union[str, Path]] = None,
    c_content: float = 0.47,
    random_runs: int = 1,
) -> np.ndarray:
    # ------- Step 1 - Load maps ----------
    # loading both rasters
    grassland_path = _resolve_optional_path(grass_lu_fp, "land_use", "lu_Grassland.tif")
    with rasterio.open(grassland_path) as src:
        lu = src.read(1)
        lu_nd = src.nodata

    climate_path = _resolve_optional_path(
        fao_climate_map, "soil_weather", "uhth_thermal_climates.tif"
    )
    with rasterio.open(climate_path) as clim:
        clim_data = clim.read(1)
        clim_nd = clim.nodata

    # Creates a mask for land use
    lu_valid = ~np.isnan(lu) & (lu == 1) & (lu != lu_nd)

    # Mask for valud climate data
    climate_valid = (~np.isnan(clim_data)) & (clim_data != clim_nd)

    # Final mask
    grass_clim_valid = lu_valid & climate_valid
    
    # ------- Step 2 - Build Residue Data for Numpy ----------
    # Getting the data
    grassland_residue_table = _load_grassland_residue_table()
    climate_ids = grassland_residue_table["FAO_ID"].to_numpy().astype("int")
    means_above = grassland_residue_table["Residue_Above"].to_numpy().astype("float32")
    means_below = grassland_residue_table["Residue_Below"].to_numpy().astype("float32")
    ses_abv = grassland_residue_table["Res_Err_Abv"].to_numpy().astype("float32")
    ses_below = grassland_residue_table["Res_Err_Below"].to_numpy().astype("float32")

    # Building lookup table
    mean_lut_abv = np.full(13, np.nan, dtype='float32')
    se_lut_abv = np.full(13, np.nan, dtype='float32')
    mean_lut_blw = np.full(13, np.nan, dtype='float32')
    se_lut_blw = np.full(13, np.nan, dtype='float32')

    # Creating lookup tables (lut)
    mean_lut_abv[climate_ids] = means_above
    se_lut_abv[climate_ids] = ses_abv

    mean_lut_blw[climate_ids]   = means_below
    se_lut_blw[climate_ids]     = ses_below

    # ------- Step 3 - Assign residue values ----------
    # Building the arrays
    clim_raster_id = clim_data.astype(int)
    
    pixel_means_above = mean_lut_abv[clim_raster_id]
    pixel_es_above    = se_lut_abv[clim_raster_id]
    
    pixel_means_below = mean_lut_blw[clim_raster_id]
    pixel_es_below    = se_lut_blw[clim_raster_id]

    # Creating the residues including standard error. Assumes normal distribution
    if int(random_runs) <= 1:
        # Deterministic: use the mean only (no SE)
        res_pixel_above = pixel_means_above * 0.5
        res_pixel_below = pixel_means_below * 0.5
    else:
        # Stochastic: average of random_runs normal draws per pixel
        n_runs = int(random_runs)
        draws_above = np.random.normal(
            loc=pixel_means_above,
            scale=pixel_es_above,
            size=(n_runs, *pixel_means_above.shape)  # Construct an array of n_runs, pixel_means same shape
        )

        draws_below = np.random.normal(
            loc=pixel_means_below,
            scale=pixel_es_below,
            size=(n_runs, *pixel_means_below.shape)  # Construct an array of n_runs, pixel_means same shape
        )

        res_pixel_above = draws_above.mean(axis=0) * 0.5
        res_pixel_below = draws_below.mean(axis=0) * 0.5

    res_pixel_total = res_pixel_above + res_pixel_below

    # Finally asigning them to pixels
    grassland_residue = np.full_like(lu, fill_value=np.nan, dtype='float32')
    grassland_residue[grass_clim_valid] = res_pixel_total[grass_clim_valid] * c_content

    return grassland_residue

# -------- Dung Calculations --------------
# Class to handle Dung Calculations
class DungBundle(NamedTuple):
    array: np.ndarray
    nodata: Optional[int]
    dtype: str
    transform: rasterio.Affine
    crs: rasterio.crs.CRS
    profile: dict
    path: str

# Class to store output of dung calculations
from dataclasses import dataclass
@dataclass
class RasterResult:
    array: np.ndarray
    profile: dict
    name: str = ""
    
    def write(self, path: Union[str, Path]) -> None:
        prof = self.profile.copy()
        # Ensure float32 + nodata
        prof.update(dtype="float32", count=1)
        with rasterio.open(path, "w", **prof) as dst:
            dst.write(self.array.astype("float32"), 1)
            if self.name:
                dst.set_band_description(1, self.name)

ANIMAL_DENSITY_REGISTRY = {
    "cattle_other": data_path("grasslands", "livestock", "grassland_cattle.tif"),
    "cattle_dairy": data_path("grasslands", "livestock", "grassland_cattle.tif"),
    "goat": data_path("grasslands", "livestock", "grassland_goat.tif"),
    "sheep": data_path("grasslands", "livestock", "grassland_sheep.tif"),
}

grassland_dung_regions_raster_fp = data_path(
    "grasslands", "livestock", "grassland_dung_regions.tif"
)


@lru_cache(maxsize=1)
def _load_dung_data() -> pl.DataFrame:
    """Load the dung deposition lookup table."""

    dung_path = data_path("grasslands", "Animals_Dung_IPCC.xlsx")
    try:
        return pl.read_excel(
            dung_path, sheet_name="C_Excr_Animals_tCpheadpyr"
        )
    except FileNotFoundError as exc:  # pragma: no cover - optional input tables
        raise FileNotFoundError(
            "Dung excretion lookup table not found at"
            f" {dung_path}. Provide an explicit file path or override the input."
        ) from exc
raster_id = [1,2,3,4,5,6,7,8,9]
dung_regions_mean = [
    "India - Mean",
    "Eastern Europe",
    "Western Europe",
    "Middle East - Mean",
    "North America",
    "LATAM - Mean",
    "Asia - Mean",
    "Africa - Mean",
    "Oceania"  
]
dung_regions_hps = [
    "India - High PS",
    "Eastern Europe",
    "Western Europe",
    "Middle East - High PS",
    "North America",
    "LATAM - High PS",
    "Asia - High PS",
    "Africa - High PS",
    "Oceania" 
]
dung_regions_lps = [
    "India - Low PS",
    "Eastern Europe",
    "Western Europe",
    "Middle East - Low PS",
    "North America",
    "LATAM - Low PS",
    "Asia - Low PS",
    "Africa - Low PS",
    "Oceania" 
]

dung_mean_names =pl.DataFrame(
    {"raster_id":  raster_id,
    "region": dung_regions_mean}
)
dung_hps_names =pl.DataFrame(
    {"raster_id":  raster_id,
    "region": dung_regions_hps}
)
dung_lps_names =pl.DataFrame(
    {"raster_id":  raster_id,
    "region": dung_regions_lps}
)

def calculate_carbon_dung(animals: Union[str, List[str]], cattle_dw_productivity: str = "average"):
    # Check if animal type is valid
    # normalize input
    if isinstance(animals, str): #transform into a list if needed
        animals = [animals]
    animals = [a.lower() for a in animals] # put everything in lower case

     # validate
    unknown = [a for a in animals if a not in ANIMAL_DENSITY_REGISTRY]
    if unknown:
        raise ValueError(f"Unknown animals: {unknown}. Choose among: {list(ANIMAL_DENSITY_REGISTRY)}")
    
    # Checks if developing world pathway is correct
    dw_pathway = ["average", "high", "low"]
    if cattle_dw_productivity not in dw_pathway:
        raise ValueError(f"{cattle_dw_productivity} not valid. Choose between {dw_pathway}")
    
    # Opens dung raster regions
    with rasterio.open(grassland_dung_regions_raster_fp) as region_src:
        dung_regions = region_src.read(1, masked=True)
        valid_regions = ~np.ma.getmaskarray(dung_regions)
        dung_reg_nd  = region_src.nodata
        profile      = region_src.profile

    # --------- Step 1 - Loading data --------------
    animals_density: Dict[str, DungBundle] = {}
    dung_data = _load_dung_data()
    for a in animals:
        fp = ANIMAL_DENSITY_REGISTRY[a]
        with rasterio.open(fp) as src:
            arr = src.read(1)
            nd  = src.nodata
            animals_density[a] = DungBundle(
                array=arr,
                nodata=nd,
                dtype=str(arr.dtype),
                transform=src.transform,
                crs=src.crs,
                profile=src.profile,
                path=fp,
            )

    # --------- Step 2 - calculating dung --------------
    # Loads matching table:
    if cattle_dw_productivity == "average":
        dung_scenario_regions = dung_mean_names
    elif cattle_dw_productivity == "high":
        dung_scenario_regions = dung_hps_names
    else:  # low productivity
        dung_scenario_regions = dung_lps_names 
    

    carbon_out: Dict[str, RasterResult] = {}
    for a in animals:
        # load dung region values for the given scenario
        dung_region_values = dung_data.filter(
            pl.col("region").is_in(dung_scenario_regions["region"].to_list())
        ).select("region", a)

        # create a look up table (lut) to link zones to annual c excretion rates
        lut_df = (dung_scenario_regions.
                  join(dung_region_values, how="left", on="region").
                  select("raster_id", a).rename({a: "annual_c_perha"})
                  )
        
        # 2) Dense LUT array indexed by raster_id
        ids  = lut_df["raster_id"].to_numpy()
        vals = lut_df["annual_c_perha"].to_numpy()
        max_id = 9
        lut = np.full(max_id + 1, np.nan, dtype="float32")
        lut[ids] = vals.astype("float32")    

        # Valid mask
        lsu_km2 = animals_density[a].array
        nd   = animals_density[a].nodata
        has_animals = (lsu_km2 != nd) & ~np.isnan(lsu_km2) & (lsu_km2 >= 0)

        # Masking regions just in case
        dung_has_animals = np.where(has_animals, dung_regions, -999)

        # Create the dung values array
        dung_region_values_array = np.full_like(lsu_km2, fill_value=np.nan, dtype="float32")
        rid = dung_regions.filled(-1).astype(np.int16)

        valid = has_animals & valid_regions

        out_of_bounds = (rid >= 0) & (rid > max_id)
        if np.any(out_of_bounds):
            logging.getLogger(__name__).warning(
                "Encountered dung region ids outside expected range 0-%s; masking %s cells.",
                max_id,
                int(out_of_bounds.sum()),
            )
            valid = valid & ~out_of_bounds

        dung_region_values_array[valid] = lut[rid[valid]]

        # Finally calculate the carbon output
        animal_carbon = np.where(
            valid,
            lsu_km2 / 100 * dung_region_values_array,
            np.nan,
        ) # 100 ha is 1 km2
        
        # Storing results in return array
        out_profile = animals_density[a].profile.copy()
        out_profile.update(
            dtype = "float32",
            nodata = np.nan,
            count = 1
        )
        carbon_out[a] = RasterResult(
            array = animal_carbon.astype("float32"),
            profile = out_profile,
            name = f"{a}_annual_c_t_per_ha"
        )

    
    return carbon_out
