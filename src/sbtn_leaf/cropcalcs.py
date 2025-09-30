############## Crop Calcs ##############
# Script to calculate yields and residues for different crops

#### MODULES ####
from pathlib import Path
import pandas as pd
import polars as pl
import numpy as np
import calendar
import rasterio.fill
import xarray as xr
import os

from collections import defaultdict  # Imports a special dict subclass that auto‐initializes missing keys.
from typing import Optional

import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
import rioxarray as rxr

from sbtn_leaf.PET import monthly_KC_curve, calculate_crop_based_PET_raster_vPipeline
from sbtn_leaf.map_calculations import resample_to_match_noSaving


from rasterio.enums import Resampling


from rasterio.fill import fillnodata
from scipy import ndimage
from typing import Optional, Tuple, Dict


##### DATA ####
er_17 = gpd.read_file("../data/Ecoregions2017/Ecoregions2017.shp")
crops_name_table = pl.read_csv('../data/crops/crop_naming_index.csv')
fao_stats = pl.read_csv('../data/crops/Production_Crops_Livestock_E_All_Data.csv')
fao_crop_yields_1423 = pl.read_csv('../data/crops/fao_crop_yields_1423.csv', separator=';')
country_shp = gpd.read_file('../data/CountryLayers/Country_Level0/g2015_2014_0.shp')
crop_ag_res_table = pl.read_excel("../data/crops/crop_residue_data.xlsx", sheet_name="crop_ABG_Res")
crop_res_table = pl.read_excel("../data/crops/crop_residue_data.xlsx", sheet_name="crop_res_ratios")
rain_monthly_fp = "../data/soil_weather/uhth_monthly_avg_precip.tif"
uhth_climates_fp = "../data/soil_weather/uhth_thermal_climates.tif"

K_Crops = pl.read_csv("../data/crops/K_Crop_Data.csv")
abs_date_table = pl.read_csv("../data/crops/AbsoluteDayTable.csv")

days_in_month_table = pl.DataFrame({
    "Month": list(range(1, 13)),
    "Days_in_Month": [calendar.monthrange(2023, month)[1] for month in range(1, 13)]
})


thermal_climates_table = pl.DataFrame({
    "id": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    "TC_Name": [
        "Tropics, lowland",
        "Tropics, highland",
        "Subtropics, summer rainfall",
        "Subtropics, winter rainfall",
        "Subtropics, low rainfall",
        "Temperate, oceanic",
        "Temperate, sub-continental",
        "Temperate, continental",
        "Boreal, oceanic",
        "Boreal, sub-continental",
        "Boreal, continental",
        "Arctic"
    ],
    "TC_Group": [
        "Tropics",
        "Tropics",
        "Subtropics summer rainfall",
        "Subtropics winter rainfall",
        "Subtropics winter rainfall",
        "Oceanic temperate",
        "Sub-continental temperate and continental temperate",
        "Sub-continental temperate and continental temperate",
        "Sub-continental boreal, continental boreal and polar/arctic",
        "Sub-continental boreal, continental boreal and polar/arctic",
        "Sub-continental boreal, continental boreal and polar/arctic",
        "Sub-continental boreal, continental boreal and polar/arctic"
    ]
})


# Dictionaries
# Create a lookup dictionary for thermal climates
# This will map thermal climate IDs to their TC_Group names.
climate_zone_lookup = dict(
    zip(thermal_climates_table["id"].to_list(),
        thermal_climates_table["TC_Group"].to_list())
)

# Dictionary to group climate zones by their TC_Group
# This will help in quickly accessing all zone IDs belonging to a specific group.
zone_ids_by_group = defaultdict(list)
for zone_id, group in climate_zone_lookup.items():
    zone_ids_by_group[group].append(zone_id)  # Append zone ID to the list for its group


#### FUNCTIONS ####
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

def create_crop_yield_shapefile(fao_crop: str):
    # Checks if the crop is in the list
    if fao_crop not in crops_name_table['FAO_Crop'].unique():
        raise ValueError(f'{fao_crop} not found or has no data')

    # Extract needed data
    yields_df = (
        fao_crop_yields_1423.filter(pl.col("Item") == fao_crop)
        .select(["Area", "Unit", "avg_yield_1423", "ratio_yield_20_toavg"])
        .to_pandas()
    )

    # rename to shorter namies
    yields_df = yields_df.rename(columns={"avg_yield_1423": "avg_yield", "ratio_yield_20_toavg": "yld_ratio"})

    # Merges with shapefile
    yield_shp = country_shp.merge(yields_df, how='left', left_on='ADM0_NAME', right_on='Area').drop(columns='Area')

    return yield_shp


def create_crop_yield_raster(croplu_grid_raster: str, fao_crop_shp: gpd.GeoDataFrame, spam_crop_raster: str, output_rst_path: str, 
                             spam_band=1, resampling_method = Resampling.bilinear):

    # 1--- ) Opens the master_grid_raster and extracts metadata, transform, crs, and dimensions ---
    with rasterio.open(croplu_grid_raster) as crop_lu:
        lu_meta         = crop_lu.meta.copy()  # Gets metadata
        lu_crs          = crop_lu.crs
        lu_transform    = crop_lu.transform
        lu_width        = crop_lu.width
        lu_height       = crop_lu.height
        lu_data         = crop_lu.read(1)

    # Create a land use mask
    lu_mask =  lu_data == 1
        
    # --- 2) Opens the spam_crop_raster ---
    with rasterio.open(spam_crop_raster) as spam:
        spam_data = spam.read(spam_band)
        spam_nodata = spam.nodata

        # prepare empty array for reprojected SPAM with the same dimensions as the land use grid
        spam_on_lu = np.full(shape=(lu_height, lu_width), fill_value=np.nan, dtype='float32')

        # --- 3) fills the array with the reprojected spam values ---
        reproject(
         source=spam_data,
         destination=spam_on_lu,
         src_transform=spam.transform,
         src_crs=spam.crs,
         src_nodata=spam_nodata,
         dst_transform=lu_transform,
         dst_crs=lu_crs,
         dst_nodata=np.nan,
         resampling=resampling_method
        )

    # --- 4) Load FAO shapefile, rasterize zones ---
    fao_yields_gdf = fao_crop_shp
    fao_yields_gdf = fao_yields_gdf.to_crs(lu_crs)  # Reprojects to land use raster
    
    # Checks if the data columns are there 
    for field in ("avg_yield_1423", "ratio_yield_20_toavg"):
        if field not in fao_yields_gdf.columns:
            raise KeyError(f"Shapefile missing '{field}'")
        
    # Now rasterizing    
    fao_yields_gdf = fao_yields_gdf.reset_index(drop=True)  # Eliminates the old index and assigned 0, 1, 2, ...
    fao_yields_gdf["zone_id"] = fao_yields_gdf.index.astype("int32")  # Creates a zone_ide for the rasterization
    global_fao_ratio = fao_yields_gdf["ratio_yield_20_toavg"].dropna().mean()  # Calcualtes a global average, disregarding nans
    fao_yields_gdf["avg_yield_1423"] = fao_yields_gdf["avg_yield_1423"]/1000          # Transform into tons
    
    # Rasterizes the shapefile
    shapes = ((geom, zid) for geom, zid in zip(fao_yields_gdf.geometry, fao_yields_gdf.zone_id))
    zone_array = rasterize(
        shapes=shapes,
        out_shape=(lu_height, lu_width),
        transform=lu_transform,
        fill=-1,
        dtype="int32"
    )

    # --- 4) Build result, two‐stage fill ---
    # Creates an array of the same shape as land use
    result = np.full((lu_height, lu_width), np.nan, dtype="float32")
    
    for _, row in fao_yields_gdf.iterrows():
        zid       = row["zone_id"]
        ratio     = row["ratio_yield_20_toavg"]
        avg_yield = row["avg_yield_1423"]/1000  # Dividing by a 1,000 to transform to tonnes

        mask_zone       = zone_array == zid  # Creates a mask for each zone id
        spam_scaled     = spam_on_lu * ratio  # Calculates the scaled yield
        
        # Assigns spam values for areas that have a spam value
        mask_spam_valid = mask_zone & ~np.isnan(spam_scaled)  # Creates a mask for values for each zone id that is not empty on the spam raster
        result[mask_spam_valid] = spam_scaled[mask_spam_valid]  # Assings into the results raster the value of the spam zone
        
        # Assigns fao values when spam doesn't 
        mask_need_avg          = mask_zone & np.isnan(result)  # Creates a mask for the mask zone where the results is empty
        result[mask_need_avg]   = avg_yield  # Assigns the average yield for the empty pixels

    # There could be pixels where FAO had no values, but SPAM had them
    mask_spam_fallback = (~np.isnan(spam_on_lu)) & np.isnan(result) & lu_mask
    result[mask_spam_fallback] = spam_on_lu[mask_spam_fallback] * global_fao_ratio
    
    # Mask out where there's no land use
    result[~lu_mask] = np.nan

    # --- 5) Write out ---
    lu_meta.update(dtype="float32", count=1, nodata=np.nan)
    with rasterio.open(output_rst_path, "w", **lu_meta) as dst:
        dst.write(result, 1)

    print(f"Yield raster written to {output_rst_path}")


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
    fao_yield_ratio_name: str = "yld_ratio"
):
    """
    Create a crop yield raster by combining SPAM data with FAO shapefile yields,
    optionally applying irrigation/rainfed scaling to FAO averages.
    """
    # 1) Open cropland LU raster
    with rasterio.open(croplu_grid_raster) as crop_lu:
        lu_meta      = crop_lu.meta.copy()
        lu_crs       = crop_lu.crs
        lu_transform = crop_lu.transform
        lu_height    = crop_lu.height
        lu_width     = crop_lu.width
        lu_data      = crop_lu.read(1)

    lu_mask = (lu_data == 1)

    # 2) Reproject SPAM onto LU grid
    with rasterio.open(spam_crop_raster) as spam:
        spam_data   = spam.read(spam_band)
        spam_on_lu  = np.full((lu_height, lu_width), np.nan, dtype="float32")
        reproject(
            source=spam_data,
            destination=spam_on_lu,
            src_transform=spam.transform,
            src_crs=spam.crs,
            src_nodata=spam.nodata,
            dst_transform=lu_transform,
            dst_crs=lu_crs,
            dst_nodata=np.nan,
            resampling=resampling_method
        )

    # 3) Rasterize FAO yields & ratios
    fao_gdf = fao_crop_shp.to_crs(lu_crs).reset_index(drop=True)
    for field in (fao_avg_yield_name, fao_yield_ratio_name):
        if field not in fao_gdf.columns:
            raise KeyError(f"Missing '{field}' in FAO shapefile")

    # convert avg_yield to tons
    fao_gdf[fao_avg_yield_name] = fao_gdf[fao_avg_yield_name] / 1000.0
    global_fao_ratio = fao_gdf[fao_yield_ratio_name].dropna().mean()

    fao_gdf["zone_id"] = fao_gdf.index.astype("int32")
    shapes = ((geom, zid) for geom, zid in zip(fao_gdf.geometry, fao_gdf.zone_id))
    zone_array = rasterize(
        shapes=shapes,
        out_shape=(lu_height, lu_width),
        transform=lu_transform,
        fill=-1,
        dtype="int32"
    )

    # build per-zone FAO arrays
    fao_yields_array = np.full((lu_height, lu_width), np.nan, dtype="float32")
    fao_ratios_array = np.full((lu_height, lu_width), np.nan, dtype="float32")
    for _, row in fao_gdf.iterrows():
        zid = int(row["zone_id"])
        zid_mask = (zone_array == zid)
        fao_yields_array[zid_mask] = row[fao_avg_yield_name]
        fao_ratios_array[zid_mask] = row[fao_yield_ratio_name]

    valid_fao = ~np.isnan(fao_yields_array)

    # --- Optional: apply irrigation/rainfed scaling to FAO yields ---
    if irr_yield_scaling is not None:
        if any(p is None for p in (all_fp, irr_fp, rf_fp)):
            raise ValueError("Need all_fp, irr_fp and rf_fp for irrigation scaling")

        # resample modifiers
        all_fp_on_lu = resample_to_match_noSaving(all_fp, croplu_grid_raster, dst_nodata = np.nan)
        irr_fp_on_lu = resample_to_match_noSaving(irr_fp,croplu_grid_raster, dst_nodata = np.nan)
        rf_fp_on_lu  = resample_to_match_noSaving(rf_fp, croplu_grid_raster,dst_nodata = np.nan)

        # compute ratios
        irr_ratios, rf_ratios = calculate_SPAM_yield_modifiers(all_fp_on_lu, irr_fp_on_lu, rf_fp_on_lu)
        watering_ratio = irr_ratios if irr_yield_scaling == 'irr' else rf_ratios

        valid_wat = ~np.isnan(watering_ratio)
        avg_wat   = np.nanmean(watering_ratio)

        # scale FAO yields
        scaled = np.where(valid_wat, fao_yields_array * watering_ratio, np.nan)
        scaled = fillnodata(scaled, mask=np.isnan(scaled), max_search_distance=1, smoothing_iterations=2)
        scaled = np.where(np.isnan(scaled) & valid_fao,
                          fao_yields_array * avg_wat,
                          scaled)
        fao_yields_array = scaled

    # 4) Build result: SPAM first, then FAO, then SPAM‐fallback
    result = np.full((lu_height, lu_width), np.nan, dtype="float32")

    # a) per-zone SPAM scaling
    for _, row in fao_gdf.iterrows():
        zid   = int(row["zone_id"])
        zid_mask  = (zone_array == zid)
        ratio = row[fao_yield_ratio_name]
        spam_scaled = spam_on_lu * ratio
        valid_mask = zid_mask & ~np.isnan(spam_scaled)
        result[valid_mask] = spam_scaled[valid_mask]

    # a.5) all‐SPAM irrigation scaling (only if requested)
    if irr_yield_scaling is not None:
        print("  → Applying irrigation scaling to all‐SPAM yields…")
        mask_all = (
            lu_mask
            & np.isnan(result)                       # still empty
            & ~np.isnan(all_fp_on_lu)                # has an “all” spam value
        )
        result[mask_all] = all_fp_on_lu[mask_all] * avg_wat

    # b) FAO fill for remaining
    mask_fao = lu_mask & np.isnan(result) & valid_fao
    result[mask_fao] = fao_yields_array[mask_fao]

    # c) SPAM fallback
    mask_spam = (~np.isnan(spam_on_lu)) & np.isnan(result) & lu_mask
    result[mask_spam] = spam_on_lu[mask_spam] * global_fao_ratio

    # compute region & biome averages
    ecoregion_avg, biome_avg, zone_array, biome_name_map = \
        calculate_average_yield_by_ecoregion_and_biome(
            result, croplu_grid_raster
        )

    # d) assign per-pixel ecoregion or biome average for remaining holes
    remaining = lu_mask & np.isnan(result)
    if np.any(remaining):
        ys, xs = np.where(remaining)
        for y, x in zip(ys, xs):
            zid = int(zone_array[y, x])
            if zid in ecoregion_avg:
                result[y, x] = ecoregion_avg[zid]
            else:
                biome = biome_name_map.get(zid)
                # Ensure biome is a string before using as key
                if isinstance(biome, str):
                    result[y, x] = biome_avg.get(biome, global_fao_ratio)
                else:
                    result[y, x] = global_fao_ratio

    # e) final hole fill by nearest neighbor
    remaining = lu_mask & np.isnan(result)
    if np.any(remaining):
        valid = ~np.isnan(result)
        dist, (iy, ix) = ndimage.distance_transform_edt(
            ~valid, return_distances=True, return_indices=True
        )
        filled = result[iy, ix]
        result[remaining] = filled[remaining]


    # mask non‐cropland
    result[~lu_mask] = np.nan

    # 5) write out
    lu_meta.update(dtype="float32", count=1, nodata=np.nan)
    with rasterio.open(output_rst_path, "w", **lu_meta) as dst:
        dst.write(result, 1)

    print(f"Yield raster written to {output_rst_path}")


def calculate_SPAM_yield_modifiers(
    all_yields: np.ndarray,
    irr_yields: np.ndarray,
    rf_yields: np.ndarray,
    save_ratios: bool = False,
    all_rasters_fp: Optional[str] = None,
    irr_ratios_fp: Optional[str] = None,
    rf_ratios_fp: Optional[str] = None
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
    # pre‑fill result with NaNs
    irr_ratios = np.full_like(all_yields, np.nan, dtype=float)
    rf_ratios  = np.full_like(all_yields, np.nan, dtype=float)
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

    print(f"Average irrigated ratio: {np.nanmean(irr_ratios)}")
    print(f"Average rainfed ratio: {np.nanmean(rf_ratios)}")

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
    er_shapefile = er_17
) -> Tuple[Dict[int, float], Dict[str, float], np.ndarray, Dict[int, str]]:
    """
    Calculate average yields per ecoregion and per biome.
    Returns:
      - ecoregion_avg: zone_id -> average yield
      - biome_avg: biome_name -> average yield
      - zone_array: rasterized zone_id array
      - biome_name_map: zone_id -> biome_name mapping
    """
    # load LU raster for transform & CRS
    with rasterio.open(croplu_grid_raster) as src:
        transform, crs = src.transform, src.crs
        height, width = src.height, src.width

    # load and reproject ecoregions
    er_gdf = er_17.to_crs(crs).reset_index(drop=True)
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
    
    if crop not in (crops_name_table["Crops"].to_list() + crops_name_table["IPCC_Crop"].to_list()):
        raise ValueError(f"{crop} not found in data table")
    
    # Initialize crop amounts
    ABG = 0.0
    BG = 0.0
    Res = ABG + BG

    # Getting IPCC name and calculating belowground
    ipcc_crop = crops_name_table.filter(pl.col('Crops') == crop).select('IPCC_Crop').item()
    res_crop_data = crop_res_table.filter(pl.col('Crop') == ipcc_crop)
    dry = res_crop_data.select('DRY').item()
    RS = res_crop_data.select("RS").item()

    # Checkes if ABG can be calculated with line equation
    if crop in crop_ag_res_table["Crop"].to_list():
        AG_crop_data = crop_ag_res_table.filter(pl.col("Crop") == crop)
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
    ipcc_crop = (
        crops_name_table
        .filter(pl.col("Crops") == crop)
        .select("IPCC_Crop")
        .to_series()
        .item()
    )

    # 3) Pull core residue parameters
    res_row = crop_res_table.filter(pl.col("Crop") == ipcc_crop)
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
    ag_row = crop_ag_res_table.filter(pl.col("Crop") == crop)
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


def create_plant_cover_monthly_curve(crop: str, climate: str):
    # Check the crops exists
    if crop not in K_Crops['Crop'].unique():
        raise ValueError(f"Crop '{crop}' not found in K_Crops table.")
    
    if climate not in K_Crops["Climate_Zone"]:
        raise ValueError(f"Climate zone '{climate}' not found in K_Crops table.")
    
    # Retrieve plant cover data for the specified crop
    pc_starts = K_Crops.filter((pl.col('Crop') == crop) & (pl.col('Climate_Zone') == climate)).select('SCP_Starts').item()
    pc_ends = K_Crops.filter((pl.col('Crop') == crop) & (pl.col('Climate_Zone') == climate)).select('SCP_End').item()

    # Create a DataFrame for the plant cover curve
    plant_cover_curve = pl.DataFrame(
        {
            "Month": list(range(1, 13)),
            "Plant_Cover": [0] * 12
        }
    )

    # Fill the plant cover curve based on start and end dates
    plant_cover_curve = plant_cover_curve.with_columns(
        pl.when((pl.col('Month')>=pc_starts), (pl.col('Month')<=pc_ends)).then(1).otherwise(0).alias('Plant_Cover')
    )

    return plant_cover_curve


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
    output_nodata: int = 255
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

    # 4. Loop over each unique climate ID
    unique_ids = np.unique(ids[~np.isnan(ids)]).astype(int)
    for cid in unique_ids:
        group = climate_zone_lookup.get(cid)
        if group is None:
            continue  # or raise
        # Get the 12-month vector (0/1) from your existing function
        pc_df = create_plant_cover_monthly_curve(crop, group)
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

def compute_monthly_residue_raster(
    crop: str,
    climate_raster_path: str,
    plant_residue: xr.DataArray,
    save_path: str,
    output_nodata: float = np.nan
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
    # 1. Load & squeeze the climate IDs
    clim = rxr.open_rasterio(climate_raster_path, masked=False)
    if "band" in clim.dims:
        clim = clim.isel(band=0)
    ids = clim.values  # now strictly 2D (y,x)
    y, x = clim.coords["y"], clim.coords["x"]

    # 2. Squeeze the plant-residue to 2D, always taking the first band if there are multiples
    pr_da = plant_residue
    if "band" in pr_da.dims:
        pr_da = pr_da.isel(band=0)
    if pr_da.ndim != 2:
        raise ValueError(f"Expected plant_residue to be 2D after squeezing, got {pr_da.shape}")
    pr = pr_da.values

    # 2. Prepare output array
    n_months = 12
    y, x = pr_da.coords["y"], pr_da.coords["x"]
    out = np.full((n_months, y.size, x.size),
                  fill_value=output_nodata,
                  dtype="float32")

    # 3. Loop over climate IDs
    unique_ids = np.unique(ids[~np.isnan(ids)]).astype(int)
    for cid in unique_ids:
        group = climate_zone_lookup.get(cid)
        if group is None:
            continue

        # 3a. K‐curve fraction
        kc_df = monthly_KC_curve(crop, group).select(["Month", "Kc"])
        k_vals = np.array(kc_df["Kc"].to_list(), dtype=float)
        k_frac = k_vals / k_vals.sum()

        # 3b. Pixel positions for this ID
        rows, cols = np.where(ids == cid)  # now ids is 2D

        # 3c. Annual residue values at those pixels
        pr_vals = pr[rows, cols]  # now works, pr is 2D

        # 3d. Assign monthly: shape (12, Npix)
        out[:, rows, cols] = k_frac[:, None] * pr_vals[None, :]

    # 4. Wrap as xarray and set metadata
    da = xr.DataArray(
        out,
        dims=("month", "y", "x"),
        coords={
            "month": np.arange(1, n_months+1),
            "y": y,
            "x": x
        },
        name=f"{crop}_residue_monthly"
    ).astype("float32")

    da = da.rio.write_crs(clim.rio.crs)
    da = da.rio.write_transform(clim.rio.transform())
    da = da.rio.write_nodata(output_nodata)
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    # 5. Optionally save as a 12-band GeoTIFF
    da.rio.to_raster(save_path)


def compute_monthly_residue_raster_fromAnnualRaster(
    crop: str,
    plant_residue: str,
    save_path: str,
    climate_raster_path: str = uhth_climates_fp,
    output_nodata: float = np.nan
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
    # 1. Load & squeeze the climate IDs
    clim = rxr.open_rasterio(climate_raster_path, masked=False)
    if "band" in clim.dims:
        clim = clim.isel(band=0)
    ids = clim.values  # now strictly 2D (y,x)
    y, x = clim.coords["y"], clim.coords["x"]

    # 2. Squeeze the plant-residue to 2D, always taking the first band if there are multiples
    pr_da = rxr.open_rasterio(plant_residue, masked=True)
    if "band" in pr_da.dims:
        pr_da = pr_da.isel(band=0)
    if pr_da.ndim != 2:
        raise ValueError(f"Expected plant_residue to be 2D after squeezing, got {pr_da.shape}")
    pr = pr_da.values

    # 2. Prepare output array
    n_months = 12
    y, x = pr_da.coords["y"], pr_da.coords["x"]
    out = np.full((n_months, y.size, x.size),
                  fill_value=output_nodata,
                  dtype="float32")

    # 3. Loop over climate IDs
    unique_ids = np.unique(ids[~np.isnan(ids)]).astype(int)
    for cid in unique_ids:
        group = climate_zone_lookup.get(cid)
        if group is None:
            continue

        # 3a. K‐curve fraction
        kc_df = monthly_KC_curve(crop, group).select(["Month", "Kc"])
        k_vals = np.array(kc_df["Kc"].to_list(), dtype=float)
        k_frac = k_vals / k_vals.sum()

        # 3b. Pixel positions for this ID
        rows, cols = np.where(ids == cid)  # now ids is 2D

        # 3c. Annual residue values at those pixels
        pr_vals = pr[rows, cols]  # now works, pr is 2D

        # 3d. Assign monthly: shape (12, Npix)
        out[:, rows, cols] = k_frac[:, None] * pr_vals[None, :]

    # 4. Wrap as xarray and set metadata
    da = xr.DataArray(
        out,
        dims=("month", "y", "x"),
        coords={
            "month": np.arange(1, n_months+1),
            "y": y,
            "x": x
        },
        name=f"{crop}_residue_monthly"
    ).astype("float32")

    da = da.rio.write_crs(clim.rio.crs)
    da = da.rio.write_transform(clim.rio.transform())
    da = da.rio.write_nodata(output_nodata)
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    # 5. Optionally save as a 12-band GeoTIFF
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
    crop_practice_string: str,
    lu_data_path: str,
    spam_crop_raster: str,
    output_data_folder: str,
    irr_yield_scaling: str,
    spam_all_fp: str,
    spam_irr_fp: str,
    spam_rf_fp: str,
    all_new_files=False,
):

    # Output saving string bases
    output_crop_based = f"{output_data_folder}/{crop_name}"
    output_practice_based = f"{output_crop_based}_{crop_practice_string}"

    # Step 0 - rasterize input path
    lu_bin_output = f"{output_practice_based}_lu.tif"
    if all_new_files or not os.path.exists(lu_bin_output):
        print("Creating lu raster...")
        lu_array = binarize_raster_pipeline(lu_data_path, lu_bin_output)
    else:
        print("Land use binary raster already exist. Skipping...")
        lu_array = rxr.open_rasterio(lu_bin_output, masked=False).squeeze()

    # Step 1 - Prepare PET and irrigation
    # Step 1.1 - PET
    pet_monthly_output_path = f"{output_crop_based}_pet_monthly.tif"
    if all_new_files or not os.path.exists(pet_monthly_output_path):
        print("Creating PET raster...")    
        monthly_pet = calculate_crop_based_PET_raster_vPipeline(
            crop_name=crop_name,
            landuse_array=lu_array,
            output_monthly_path = pet_monthly_output_path
        )
    else:
        print("PET raster already exists — skipping computation.")

    # Step 1.2 - Irrigation
    irr_monthly_output_path = f"{output_crop_based}_irr_monthly.tif"
    if all_new_files or not os.path.exists(irr_monthly_output_path):
        print("Creating irrigation raster...")
        irr = calculate_irrigation_vPipeline(
            evap=monthly_pet,
            output_path=irr_monthly_output_path
        )
    else:
        print("Irrigation raster already exists — skipping computation.")

    # Step 2 - Calculate yields
    # preparing fao yield shapefile
    fao_crop_name = crops_name_table.filter(pl.col("Crops")== crop_name).select(pl.col("FAO_Crop")).item()
    print(f"Creating {fao_crop_name} shapefile...")
    fao_yield_shp = create_crop_yield_shapefile(fao_crop_name)

    # Create irrigation adjusted yields
    yield_output_path = f"{output_practice_based}_yield.tif"
    if all_new_files or not os.path.exists(yield_output_path):
        print("Creating yield raster...")
        create_crop_yield_raster_withIrrigationPracticeScaling_vPipeline(
            croplu_grid_raster=lu_bin_output,
            fao_crop_shp=fao_yield_shp,
            spam_crop_raster=spam_crop_raster,
            output_rst_path=yield_output_path,
            irr_yield_scaling=irr_yield_scaling,
            all_fp = spam_all_fp,
            irr_fp = spam_irr_fp,
            rf_fp= spam_rf_fp
        )
    else:
        print("Yields raster already exists — skipping computation.")

    # Step 3 - Create plant cover raster
    plantcover_output_path = f"{output_crop_based}_pc_monthly.tif"
    if all_new_files or not os.path.exists(plantcover_output_path):
        print("Creating plant cover raster...")
        create_plant_cover_monthly_raster(crop_name, plantcover_output_path)
    else:
        print("Plant Cover raster already exists — skipping computation.")

    # Step 4 - Create plant residue raster
    plant_residue_output_path = f"{output_practice_based}_residues_monthly.tif"
    if all_new_files or not os.path.exists(plant_residue_output_path):
        print("Creating plant residue raster...")
        create_monthly_residue_vPipeline(crop_name, yield_raster_path=yield_output_path, output_path=plant_residue_output_path)
    else:
        print("Plant Residues raster already exists — skipping computation.")

    print(f"All data created for {crop_name}, {crop_practice_string}!!!")


#### PIPELINE SUPPORTING FUNCTIONS
def binarize_raster_pipeline(
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

    # Transform data into array if needed
    with rasterio.open(rain_fp) as src:
        rain = src.read().astype("float32")
        rain_crs = src.crs
        rain_profile = src.profile

    rain = np.asarray(rain, dtype=float)
    evap = np.asarray(evap, dtype=float)

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


def create_crop_yield_raster_withIrrigationPracticeScaling_vPipeline(
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
    fao_yield_ratio_name: str = "yld_ratio"
):
    """
    Create a crop yield raster by combining SPAM data with FAO shapefile yields,
    optionally applying irrigation/rainfed scaling to FAO averages.
    """
    # 1) Open cropland LU raster
    with rasterio.open(croplu_grid_raster) as crop_lu:
        lu_meta      = crop_lu.meta.copy()
        lu_crs       = crop_lu.crs
        lu_transform = crop_lu.transform
        lu_height    = crop_lu.height
        lu_width     = crop_lu.width
        lu_data      = crop_lu.read(1)

    lu_mask = (lu_data == 1)

    # 2) Reproject SPAM onto LU grid
    with rasterio.open(spam_crop_raster) as spam:
        spam_data   = spam.read(spam_band)
        spam_on_lu  = np.full((lu_height, lu_width), np.nan, dtype="float32")
        reproject(
            source=spam_data,
            destination=spam_on_lu,
            src_transform=spam.transform,
            src_crs=spam.crs,
            src_nodata=spam.nodata,
            dst_transform=lu_transform,
            dst_crs=lu_crs,
            dst_nodata=np.nan,
            resampling=resampling_method
        )

    # 3) Rasterize FAO yields & ratios
    fao_gdf = fao_crop_shp.to_crs(lu_crs).reset_index(drop=True)
    for field in (fao_avg_yield_name, fao_yield_ratio_name):
        if field not in fao_gdf.columns:
            raise KeyError(f"Missing '{field}' in FAO shapefile")

    # convert avg_yield to tons
    fao_gdf[fao_avg_yield_name] = fao_gdf[fao_avg_yield_name] / 1000.0
    global_fao_ratio = fao_gdf[fao_yield_ratio_name].dropna().mean()

    fao_gdf["zone_id"] = fao_gdf.index.astype("int32")
    shapes = ((geom, zid) for geom, zid in zip(fao_gdf.geometry, fao_gdf.zone_id))
    zone_array = rasterize(
        shapes=shapes,
        out_shape=(lu_height, lu_width),
        transform=lu_transform,
        fill=-1,
        dtype="int32"
    )

    # build per-zone FAO arrays
    fao_yields_array = np.full((lu_height, lu_width), np.nan, dtype="float32")
    fao_ratios_array = np.full((lu_height, lu_width), np.nan, dtype="float32")
    for _, row in fao_gdf.iterrows():
        zid = int(row["zone_id"])
        zid_mask = (zone_array == zid)
        fao_yields_array[zid_mask] = row[fao_avg_yield_name]
        fao_ratios_array[zid_mask] = row[fao_yield_ratio_name]

    valid_fao = ~np.isnan(fao_yields_array)

    # --- Optional: apply irrigation/rainfed scaling to FAO yields ---
    if irr_yield_scaling is not None:
        if any(p is None for p in (all_fp, irr_fp, rf_fp)):
            raise ValueError("Need all_fp, irr_fp and rf_fp for irrigation scaling")

        # resample modifiers
        all_fp_on_lu = resample_to_match_noSaving(all_fp, croplu_grid_raster, dst_nodata = np.nan)
        irr_fp_on_lu = resample_to_match_noSaving(irr_fp,croplu_grid_raster, dst_nodata = np.nan)
        rf_fp_on_lu  = resample_to_match_noSaving(rf_fp, croplu_grid_raster,dst_nodata = np.nan)

        # compute ratios
        irr_ratios, rf_ratios = calculate_SPAM_yield_modifiers(all_fp_on_lu, irr_fp_on_lu, rf_fp_on_lu)
        watering_ratio = irr_ratios if irr_yield_scaling == 'irr' else rf_ratios

        valid_wat = ~np.isnan(watering_ratio)
        avg_wat   = np.nanmean(watering_ratio)

        # scale FAO yields
        scaled = np.where(valid_wat, fao_yields_array * watering_ratio, np.nan)
        scaled = fillnodata(scaled, mask=np.isnan(scaled), max_search_distance=1, smoothing_iterations=2)
        scaled = np.where(np.isnan(scaled) & valid_fao,
                          fao_yields_array * avg_wat,
                          scaled)
        fao_yields_array = scaled

    # 4) Build result: SPAM first, then FAO, then SPAM‐fallback
    result = np.full((lu_height, lu_width), np.nan, dtype="float32")

    # a) per-zone SPAM scaling
    for _, row in fao_gdf.iterrows():
        zid   = int(row["zone_id"])
        zid_mask  = (zone_array == zid)
        ratio = row[fao_yield_ratio_name]
        spam_scaled = spam_on_lu * ratio
        valid_mask = zid_mask & ~np.isnan(spam_scaled)
        result[valid_mask] = spam_scaled[valid_mask]

    # a.5) all‐SPAM irrigation scaling (only if requested)
    if irr_yield_scaling is not None:
        irrigation_scaling_type = "irrigation" if irr_yield_scaling == "irr" else "rainfed"
        print(f"  → Applying {irrigation_scaling_type} scaling to all‐SPAM yields…")
        mask_all = (
            lu_mask
            & np.isnan(result)                       # still empty
            & ~np.isnan(all_fp_on_lu)                # has an “all” spam value
        )
        result[mask_all] = all_fp_on_lu[mask_all] * avg_wat

    # b) FAO fill for remaining
    mask_fao = lu_mask & np.isnan(result) & valid_fao
    result[mask_fao] = fao_yields_array[mask_fao]

    # c) SPAM fallback
    mask_spam = (~np.isnan(spam_on_lu)) & np.isnan(result) & lu_mask
    result[mask_spam] = spam_on_lu[mask_spam] * global_fao_ratio

    # compute region & biome averages
    ecoregion_avg, biome_avg, zone_array, biome_name_map = \
        calculate_average_yield_by_ecoregion_and_biome(
            result, croplu_grid_raster
        )


    # d) assign per-pixel ecoregion or biome average for remaining holes
    remaining = lu_mask & np.isnan(result)
    if np.any(remaining):
        ys, xs = np.where(remaining)
        for y, x in zip(ys, xs):
            zid = int(zone_array[y, x])
            if zid in ecoregion_avg:
                result[y, x] = ecoregion_avg[zid]
            else:
                biome = biome_name_map.get(zid)
                # Ensure biome is a string before using as key
                if isinstance(biome, str):
                    result[y, x] = biome_avg.get(biome, global_fao_ratio)
                else:
                    result[y, x] = global_fao_ratio

    # e) final hole fill by nearest neighbor
    remaining = lu_mask & np.isnan(result)
    if np.any(remaining):
        valid = ~np.isnan(result)
        dist, (iy, ix) = ndimage.distance_transform_edt(
            ~valid, return_distances=True, return_indices=True
        )
        filled = result[iy, ix]
        result[remaining] = filled[remaining]


    # mask non‐cropland
    result[~lu_mask] = np.nan

    # 5) write out
    lu_meta.update(dtype="float32", count=1, nodata=np.nan)
    with rasterio.open(output_rst_path, "w", **lu_meta) as dst:
        dst.write(result, 1)

    print(f"Yield raster written to {output_rst_path}")


def create_monthly_residue_vPipeline(
    crop: str,
    yield_raster_path: str,
    output_path: str,
    output_nodata = np.nan,
    climate_raster_path: str = uhth_climates_fp,
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
        shape   = src.shape
        yields  = src.read(1).astype("float32")
        nodata  = src.nodata
        src_crs = src.crs
        src_transform = src.transform

    # mask nodata → NaN
    valid = (yields != nodata)
    yld_arr = np.where(valid, yields, np.nan)

    # 2) Map user crop → IPCC crop key
    ipcc_crop = (
        crops_name_table
        .filter(pl.col("Crops") == crop)
        .select("IPCC_Crop")
        .to_series()
        .item()
    )

    # 3) Pull core residue parameters
    res_row = crop_res_table.filter(pl.col("Crop") == ipcc_crop)
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
    ag_row = crop_ag_res_table.filter(pl.col("Crop") == crop)
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
    
    # Prepare output array
    n_months = 12
    y, x = shape
    out = np.full((n_months, y, x),
                  fill_value=output_nodata,
                  dtype="float32")


    # 1. Load & squeeze the climate IDs
    clim = rxr.open_rasterio(climate_raster_path, masked=False)
    if "band" in clim.dims:
        clim = clim.isel(band=0)
    ids = clim.values  # now strictly 2D (y,x)

    # 3. Loop over climate IDs
    unique_ids = np.unique(ids[~np.isnan(ids)]).astype(int)
    for cid in unique_ids:
        group = climate_zone_lookup.get(cid)
        if group is None:
            continue

        # 3a. K‐curve fraction
        kc_df = monthly_KC_curve(crop, group).select(["Month", "Kc"])
        k_vals = np.array(kc_df["Kc"].to_list(), dtype=float)
        k_frac = k_vals / k_vals.sum()

        # 3b. Pixel positions for this ID
        rows, cols = np.where(ids == cid)  # now ids is 2D

        # 3c. Annual residue values at those pixels
        pr_vals = Res[rows, cols]  # now works, pr is 2D

        # 3d. Assign monthly: shape (12, Npix)
        out[:, rows, cols] = k_frac[:, None] * pr_vals[None, :]

    # 4. Wrap as xarray and set metadata
    da = xr.DataArray(
        out,
        dims=("month", "y", "x"),
        coords={
            "month": np.arange(1, n_months+1),
            "y":   clim.coords["y"],
            "x":   clim.coords["x"],
        },
        name=f"{crop}_residue_monthly",
    ).astype("float32")

    # 1️⃣ Tell rioxarray which dims are spatial
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    # 2️⃣ Write CRS, transform, nodata in order
    da.rio.to_raster(
        output_path,
        driver="GTiff",
        crs=src_crs,
        transform=src_transform,
        dtype="float32",
        nodata=output_nodata,
    )
