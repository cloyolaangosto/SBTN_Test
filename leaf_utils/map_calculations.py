# Methods for map calculations

# MODULES
# GIS
import geopandas as gpd
import logging
import pyproj
import rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from rasterio.enums import Resampling
from rasterio.mask import mask
from rasterio.transform import xy
import rioxarray
import tempfile
from typing import Optional
from joblib import Parallel, delayed

# Data Analysis
import numpy as np
import pandas as pd

# My modules
import leaf_utils.map_plotting as mp

############
### DATA ###
############
# FAO shapefiles paths
country_shp = gpd.read_file("data/CountryLayers/Country_Level0/g2015_2014_0.shp")
subcountry_shp= gpd.read_file("data/CountryLayers/SubCountry_Level1/g2015_2014_1.shp")

# Ecoregions - One map with all ecoregions
er_2017_fp = "data/Ecoregions2017/Ecoregions2017.shp"
er_2017_shp = gpd.read_file(er_2017_fp)


###############
### LOGGERS ###
###############

# Set up the global logger only once.
raster_logger = logging.getLogger("calculate_raster_cf")
if not raster_logger.hasHandlers():
    raster_logger.setLevel(logging.DEBUG)
    # FileHandler logs all details to a single file.
    fh = logging.FileHandler("calculate_cf.log")
    fh.setLevel(logging.DEBUG)
    # StreamHandler prints only key messages (INFO level and above) to the console.
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    raster_logger.addHandler(fh)
    raster_logger.addHandler(ch)

shape_logger = logging.getLogger("calculate_shape_cf")
if not shape_logger.hasHandlers():
    shape_logger.setLevel(logging.DEBUG)
    # FileHandler logs all details to a single file.
    fh = logging.FileHandler("calculate_cf.log")
    fh.setLevel(logging.DEBUG)
    # StreamHandler prints only key messages (INFO level and above) to the console.
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    shape_logger.addHandler(fh)
    shape_logger.addHandler(ch)

#################
### FUNCTIONS ###
#################

# Ecoregion test
def run_categorical_map_plot_ecoregions_2017():
    er_2017_fp = 'data/Ecoregions2017/Ecoregions2017.shp'
    er_2017_shp = gpd.read_file(er_2017_fp)
    er_2017_shp.columns
    er_2017_shp['ECO_NAME'].head()

    mp.plotly_shapefile_categorical(er_2017_shp,'BIOME_NAME', "WWF Ecoregions 2017")
    
#######################
### CF Calculations ###
#######################

# CF Ecoregion calculations
def calculation_ecoregion_cfs_from_shp(cf_shp: gpd.GeoDataFrame, value_column: str, calculated_variable_name = 'area_weighted_cf', er_shp: gpd.GeoDataFrame = er_2017_shp):
    """
    Calculate area-weighted values for each ecoregion in the er_shp based on the values from cf_shp.

    Parameters:
        cf_shp (GeoDataFrame): GeoDataFrame containing polygons and values.
        er_shp (GeoDataFrame): GeoDataFrame containing ecoregions and their shapes.
        value_column (str): The column name in cf_shp that holds the values to be weighted.
        calculated_variable_name (str): Name of the calculated value. Default: area_weighted_cf

    Returns:
        pd.DataFrame: DataFrame with ecoregion names and their corresponding area-weighted values.
        GeoDataFrame: er_shp with a new column containing the area-weighted values.
    """

    # Ensure the CRS of both shapefiles are the same
    if cf_shp.crs != er_shp.crs:
        if er_shp.crs is None:
            raise ValueError("er_shp.crs is None. Cannot reproject cf_shp to a None CRS.")
        cf_shp = cf_shp.to_crs(er_shp.crs)

    # Reprojecting for more accurate calculations. Using a projected CRS (e.g., EPSG:3857 or UTM)
    if cf_shp.crs is not None and getattr(cf_shp.crs, "is_geographic", False):
        cf_shp = cf_shp.to_crs("EPSG:3857")  # Reproject to Web Mercator (or use a UTM zone)

    if er_shp.crs is not None and getattr(er_shp.crs, "is_geographic", False):
        er_shp = er_shp.to_crs("EPSG:3857")  # Reproject ecoregions if needed

    # Perform a spatial join between cf_shp and er_shp based on intersection of polygons in both shapefiles 
    joined_gdf = gpd.sjoin(cf_shp, er_shp, how="inner", predicate="intersects")

    # Calculate the area of each cf_shp polygon
    joined_gdf['area'] = joined_gdf.geometry.area

    # Calculate the weighted value for each polygon (value * area)
    joined_gdf['weighted_value'] = joined_gdf[value_column] * joined_gdf['area']

    # Group by ecoregion and calculate the total area and total weighted value for each ecoregion
    area_weighthed_cf = joined_gdf.groupby('ECO_NAME')['weighted_value'].sum() / joined_gdf.groupby('ECO_NAME')['area'].sum()

    # Create a DataFrame with ecoregions and their corresponding area-weighted values
    result_df = pd.DataFrame(area_weighthed_cf).reset_index()
    result_df.columns = ['ECO_NAME', calculated_variable_name]
    
    # Add the area-weighted values to the er_shp GeoDataFrame
    er_shp = er_shp.merge(result_df, on='ECO_NAME', how='left')

    # Return the result: a DataFrame with ecoregion area-weighted values and the updated er_shp
    return result_df, er_shp

# CF Ecoregion calculations v2
def calculation_ecoregion_cfs_from_shp_v2_merge_OBJECTID(cf_shp: gpd.GeoDataFrame, value_column: str, calculated_variable_name = 'area_weighted_cf', er_shp: gpd.GeoDataFrame = er_2017_shp):
    """
    Calculate area-weighted characterization factors for each ecoregion in the er_shp based on the values from cf_shp. In v2, shapefile results are based on polygon merger. Dataframe still delivers table based on ecoregion calculations. 

    Parameters:
        cf_shp (GeoDataFrame): GeoDataFrame containing polygons and values.
        er_shp (GeoDataFrame): GeoDataFrame containing ecoregions and their shapes.
        value_column (str): The column name in cf_shp that holds the values to be weighted.
        calculated_variable_name (str): Name of the calculated value. Default: area_weighted_cf

    Returns:
        pd.DataFrame: DataFrame with ecoregion names and their corresponding area-weighted cfs.
        GeoDataFrame: er_shp with a new column containing the area-weighted cfs.
    """
    
    # Ensure the CRS of both shapefiles are the same
    if cf_shp.crs != er_shp.crs:
        cf_shp = cf_shp.to_crs(er_shp.crs)

    # Reprojecting for more accurate calculations. Using a projected CRS (e.g., EPSG:3857 or UTM)
    if cf_shp.crs.is_geographic:
        cf_shp = cf_shp.to_crs("EPSG:3857")  # Reproject to Web Mercator (or use a UTM zone)

    if er_shp.crs.is_geographic:
        er_shp = er_shp.to_crs("EPSG:3857")  # Reproject ecoregions if needed

    # Perform a spatial join between cf_shp and er_shp
    joined_gdf = gpd.sjoin(cf_shp, er_shp, how="inner", predicate="intersects")

    # Calculate the area of each cf_shp polygon
    joined_gdf['area'] = joined_gdf.geometry.area

    # Calculate the weighted value for each polygon (value * area)
    joined_gdf['weighted_value'] = joined_gdf[value_column] * joined_gdf['area']

    # Group by ecoregion and calculate the total area and total weighted value for each ecoregion. This is done by ecoregion and polygon, as some ecoregions are divided into several polygons.
    grouped_gdf = joined_gdf.groupby('OBJECTID_right')['weighted_value'].sum() / joined_gdf.groupby('OBJECTID_right')['area'].sum()
    grouped_df = joined_gdf.groupby('ECO_NAME')['weighted_value'].sum() / joined_gdf.groupby('ECO_NAME')['area'].sum()
    
    # Create a DataFrame with ecoregions and their corresponding area-weighted values
    result_df = pd.DataFrame(grouped_df).reset_index()
    result_df.columns = ['ECO_NAME', calculated_variable_name]

    # Add the area-weighted values to the er_shp GeoDataFrame. For this we use the calculations made at the polygon level
    result_df_polygons = pd.DataFrame(grouped_gdf).reset_index()
    result_df_polygons.columns = ['OBJECTID', calculated_variable_name]
    er_shp_result = er_shp.merge(result_df_polygons, on='OBJECTID', how='left')

    # Return the result: a DataFrame with ecoregion area-weighted values and the updated er_shp
    return result_df, er_shp_result


def calculate_area_weighted_cfs_from_shp_with_std_and_median(cf_shp: gpd.GeoDataFrame, value_column: str, brd_shp: gpd.GeoDataFrame = er_2017_shp, brdr_name: str = "ECO_NAME", calculated_variable_name='area_weighted_cf', return_invalid=False, skip_brd_chck=False):
    """
    Calculate area-weighted characterization factors for a given border shapefile (ecoregion global shapefile by default), as well as its standard deviation, and median for each unit. It returns a DataFrame with the results and shapefile with the new characterization factors and related statistics.

    Parameters:
        cf_shp (GeoDataFrame): GeoDataFrame containing polygons and values.
        brd_shp (GeoDataFrame): GeoDataFrame containing ecoregions and their shapes. Default: er_2017_shp (global ecoregion shapefile)
        value_column (str): The column name in cf_shp that holds the values to be weighted.
        calculated_variable_name (str): Name of the calculated value. Default: area_weighted_cf

    Returns:
        pd.DataFrame: DataFrame with ecoregion names and their corresponding area-weighted values, standard deviation, and median.
        GeoDataFrame: brd_shp with new columns containing these statistics.
    """

    print("Calculating area-weighted values for {value_column} with brd_shp shapefile...")
    
    # Checking missing geometries
    print("Checking cf_shp for missing geometries, invalid latitudes...")

    # Drop missing geometries, fixing invalid geometries for cf_shp
    cf_shp = cf_shp[cf_shp.geometry.notna()]
    cf_shp["geometry"] = cf_shp["geometry"].buffer(0)

    # Filtering latitudes for cf_shp
    cf_shp, cf_shp_invalid = filter_invalid_latitudes(cf_shp)    

    # Drop missing geometries, fixing invalid geometries for brd_shp
    if skip_brd_chck:
        print("Skipping brd_shp check.")
    else:
        print("Checking brd_shp for missing geometries, invalid latitudes...")
        brd_shp = brd_shp[brd_shp.geometry.notna()]
        brd_shp["geometry"] = brd_shp["geometry"].buffer(0)

    # Ensure the CRS of both shapefiles are the same
    print("Checking CRS and ensuring they're the same...")
    if cf_shp.crs != brd_shp.crs:
        cf_shp = cf_shp.to_crs(brd_shp.crs)

    # Reprojecting for more accurate calculations using a projected CRS (e.g., EPSG:3857)
    if cf_shp.crs.is_geographic:
        cf_shp = cf_shp.to_crs("EPSG:3857")  

    if brd_shp.crs.is_geographic:
        brd_shp = brd_shp.to_crs("EPSG:3857")

    # Remove NaN and Inf values from the value column
    print("Removing NaN and Inf values for cf_shp...")
    cf_shp = cf_shp[np.isfinite(cf_shp[value_column])]

    print("Data cleaned. Starting calculations...")

    # Perform a spatial join between cf_shp and brd_shp based on intersection of polygons
    print("Performing spatial join...")
    joined_gdf = gpd.sjoin(cf_shp, brd_shp, how="inner", predicate="intersects")  # Spatial join, how = "inner" means only the intersecting polygons are kept, predicate = "intersects" means the polygons intersect

    # Calculate the area of each cf_shp polygon
    joined_gdf['area'] = joined_gdf.geometry.area

    print("Calculating area-weighted values...")
    # Calculate the weighted value for each polygon (value * area)
    joined_gdf['weighted_value'] = joined_gdf[value_column] * joined_gdf['area']

    # Group by ecoregion and calculate total area and total weighted value for each ecoregion
    total_weighted_value = joined_gdf.groupby(brdr_name)['weighted_value'].sum()
    total_area = joined_gdf.groupby(brdr_name)['area'].sum()
    area_weighted_avg = total_weighted_value / total_area

    # Calculate area-weighted standard deviation
    def weighted_std(group):
        mean_value = np.average(group[value_column], weights=group['area'])
        variance = np.average((group[value_column] - mean_value) ** 2, weights=group['area'])
        return np.sqrt(variance)

    area_weighted_std = joined_gdf.groupby(brdr_name).apply(weighted_std)

    # Calculate median
    area_weighted_median = joined_gdf.groupby(brdr_name)[value_column].median()

    print("Calculations completed. Creating results dataframes and shapefiles...")

    # Create a DataFrame with results
    result_df = pd.DataFrame({
        brdr_name: area_weighted_avg.index,
        calculated_variable_name: area_weighted_avg.values,
        f"{calculated_variable_name}_std": area_weighted_std.values,
        f"{calculated_variable_name}_median": area_weighted_median.values
    })

    # Add results to the brd_shp GeoDataFrame
    brd_shp_result = brd_shp.merge(result_df, on=brdr_name, how='left')

    # Return the DataFrame with area-weighted values and the updated brd_shp
    print("All done!")

    if return_invalid:
        return result_df, brd_shp_result, cf_shp_invalid
    else:
        return result_df, brd_shp_result


def calculate_area_weighted_cfs_from_raster_with_std_and_median(raster_input_filepath: str, cf_name: str, cf_unit: str, flow_name: str, area_type: str, run_test = False):
    """
    Calculate area-weighted characterization factors (CFs) from a raster file, including the mean, standard deviation, and median, for specified geographic regions (Ecoregion, Country, or Subcountry).

    Parameters
    -----------
    raster_input_filepath : str
        Filepath to the input raster file containing CF values.
    cf_name : str
        Name of the characterization factor (e.g., "Global Warming Potential").
    cf_unit : str
        Unit of the characterization factor (e.g., "kg CO2-eq").
    flow_name : str
        Name of the flow associated with the CF (e.g., "Carbon Dioxide").
    area_type : str
        Type of geographic area to calculate CFs for. Valid values are "Ecoregion", "Country", or "Subcountry".
    run_test : bool, optional
        If True, the function will only process the first 5 geometries in the shapefile for testing purposes. 
        Default is False.

    Returns
    --------
    final_results : list
        A list containing two elements:
        1. results_df : pandas.DataFrame
            A DataFrame with the calculated CFs (mean, median, and standard deviation) for each region.
        2. final_gdf : geopandas.GeoDataFrame
            A GeoDataFrame with the calculated CFs merged with the original shapefile for spatial data output.
    
    Notes
    ------
    - The function assumes that global variables `er_2017_shp`, `country_shp`, and `subcountry_shp` are defined and contain the shapefiles for Ecoregions, Countries, and Subcountries, respectively.
    - The raster file must have a defined coordinate reference system (CRS) that matches the CRS of the shapefile.
    - The function calculates area-weighted statistics using the pixel dimensions of the raster in the projected CRS.
    - If no valid data is found for a region, that region is skipped.
    
    Raises
    -------
    rioxarray.exceptions.NoDataInBounds
        If no data is found within the bounds of a region during the raster clipping process.
    """
    
    # Use the preconfigured global logger.
    global raster_logger


    # Check if a valid area type is provided
    valid_area_type = {"Ecoregion", "Country", "Subcountry"}
    if area_type is None or area_type not in valid_area_type:
        raster_logger.info("Need to define area type. Valid values are Ecoregion, Country, or Subcountry")
        return

    # Determine which shapefile to use (assume these globals are defined: er_shp, country_shp, subcountry_shp)
    if area_type == "Ecoregion":
        shp = er_2017_shp 
    elif area_type == "Country":
        shp = country_shp
    else:
        shp = subcountry_shp

    #Testing only goes through the first 5 geometries
    if run_test:
        shp = shp.head()
    
    # Opening the raster
    raster = rioxarray.open_rasterio(raster_input_filepath, masked=True)
    raster_crs = raster.rio.crs  # Get the coordinate reference system of the raster.

    raster_logger.info(f"Starting: Calculating {area_type} weighted average CF for {flow_name}")
    
    # Ensure CRS matches
    shp = shp.to_crs(raster_crs)

    # Prepare an empty list for results
    results = []

    #Looping through all shapefile geometries
    for idx, region in shp.iterrows():
        geom = [region["geometry"]]  # Single geometry as list
        
        if area_type == "Ecoregion":
            region_text = "Object # " + str(region["OBJECTID"]) + " - " + region['ECO_NAME']
        elif area_type == "Country":
            region_text = region["ADM0_NAME"]
        else:
            region_text = region["ADM0_NAME"] + " - " + region['ADM1_NAME']

        try:
            # Printing current Ecoregion loop
            raster_logger.debug(f"Calculating  {region_text}")

            # Attempt to mask the raster with the polygon (clip operation)
            masked_raster = raster.rio.clip(geom, drop=True)
            masked_raster_nodata = masked_raster.rio.nodata
            
            # Extract data values
            data = masked_raster.values[0]  # Assuming single band

            if masked_raster_nodata is not None:
                valid_mask = ~np.isnan(data) & (data != masked_raster_nodata)
            else:
                valid_mask = ~np.isnan(data)
            valid_data = data[valid_mask].astype(np.float64)

            # Check if all values are NaN or equal to nodata
            if masked_raster_nodata is not None:
                no_data_condition = np.logical_or(np.isnan(valid_data), valid_data == masked_raster_nodata)
            else:
                no_data_condition = np.isnan(valid_data)
                
            if np.all(no_data_condition):
                raster_logger.debug(f"No overlap found for {region_text}. Skipping...")
                continue

            # Calculate cell areas (using pixel dimensions in projected CRS)
            transform = masked_raster.rio.transform()
            cell_area = abs(transform[0] * transform[4])  # Pixel width * height

            # Compute area-weighted statistics
            if valid_data.size > 0:
                area_weights = np.full_like(a=valid_data, fill_value=cell_area)  # Creates an array of the same shape as valid data with cell area values
                
                # meand and median
                weighted_mean = np.average(valid_data, weights=area_weights)
                weighted_median_value  = weighted_median(valid_data, area_weights)
                
                # Compute weighted variance.
                weighted_variance = np.average((valid_data - weighted_mean)**2,     weights=area_weights)
                weighted_std = np.sqrt(weighted_variance)
            else:
                weighted_mean = np.nan  # No valid data
                weighted_median_value = np.nan  # No valid data
                weighted_std = np.nan  # No valid data

            
            # Append results
            if area_type == "Ecoregion":
                results.append(
                    {
                        "ecoregion_geom_id": region["OBJECTID"],
                        "ecoregion_name": region["ECO_NAME"],
                        "Biome": region["BIOME_NAME"],
                        "impact_category": cf_name,
                        "flow_name": flow_name,
                        "Unit": cf_unit,
                        "cf": weighted_mean,
                        "cf_median": weighted_median_value,
                        "cf_std": weighted_std
                    }
                )
            elif area_type == "Country":
                results.append(
                    {
                        "country": region["ADM0_NAME"],
                        "impact_category": cf_name,
                        "flow_name": flow_name,
                        "Unit": cf_unit,
                        "cf": weighted_mean,
                        "cf_median": weighted_median_value,
                        "cf_std": weighted_std
                    }
                )
            else:
                results.append(
                    {
                        "country": region["ADM0_NAME"],
                        "subcountry (adm1)": region["ADM1_NAME"],
                        "impact_category": cf_name,
                        "flow_name": flow_name,
                        "Unit": cf_unit,
                        "cf": weighted_mean,
                        "cf_median": weighted_median_value,
                        "cf_std": weighted_std
                    }
                )

            raster_logger.debug(f"Calculations finished for region {region_text}. Next!\n")

        except rioxarray.exceptions.NoDataInBounds:
            raster_logger.debug(f"No overalp for region {region_text}, skipping...\n")
            continue  # Skip to the next iteration if clipping fails

    raster_logger.info(f"Calculations complete for {raster_input_filepath}! Found matches for {len(results)} regions.")
    
    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)

    # Merge with shapefile for spatial data output
    if area_type == "Ecoregion":
        final_gdf = shp.merge(results_df, how="left", left_on="OBJECTID", right_on="ecoregion_geom_id")
        final_gdf = final_gdf.drop(columns=["ecoregion_geom_id", "ecoregion_name", "Biome"])
    elif area_type == "Country":
        final_gdf = shp.merge(results_df, how="left", left_on="ADM0_NAME", right_on="country")
        final_gdf = final_gdf.drop(columns=["country"])
    else:
        final_gdf = shp.merge(results_df, how="left", left_on="ADM1_NAME", right_on="subcountry (adm1)")
        final_gdf = final_gdf.drop(columns=["country", "subcountry (adm1)"])
    

    final_results = [results_df, final_gdf]
    return final_results


def weighted_median(values, weights):
    """
    Compute the weighted median of values given corresponding weights.
    
    Parameters:
        values (np.array): Array of data values.
        weights (np.array): Array of weights, same shape as values.
        
    Returns:
        The weighted median.
    """
    # Sort values and weights based on the values.
    sorted_indices = np.argsort(values)  # return indices that would sort an array
    sorted_values = values[sorted_indices]
    sorted_weights = weights[sorted_indices]
    
    # Compute the cumulative sum of weights.
    cumulative_weights = np.cumsum(sorted_weights)
    
    # Find the index where cumulative weight exceeds half the total weight.
    cutoff = cumulative_weights[-1] / 2.0
    median_index = np.where(cumulative_weights >= cutoff)[0][0]
    
    return sorted_values[median_index]


def run_diagnostic(cf_shp, value_column:str, brd_shp: gpd.GeoDataFrame = er_2017_shp):
    # Checking missing geometries
    print(f"Missing geometries in cf_shp: {cf_shp.geometry.isna().sum()}")
    print(f"Missing geometries in brd_shp: {brd_shp.geometry.isna().sum()}")

    # Checking invalid geometries
    print(f"Invalid geometries in cf_shp: {cf_shp[~cf_shp.is_valid].shape[0]} / {cf_shp.shape[0]}")
    print(f"Invalid geometries in brd_shp: {brd_shp[~brd_shp.is_valid].shape[0]} / {brd_shp.shape[0]}")

    # Checking CRS
    print(f"cf_shp CRS: {cf_shp.crs}")
    print(f"brd_shp CRS: {brd_shp.crs}")

    # Checking NaN and Inf values
    cf_shp = cf_shp[np.isfinite(cf_shp[value_column])]

    # Exploring first few geometries
    print(f"cf_shp geometries head: {cf_shp.geometry.head()}")
    print(f"brd_shp geometries head: {brd_shp.geometry.head()}")


def filter_invalid_latitudes(gdf):
    """ 
    Splits a GeoDataFrame into valid and invalid latitude geometries.
    
    Returns:
        valid_gdf (GeoDataFrame): Contains only valid geometries (-90 ≤ latitude ≤ 90).
        invalid_gdf (GeoDataFrame): Contains only invalid geometries (latitude out of bounds).
    """
    def is_valid_latitude(geom):
        minx, miny, maxx, maxy = geom.bounds
        return -90 <= miny <= 90 and -90 <= maxy <= 90  # Check if within valid latitudes
    
    # Filter valid and invalid geometries
    valid_gdf = gdf[gdf.geometry.apply(is_valid_latitude)]
    invalid_gdf = gdf[~gdf.geometry.apply(is_valid_latitude)]

    print(f"Valid geometries: {valid_gdf.shape[0]}")
    print(f"Invalid geometries: {invalid_gdf.shape[0]}")
    
    return valid_gdf, invalid_gdf


#################
def process_region(region, raster_input_filepath, cf_name, cf_unit, flow_name, area_type, raster_crs, logger):
    """
    Process a single region: clip the raster to its geometry and compute area-weighted statistics.
    Returns a dictionary with the results (or None if the region is skipped).
    """
    # Each worker re-opens the raster file.
    try:
        raster = rioxarray.open_rasterio(raster_input_filepath)
    except Exception as e:
        logger.error(f"Error opening raster in worker: {e}")
        return None

    geom = [region["geometry"]]
    if area_type == "Ecoregion":
        region_text = "Object # " + str(region["OBJECTID"]) + " - " + region['ECO_NAME']
    elif area_type == "Country":
        region_text = region["ADM0_NAME"]
    else:
        region_text = region["ADM0_NAME"] + " - " + region['ADM1_NAME']

    try:
        logger.debug(f"Processing region: {region_text}")
        # Clip the raster using rioxarray's clip function.
        masked_raster = raster.rio.clip(geom, drop=True)
        masked_raster_nodata = masked_raster.rio.nodata

        data = masked_raster.values[0]  # Assuming single band

        # Create a valid data mask.
        valid_mask = ~np.isnan(data)
        if masked_raster_nodata is not None:
            valid_mask &= (data != masked_raster_nodata)
        valid_data = data[valid_mask].astype(np.float64)

        # If all data is nodata, skip this region.
        if masked_raster_nodata is not None:
            no_data_condition = np.logical_or(np.isnan(valid_data), valid_data == masked_raster_nodata)
        else:
            no_data_condition = np.isnan(valid_data)
        if np.all(no_data_condition):
            logger.debug(f"No valid data for {region_text}. Skipping...")
            return None

        # Compute cell area from the affine transform.
        transform = masked_raster.rio.transform()
        cell_area = abs(transform[0] * transform[4])

        if valid_data.size > 0:
            area_weights = np.full_like(valid_data, fill_value=cell_area)
            weighted_mean = np.average(valid_data, weights=area_weights)
            weighted_median_value = weighted_median(valid_data, area_weights)
            weighted_variance = np.average((valid_data - weighted_mean) ** 2, weights=area_weights)
            weighted_std = np.sqrt(weighted_variance)
        else:
            weighted_mean = weighted_median_value = weighted_std = np.nan

        # Build the result dictionary.
        if area_type == "Ecoregion":
            result = {
                "ecoregion_geom_id": region["OBJECTID"],
                "ecoregion_name": region["ECO_NAME"],
                "Biome": region["BIOME_NAME"],
                "impact_category": cf_name,
                "flow_name": flow_name,
                "Unit": cf_unit,
                "cf": weighted_mean,
                "cf_median": weighted_median_value,
                "cf_std": weighted_std
            }
        elif area_type == "Country":
            result = {
                "country": region["ADM0_NAME"],
                "impact_category": cf_name,
                "flow_name": flow_name,
                "Unit": cf_unit,
                "cf": weighted_mean,
                "cf_median": weighted_median_value,
                "cf_std": weighted_std
            }
        else:
            result = {
                "country": region["ADM0_NAME"],
                "subcountry (adm1)": region["ADM1_NAME"],
                "impact_category": cf_name,
                "flow_name": flow_name,
                "Unit": cf_unit,
                "cf": weighted_mean,
                "cf_median": weighted_median_value,
                "cf_std": weighted_std
            }
        logger.debug(f"Finished processing {region_text}.")
        return result

    except rioxarray.exceptions.NoDataInBounds:
        logger.debug(f"No overlap for region {region_text}, skipping...")
        return None
    except Exception as e:
        logger.error(f"Error processing region {region_text}: {e}")
        return None

def calculate_area_weighted_cfs_from_raster_with_std_and_median_MULTICORE(
    raster_input_filepath: str, 
    cf_name: str, 
    cf_unit: str, 
    flow_name: str, 
    area_type: str, 
    run_test=False,
    n_jobs=-1  # Use all available cores by default.
):
    """
    Optimized function that parallelizes per-region processing.
    Detailed logs are recorded via the global logger; only key messages appear on the console.
    """
    global raster_logger  # Assuming a global logger is configured

    # Validate area type.
    valid_area_type = {"Ecoregion", "Country", "Subcountry"}
    if area_type is None or area_type not in valid_area_type:
        raster_logger.info("Need to define area type. Valid values are Ecoregion, Country, or Subcountry")
        return

    # Select the appropriate shapefile.
    if area_type == "Ecoregion":
        shp = er_2017_shp
    elif area_type == "Country":
        shp = country_shp
    else:
        shp = subcountry_shp

    if run_test:
        shp = shp.head()

    # Open the raster once to get its CRS; then close it.
    raster = rioxarray.open_rasterio(raster_input_filepath)
    raster_crs = raster.rio.crs
    raster.close()
    shp = shp.to_crs(raster_crs)

    raster_logger.info(f"Starting: Calculating {area_type} weighted average CF for {flow_name} with MULTICORE")

    # Process all regions in parallel.
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_region)(
            region, raster_input_filepath, cf_name, cf_unit, flow_name, area_type, raster_crs, raster_logger
        )
        for idx, region in shp.iterrows()
    )
    # Remove None results (skipped regions)
    results = [r for r in results if r is not None]

    raster_logger.info(f"Finished: Calculations complete! Found matches for {len(results)} regions.")

    results_df = pd.DataFrame(results)
    if area_type == "Ecoregion":
        final_gdf = shp.merge(results_df, how="left", left_on="OBJECTID", right_on="ecoregion_geom_id")
        final_gdf = final_gdf.drop(columns=["ecoregion_geom_id", "ecoregion_name", "Biome"])
    elif area_type == "Country":
        final_gdf = shp.merge(results_df, how="left", left_on="ADM0_NAME", right_on="country")
        final_gdf = final_gdf.drop(columns=["country"])
    else:
        final_gdf = shp.merge(results_df, how="left", left_on="ADM1_NAME", right_on="subcountry (adm1)")
        final_gdf = final_gdf.drop(columns=["country", "subcountry (adm1)"])

    final_results = [results_df, final_gdf]
    return final_results




############################
### SHAPEFILE OPERAIOTNS ###
############################

def rasterize_shapefile_to_target_raster(gdf: gpd.geodataframe, raster_filepath: str, value_column: str, output_path= None, band = 1, no_data = None):
    
    # Checks the value column exist
    if value_column not in gdf.columns:
        raise ValueError(f"No {value_column} in gdf")
    
    print("Opening target raster")
    # Open the high-resolution raster
    with rasterio.open(raster_filepath) as src:
        raster = src.read(band)  # Assumes it's the first band
        meta = src.meta.copy()         # Copy metadata for potential output.
        transform = src.transform
        raster_shape = raster.shape
        no_data_input = src.nodata

    # Ensure the shapefile's CRS matches the raster's CRS
    print("Projecting shapefile to target raster")
    gdf = gdf.to_crs(src.crs)

    # Rasterize the polygon values to match the raster resolution
    # Each polygon is assigned its "value" attribute.
    print("Rasterizing shapefile")
    polygon_raster = rasterize(
        [(geom, value) for geom, value in zip(gdf.geometry, gdf[value_column])],
        out_shape=raster_shape,
        transform=transform,
        fill= no_data if no_data is not None else no_data_input,           # Pixels not covered by any polygon get a fill with no_data_value
        all_touched=False         # Use all_touched=True to account for partial pixel coverage
    )

    # Saves the raster
    # If an output path is provided, write the product to a new file.
    if output_path:
        print("Saving raster")
        meta.update(dtype=polygon_raster.dtype, nodata=no_data if no_data is not None else no_data_input)
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(polygon_raster, band)
    
    else:
        return polygon_raster


def calculate_polygon_means_from_raster(raster_path: str,
                                        polygon_gdf: gpd.GeoDataFrame,
                                        id_column: str,
                                        band: int = 1,
                                        return_gdf: bool = False,
                                        col_mean_name: str = "mean_value"):
    """
    Overlaps a polygon shapefile with a raster and computes the mean raster value for each polygon.

    Parameters:
    -----------
    raster_path : str
        Path to the input raster (GeoTIFF).
    polygon_gdf : GeoDataFrame
        GeoDataFrame of polygons to overlap with the raster.
    id_column : str
        Unique column name used to identify each polygon (e.g., 'ECO_NAME').
    band : int
        Raster band to read (default is 1).
    return_gdf : bool
        If True, returns a copy of the input GeoDataFrame with a new column for mean values.

    Returns:
    --------
    results_df : pd.DataFrame
        DataFrame with [id_column, 'mean_raster_value'].
    result_gdf (optional): GeoDataFrame
        If return_gdf=True, returns polygon_gdf with the mean values added.
    """
    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
        polygon_gdf = polygon_gdf.to_crs(raster_crs)

        results = []
        for idx, row in polygon_gdf.iterrows():
            geom = [row.geometry]
            try:
                out_image, out_transform = mask(src, geom, crop=True)
                data = out_image[band - 1]
                nodata = src.nodata

                if nodata is not None:
                    valid = data[data != nodata]
                else:
                    valid = data[np.isfinite(data)]

                if valid.size > 0:
                    mean_val = valid.mean()
                else:
                    mean_val = np.nan

                results.append((row[id_column], mean_val))

            except Exception:
                results.append((row[id_column], np.nan))

    results_df = pd.DataFrame(results, columns=[id_column, col_mean_name])

    if return_gdf:
        merged = polygon_gdf.copy()
        merged = merged.merge(results_df, on=id_column, how="left")
        return results_df, merged
    else:
        return results_df

#########################
### RASTER OPERAIOTNS ###
#########################

def multiply_rasters(raster_paths, output_path=None, band=1):
    """
    Multiplies n number of rasters element-wise after verifying that they all share the same 
    coordinate system (CRS) and dimensions, while handling no_data values.
    
    Handling of no_data:
    - For each raster, if a no_data value is defined, we create/update a 'valid_mask' that is True 
      where pixel values are valid (i.e., not equal to no_data) and False where they are not.
    - During multiplication, no_data pixels are temporarily replaced with 1 (the neutral element for multiplication)
      so that they do not affect the product of valid pixels.
    - After all rasters have been multiplied, any pixel that was flagged as invalid (i.e. no_data in any raster)
      is set to the output no_data value (taken from the first raster if available, or np.nan otherwise).
    
    Parameters:
        raster_paths (list of str): List of file paths for the rasters to multiply.
        output_path (str, optional): Path to write the resulting product raster.
        band (int, optional): The band number to use from each raster (default is 1).
        
    Returns:
        product (numpy.ndarray): The resulting product array from multiplying all rasters.
        
    Raises:
        ValueError: If no raster paths are provided, if the rasters have different dimensions,
                    or if they have different CRS.
    """
    # Check if any raster file paths were provided
    if not raster_paths:
        raise ValueError("No raster file paths provided.")
    
    # Open the first raster to initialize the product, metadata, CRS, and no_data value.
    with rasterio.open(raster_paths[0]) as src:
        # Read the specified band and convert to float64 for more precise multiplication.
        data0 = src.read(band).astype("float64")
        meta = src.meta.copy()         # Copy metadata for potential output.
        base_crs = src.crs             # Get the coordinate reference system of the first raster.
        nodata0 = src.nodata           # Get the no_data value for the first raster.
    
    # Create a valid mask that tracks where the data are valid (i.e. not no_data).
    if nodata0 is not None:
        # True where the value is not equal to the no_data value.
        valid_mask = (data0 != nodata0)
        # Replace no_data values with 1 so that they do not affect the multiplication.
        data0 = np.where(data0 == nodata0, 1, data0)
    else:
        valid_mask = np.ones_like(data0, dtype=bool)  # All values are valid if no no_data is defined.
    
    # Initialize the product with the first raster's data.
    product = data0.copy()
    
    # Loop through the remaining raster files.
    for path in raster_paths[1:]:
        with rasterio.open(path) as src:
            # Verify that the CRS is the same as the first raster.
            if src.crs != base_crs:
                raise ValueError(f"CRS mismatch: {raster_paths[0]} and {path} have different coordinate systems.")
            
            # Read the current raster band as float64.
            data = src.read(band).astype("float64")
            
            # Ensure the current raster has the same dimensions.
            if data.shape != product.shape:
                raise ValueError(f"Raster shapes do not match: {raster_paths[0]} and {path}")
            
            # Get the no_data value for the current raster.
            nodata = src.nodata
            if nodata is not None:
                # Update the valid mask: a pixel remains valid only if it's valid in all rasters.
                valid_mask &= (data != nodata)
                # Replace no_data values with 1 (neutral for multiplication).
                data = np.where(data == nodata, 1, data)
            else:
                # If no no_data is defined, consider all pixels as valid.
                valid_mask &= np.ones_like(data, dtype=bool)
            
            # Multiply the current data into the product.
            product *= data
    
    # Determine the output no_data value:
    # Use the first raster's no_data if available; otherwise, set to np.nan.
    output_nodata = nodata0 if nodata0 is not None else np.nan
    # Set any pixel that was invalid (i.e., no_data in any input) to the output no_data value.
    product[~valid_mask] = output_nodata
    
    # If an output path is provided, write the product to a new file.
    if output_path:
        meta.update(dtype=product.dtype, nodata=output_nodata)
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(product, band)
    
    return product


def create_binary_mask(input_path, output_path, binary_value = 1, band=1, src_nodata=None, dst_nodata = -3.14e-8):
    """
    Reads an input GeoTIFF and creates a new GeoTIFF file where each cell is:
      - 1 if the cell contains a valid value
      - 0 if the cell is considered no_data
      
    You can optionally manually specify a no_data value (manual_nodata). This manual value takes
    precedence over any no_data value found in the input file. It can also be np.nan.

    Handling of no_data:
      - If a manual no_data value is provided:
          - If it's a numeric value, cells equal to that value are considered no_data.
          - If it's np.nan, cells with NaN values (checked via np.isnan) are considered no_data.
      - If no manual value is provided, the function uses the input file's nodata attribute.
      - If neither is provided, the function assumes that a cell "has value" if its value is non-zero.

    Parameters:
        input_path (str): Path to the input GeoTIFF file.
        output_path (str): Path to write the binary mask GeoTIFF file.
        band (int): The band number to process (default is 1).
        src_nodata (optional): Manually specified input, or source, no_data value (can be a number or np.nan). Default is None.
        dst_nodata (float): Destination or output no data value. Default is 3.14E-8
        
    Returns:
        binary (numpy.ndarray): The 2D binary mask array.
    """
    # Open the input GeoTIFF file.
    with rasterio.open(input_path) as src:
        # Read the specified band from the raster.
        data = src.read(band)
        # Use the manual_nodata if provided; otherwise use the file's nodata.
        nodata = src_nodata if src_nodata else src.nodata
        # Copy metadata from the source for writing the new GeoTIFF.
        meta = src.meta.copy()
    
    # Create the binary mask.
    # If a nodata value is defined, determine valid pixels.
    if nodata:
        # If the provided nodata value is np.nan, use np.isnan to check for NaN values.
        if np.isnan(nodata):
            # Valid cells are those that are not NaN.
            print(f"Input No Data defined as NaN. Value {binary_value} if is not NaN and {dst_nodata} if it is.")
            valid = ~np.isnan(data)
        else:
            print(f"Input No Data defined as {nodata}. Value {binary_value} if is not {nodata} and {dst_nodata} if it is.")
            # Valid cells are those that do not equal the nodata value.
            valid = (data != nodata)
        
        # Create a binary mask: 1 for valid, 0 for no_data.
        binary = np.where(valid, binary_value, dst_nodata).astype(np.float32)

    else:
        print(f"No Input nodata value exists. Value {binary_value} if cell has value and {dst_nodata} if not, disregarding which value it is.")
        # If no nodata value is defined, assume a cell has a value if it is non-zero.
        binary = np.where(data != 0, binary_value, dst_nodata).astype(np.float32)
    
    # Update metadata for output:
    # - Set the data type to uint8 (since our mask only contains 0 and 1).
    # - Ensure the number of bands is 1.
    # - Set the nodata value in the output metadata.
    meta.update({
        'dtype': 'float32',
        'count': 1,
        'nodata': dst_nodata
    })
    
    # Write the binary mask to a new GeoTIFF file.
    with rasterio.open(output_path, 'w', **meta) as dst:
        dst.write(binary, 1)
    
    return binary


def resample_to_match(src_path, target_path, output_path,
                      resampling_method=Resampling.nearest,
                      src_nodata=None, dst_nodata=None):
    """
    Resamples the source raster (src_path) so that it matches the resolution,
    extent, and geotransform of the target raster (target_path). During resampling,
    specified no_data values are excluded from interpolation.

    Parameters:
        src_path (str): File path of the source raster that needs to be resampled.
        target_path (str): File path of the target raster whose grid we want the source to match.
        output_path (str): File path to save the resampled source raster.
        resampling_method: Rasterio resampling method (e.g., Resampling.nearest, Resampling.bilinear).
                           Default is Resampling.nearest.
        src_nodata (optional): Value to treat as no_data in the source raster.
                               If not provided, the function uses the source raster's nodata.
        dst_nodata (optional): Value to use as no_data in the destination raster.
                               If not provided, the function uses src_nodata.
    
    Returns:
        None. The resampled raster is written to output_path.
    
    Explanation:
        - The function reads the target raster’s transform, width, height, and CRS.
        - It then opens the source raster and reads its data.
        - If src_nodata is not provided, it is taken from the source raster's metadata.
        - Similarly, if dst_nodata is not provided, it is set equal to src_nodata.
        - The reproject function is then used to resample the source data onto the target grid,
          while ignoring (excluding) source pixels that equal src_nodata.
        - The resulting output is written with the specified destination nodata value.
    """
    # Open the target raster to get its grid parameters.
    with rasterio.open(target_path) as target:
        target_transform = target.transform
        target_width = target.width
        target_height = target.height
        target_crs = target.crs

    # Open the source raster.
    with rasterio.open(src_path) as src:
        # Read the first band of the source raster.
        src_data = src.read(1)
        src_transform = src.transform
        src_crs = src.crs

        # Determine the source no_data value: use provided or the one in the file.
        if src_nodata is None:
            src_nodata = src.nodata

        # If destination no_data is not provided, set it equal to src_nodata.
        if dst_nodata is None:
            dst_nodata = src_nodata

        # Prepare an empty array for the resampled data.
        resampled_data = np.empty((target_height, target_width), dtype=src_data.dtype)

        # Perform the resampling with reproject.
        reproject(
            source=src_data,
            destination=resampled_data,
            src_transform=src_transform,
            src_crs=src_crs,
            dst_transform=target_transform,
            dst_crs=target_crs,
            resampling=resampling_method,
            src_nodata=src_nodata,  # Exclude these values from interpolation
            dst_nodata=dst_nodata
        )

        # Copy and update metadata for the output raster.
        meta = src.meta.copy()
        meta.update({
            'crs': target_crs,
            'transform': target_transform,
            'width': target_width,
            'height': target_height,
            'nodata': dst_nodata
        })

    # Write the resampled raster to disk.
    with rasterio.open(output_path, 'w', **meta) as dst:
        dst.write(resampled_data, 1)


def resample_to_match_multiband(src_path: str,
                                target_path: str,
                                output_path: str,
                                resampling_method=Resampling.nearest):
    """
    Resample every band in src_path to the grid of target_path and
    write out a single multi‐band GeoTIFF at output_path.
    """
    # --- 1) Read & stash metadata ---
    with rasterio.open(src_path) as src, rasterio.open(target_path) as tgt:
        # source parameters
        src_count     = src.count
        src_nodata    = src.nodata
        src_dtype     = src.dtypes[0]
        src_crs       = src.crs
        src_transform = src.transform

        # target parameters & base profile
        tgt_crs       = tgt.crs
        tgt_transform = tgt.transform
        tgt_width     = tgt.width
        tgt_height    = tgt.height
        profile       = tgt.meta.copy()

    # update profile to match source bands & nodata
    profile.update({
        "count":  src_count,
        "nodata": src_nodata,
        "dtype":  src_dtype,
    })

    # --- 2) Re‐open src for reading and dst for writing ---
    with rasterio.open(src_path) as src, rasterio.open(output_path, "w", **profile) as dst:
        for b in range(1, src_count + 1):
            src_band = src.read(b)
            out_band = np.empty((tgt_height, tgt_width), dtype=src_dtype)

            reproject(
                source=src_band,
                destination=out_band,
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=tgt_transform,
                dst_crs=tgt_crs,
                resampling=resampling_method,
                src_nodata=src_nodata,
                dst_nodata=src_nodata
            )

            dst.write(out_band, b)


def resample_to_match_noSaving(
    src_path: str,
    target_path: str,
    resampling_method=Resampling.nearest,
    dst_nodata: Optional[float] = None
) -> np.ndarray:
    """
    Resample `src_path` onto the grid of `target_path` and return the array.

    Parameters:
        src_path (str):        Source raster file.
        target_path (str):     Reference raster whose grid (CRS, transform, shape) to match.
        resampling_method:     One of rasterio.enums.Resampling.* (default: nearest).
        dst_nodata (float):    Nodata value to write into the resampled array;
                               if None, uses the source raster’s nodata.

    Returns:
        np.ndarray: 2D array of resampled data, dtype float32, with `dst_nodata` where no data.
    """


    # Load the target grid
    with rasterio.open(target_path) as tgt:
        tgt_transform = tgt.transform
        tgt_width     = tgt.width
        tgt_height    = tgt.height
        tgt_crs       = tgt.crs

    # Load the source band
    with rasterio.open(src_path) as src:
        src_data     = src.read(1).astype("float32")
        src_transform= src.transform
        src_crs      = src.crs
        src_nodata   = src.nodata
        dst_nodata   = dst_nodata if dst_nodata is not None else src_nodata

    # Prepare destination array as float32
    resampled = np.full((tgt_height, tgt_width), dst_nodata, dtype="float32")

    # Run reprojection/resampling
    reproject(
        source=src_data,
        destination=resampled,
        src_transform=src_transform,
        src_crs=src_crs,
        src_nodata=src_nodata,
        dst_transform=tgt_transform,
        dst_crs=tgt_crs,
        dst_nodata=dst_nodata,
        resampling=resampling_method
    )

    return resampled


def subtract_rasters_union(raster1_path, raster2_path, output_path, band=1, output_nodata=None, resampling_method=Resampling.nearest):
    """
    Subtracts two rasters over a common grid even when their overlapping zones are not contiguous.
    
    Steps:
      1. Opens both rasters and computes the union (combined extent) of their footprints.
      2. Uses the resolution (pixel size) from the first raster (you could choose differently)
         to define a common output grid covering the union of both extents.
      3. Reprojects each raster onto this common grid. During reprojection, the given (or file’s)
         no_data values are used so that only valid data is reprojected.
      4. Computes cell-by-cell subtraction where both reprojected arrays have valid data.
         In cells where one or both have no_data, the output is set to output_nodata.
      5. Writes the resulting raster to disk with updated metadata.
    
    Parameters:
      raster1_path (str): File path to the first raster.
      raster2_path (str): File path to the second raster.
      output_path (str): File path for the output raster.
      band (int, optional): Band number to process (default is 1).
      output_nodata (optional): User-defined nodata value for the output.
          If None, the function uses the first available nodata from raster1, then raster2, else -9999.
      resampling_method (optional): Rasterio resampling method (default is Resampling.nearest).
    
    Returns:
      result (numpy.ndarray): The resulting difference array on the common grid.
    """
    # Open both rasters to retrieve bounds and other info.
    with rasterio.open(raster1_path) as src1, rasterio.open(raster2_path) as src2:
        # Compute union of the extents (footprints)
        union_left   = min(src1.bounds.left,   src2.bounds.left)
        union_bottom = min(src1.bounds.bottom, src2.bounds.bottom)
        union_right  = max(src1.bounds.right,  src2.bounds.right)
        union_top    = max(src1.bounds.top,    src2.bounds.top)
        
        # Decide on the resolution. Here we use the resolution of the first raster.
        res_x, res_y = src1.res
        
        # Calculate output dimensions: width and height for the union extent.
        union_width  = int(np.ceil((union_right - union_left) / res_x))
        union_height = int(np.ceil((union_top - union_bottom) / res_y))
        
        # Define the transform for the union grid.
        # from_origin expects (upper left x, upper left y, x resolution, y resolution)
        union_transform = from_origin(union_left, union_top, res_x, res_y)
        
        # Determine no_data values:
        nodata1 = src1.nodata
        nodata2 = src2.nodata
        if output_nodata is None:
            if nodata1 is not None:
                output_nodata = nodata1
            elif nodata2 is not None:
                output_nodata = nodata2
            else:
                output_nodata = -3.14e-8
        
        # Create empty arrays to hold the reprojected data.
        # We use the data type from the first raster.
        data1_union = np.full((union_height, union_width), output_nodata, dtype=src1.dtypes[0])
        data2_union = np.full((union_height, union_width), output_nodata, dtype=src2.dtypes[0])
        
        # Reproject the first raster onto the union grid.
        reproject(
            source=src1.read(band),
            destination=data1_union,
            src_transform=src1.transform,
            src_crs=src1.crs,
            dst_transform=union_transform,
            dst_crs=src1.crs,
            resampling=resampling_method,
            src_nodata=nodata1,
            dst_nodata=output_nodata
        )
        
        # Reproject the second raster onto the same union grid.
        reproject(
            source=src2.read(band),
            destination=data2_union,
            src_transform=src2.transform,
            src_crs=src2.crs,
            dst_transform=union_transform,
            dst_crs=src2.crs,
            resampling=resampling_method,
            src_nodata=nodata2,
            dst_nodata=output_nodata
        )
    
    # Create a mask: valid only where both arrays are not the output nodata.
    valid_mask = (data1_union != output_nodata) & (data2_union != output_nodata)
    
    # Perform cell-by-cell subtraction where both are valid.
    # Else, assign the output nodata value.
    result = np.where(valid_mask, data1_union - data2_union, output_nodata)
    
    # Update metadata for the output file using src1's metadata as a base.
    meta = src1.meta.copy()
    meta.update({
        "driver": "GTiff",
        "height": union_height,
        "width": union_width,
        "transform": union_transform,
        "nodata": output_nodata
    })
    
    # Write the result to disk.
    with rasterio.open(output_path, "w", **meta) as dst:
        dst.write(result, band)
    
    return result


def mask_raster1_by_overlap_with_raster2(raster1_path: str, raster2_path: str, output_path: str = None):
    """
    Aligns raster2 to raster1 (if needed), and returns a raster array from raster1
    masked by the valid data of raster2. Optionally saves to disk.

    Parameters:
    -----------
    raster1_path : str
        Path to the source raster from which values are kept.
    raster2_path : str
        Path to the raster used to define overlap (non-nodata mask).
    output_path : str, optional
        If provided, saves the masked raster to this file.

    Returns:
    --------
    masked_array : np.ndarray
        Array with values from raster1 where raster2 is valid; np.nan elsewhere.
    """

    with rasterio.open(raster1_path) as src1:
        meta1 = src1.meta.copy()
        data1 = src1.read(1)

        with rasterio.open(raster2_path) as src2:
            # If CRS, transform, or shape mismatch, resample raster2 to match raster1
            if src1.crs != src2.crs or src1.transform != src2.transform or src1.shape != src2.shape:
                tmp_resampled = tempfile.NamedTemporaryFile(suffix=".tif", delete=False).name
                resample_to_match(
                    src_path=raster2_path,
                    target_path=raster1_path,
                    output_path=tmp_resampled,
                    resampling_method=Resampling.nearest
                )
                with rasterio.open(tmp_resampled) as r2_resampled:
                    data2 = r2_resampled.read(1)
                    nodata2 = r2_resampled.nodata
            else:
                data2 = src2.read(1)
                nodata2 = src2.nodata

    # Mask logic: retain where raster2 is valid
    mask = np.isfinite(data2) if nodata2 is None else data2 != nodata2
    masked_data = np.where(mask, data1, np.nan)

    # Write output if requested
    if output_path:
        meta1.update(dtype='float32', nodata=np.nan)
        with rasterio.open(output_path, 'w', **meta1) as dst:
            dst.write(masked_data.astype('float32'), 1)

    return masked_data



def get_pixel_center_latitude(raster_path: str, row: int, col: int) -> float:
    """
    Return the latitude (in degrees) of the centre of the pixel at (row, col)
    in the given raster.

    Parameters
    ----------
    raster_path : str
        Path to a GeoTIFF (or other rasterio-readable file), assumed to have
        a georeference (transform + CRS).
    row : int
        Row index of the pixel (0‐based).
    col : int
        Column index of the pixel (0‐based).

    Returns
    -------
    float
        Latitude of the pixel centre, in EPSG:4326.
    """
    with rasterio.open(raster_path) as src:
        # Compute the centre x,y in the raster's CRS
        x_center, y_center = xy(src.transform, row, col, offset='center')

        # If the raster is already in lat/lon, just return y_center
        if src.crs.to_string() == 'EPSG:4326':
            return y_center

        # Otherwise, reproject to EPSG:4326
        transformer = pyproj.Transformer.from_crs(src.crs, 'EPSG:4326', always_xy=True)
        lon, lat = transformer.transform(x_center, y_center)
        return lat


def binarize_raster(
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


#################################
### GET DATA POINTS FUNCTIONS ###
#################################

import geopandas as gpd

def extract_coordinates_values_from_shapefile(points_source, polygon_shapefile, column_names, points_crs="EPSG:4326"):
    """
    Extract values from a polygon shapefile for each point provided.

    Parameters
    -----------
    points_source : str or GeoDataFrame
        Either the path to a points shapefile or an existing GeoDataFrame containing point geometries.
    polygon_shapefile_path : str
        Shapefile from which values should be extracted.
    column_names : list of str
        List of column names from the polygon shapefile to extract.
    points_crs : str, optional
        The coordinate reference system for the points if they don't have one.
        Default is "EPSG:4326" (WGS84).

    Returns
    --------
    GeoDataFrame
        A GeoDataFrame of the points with new columns appended from the polygon shapefile.
        Points that do not fall within any polygon will have NaN for those attribute columns.
        
    Notes
    ------
    - If a point falls within more than one polygon, the spatial join will return duplicate rows.
    - The function uses the "within" predicate to match points to polygons.
    """
    # Load the points data. If a file path is provided, read it; otherwise assume it's a GeoDataFrame.
    if isinstance(points_source, str):
        points_gdf = gpd.read_file(points_source)
    else:
        points_gdf = points_source.copy()

    # If the points do not have a CRS, assign the default points_crs.
    if points_gdf.crs is None:
        points_gdf.set_crs(points_crs, inplace=True)

    # Load the polygon shapefile.
    poly_gdf = polygon_shapefile.copy()

    # Reproject the polygon shapefile to the points CRS if necessary.
    if poly_gdf.crs != points_gdf.crs:
        poly_gdf = poly_gdf.to_crs(points_gdf.crs)

    # Select only the columns we need (plus the geometry column) from the polygon data.
    poly_subset = poly_gdf[column_names + ["geometry"]]

    # Perform a spatial join: for each point, find the polygon in which it falls.
    # Using predicate="within" so points that fall within a polygon get joined.
    joined_gdf = gpd.sjoin(points_gdf, poly_subset, how="left", predicate="within")

    # Optionally, drop the index_right column generated by the spatial join.
    if "index_right" in joined_gdf.columns:
        joined_gdf = joined_gdf.drop(columns=["index_right"])

    return joined_gdf


def extract_coordinates_values_from_raster(points_source, raster_input_filepath, band=1, column_name="raster_value"):
    """
    Extract raster values for each point from the raster and append them as a new column.
    
    Parameters
    -----------
    points_source : str or GeoDataFrame
        Either the path to a point shapefile or an existing GeoDataFrame with point geometries.
    raster_input_filepath : str
        File path to the raster file.
    band : int, optional
        The raster band to extract values from (default is 1).
    column_name : str, optional
        The name of the new column to hold the raster values (default is "raster_value").
        
    Returns
    --------
    GeoDataFrame
        A GeoDataFrame with the extracted raster values appended as a new column.
        
    Notes
    -----
    - If a point falls outside the raster or corresponds to the nodata value, it will be assigned NaN.
    """
    # Load points data from file if a file path is provided.
    if isinstance(points_source, str):
        points_gdf = gpd.read_file(points_source)
    else:
        points_gdf = points_source.copy()
    
    # Open the raster file
    with rasterio.open(raster_input_filepath) as src:
        raster_crs = src.crs
        
        # Reproject points to raster CRS if necessary.
        if points_gdf.crs != raster_crs:
            points_gdf = points_gdf.to_crs(raster_crs)
        
        # Prepare a list of (x, y) tuples from each point's geometry.
        coords = [(geom.x, geom.y) for geom in points_gdf.geometry]
        
        # Retrieve the nodata value for the raster.
        nodata = src.nodata
        
        # Sample the raster at the point locations.
        # src.sample() returns an iterator of arrays; we extract the first element of each array.
        sampled = list(src.sample(coords, indexes=band))
        values = [np.nan if (nodata is not None and arr[0] == nodata) else arr[0] for arr in sampled]
    
    # Append the extracted values as a new column in the GeoDataFrame.
    points_gdf[column_name] = values
    return points_gdf