# Crop Water Requirement and Evapotranspiration

"""This script calculates the crop water requirement depending on the crop and climate it's located. It's calculated as the simple multiplifcation of potential evapotranspiration (PET) and a crop coefficient depending on it's growing stage (K). PET is calculated using the Thornthwaite equation (originally from Thornthwaite, C. W. An Approach toward a Rational Classification of Climate. Geogr. Rev. 1948, 38, 55, doi:10.2307/210739.), while K is obtained from Chapagain, A. K.; Hoekstra, A. Y. Water footprint of nations. Volume 1 : Main report; 2004."""

## PET - Thornwaite equation
"""This equation calculates PET based on duration of sunlight in hours, varying with season and latitude (L), number of days in a month (N), average monttly air temperature (T, in °C), and heat index (I_a). Detailed explanations can be found here: https://wikifire.wsl.ch/tiki-indexf125.html?page=Potential+evapotranspiration

PET is calculated as:
PET = 0 if T <0,
PET = 1.6 * (L/12) * (N/30) * (10 * T/I)**a

Where a is calcualted as:
a = (6.75 * 10^-7⋅ * I^3)−(7.71*10^−5 * I^2) + (0.01792 * I)+(0.49239)

There has been more updates since it's original publication since 1948,  but to keep consistency with the SOC publication we're replicating, we are using it in it's original format."""


"""
## Crop Coefficient (K)
Crop coefficient is, of course, crop dependent, as well as weather dependent. It is used to adjust the PET_0 that is only weather dependent to crops, depending on their growing stage. It follows a curve, which for annual crops looks like this:

For forage crops, it's different, as they have several rotations throughout a year, and looks like this:

For fruit trees and trees plantation... I need to read more.

Complete documentation can be found here: https://www.fao.org/4/X0490E/x0490e00.htm

The following part of the code tries to automatize this calculations for annual crops, as a first example.
"""

#### Modules ####
import math  # For mathematical operations
import calendar  # For month and day calculations
from calendar import monthrange  # For getting the number of days in a month
from pathlib import Path
from typing import Optional, Tuple, Union

import polars as pl  # For data manipulation and analysis
import numpy as np  # For numerical operations and array handling
import rasterio  # For raster data handling
from rasterio.warp import reproject, Resampling  # For raster reprojection and resampling
from collections import defaultdict  # Imports a special dict subclass that auto‐initializes missing keys.
from tqdm import tqdm  # For progress bar


##############
#### Data ####
##############

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

###################
#### Functions ####
###################

def get_month_from_absolute_day(abs_day: int):
    """
    Get the month from an absolute day number.

    Parameters:
    - abs_day: int, absolute day number (1 to 365)

    Returns:
    - int, month number (1 = January, ..., 12 = December)
    """
    if abs_day < 1 or abs_day > 365:
        raise ValueError("Absolute day must be between 1 and 365.")

    # Find the month by iterating through the days in each month
    month = abs_date_table.filter(pl.col('Day_Num') == abs_day).select('Month').item()

    return month

##### Potential Evapotranspiration (PET) #####
def daylight_duration(latitude_deg, month):
    """
    Estimate daylight duration in hours based on latitude and month.

    Parameters:
    - latitude_deg: float, latitude in degrees (-90 to 90)
    - month: int, month number (1 = January, ..., 12 = December)

    Returns:
    - Daylight duration in hours (approximate)
    """
    # Convert latitude to radians
    lat_rad = math.radians(latitude_deg)

    # Approximate solar declination angle δ in radians (Cooper’s formula)
    day_of_year = [15, 45, 74, 105, 135, 162, 198, 228, 258, 288, 318, 344]  # mid-month days
    n = day_of_year[month - 1]
    decl_rad = math.radians(23.44) * math.sin(math.radians((360 / 365.0) * (n - 81)))

    # Calculate hour angle ω₀
    cos_omega = -math.tan(lat_rad) * math.tan(decl_rad)
    if cos_omega >= 1:
        # Polar night (no sunrise)
        return 0.0
    elif cos_omega <= -1:
        # Midnight sun (24h daylight)
        return 24.0
    else:
        omega = math.acos(cos_omega)
        return (2 * omega * 180 / math.pi) / 15  # convert radians to hours (15° per hour)
    

def calcualte_heat_index(air_temp: float):
    """
    This is calcualted as: I = (T/5)^1.514
    """
    
    heat_index = (max(0, air_temp)/5)**1.514

    return heat_index

def calcualte_a(heat_index: float):
    exponent = (6.75 * 10**(-7) * heat_index**3) - (7.71 * 10**(-5) * heat_index**2) + (0.01792 * heat_index) + 0.49239
    
    return exponent


def calculate_PET_location_based(monthly_temps, year: int, lat: float):
    """
    Calculate the Potential Evapotranspiration (PET) in mm per month using the Thornthwaite equation. It doesn't take into consideration the crop growing, only the location and the month.

    Parameters:
    - monthly_temps: air temperature series in degrees Celsius
    - year: year
    - lat: float, latitude in degrees

    Returns:
    - PET series in mm per month
    """
    # Checks there are 12 months of temp
    if len(monthly_temps) != 12:
        raise ValueError("monthly_temps must contain exactly 12 values for each month.")

    # Step 1: Calculate I (heat index)
    I = sum([(t / 5.0) ** 1.514 for t in monthly_temps if t > 0])

    # Step 2: Calculate exponent a
    a  = calcualte_a(I)
    
    # Step 3: Calculate PET Series 
    PET = []
    for i, T in enumerate(monthly_temps):
        if T <= 0:
            PET.append(0.0)
            continue

        # Days in the month
        N = monthrange(year, i+1)[1]  # 2024 is leap year-safe

        # Day length
        L = daylight_duration(lat, i + 1)

        # Thornthwaite formula
        PET_i = 16 * (L / 12) * (N / 30) * ((10 * T / I) ** a)
        PET.append(PET_i)

    return PET


##### Crop Coefficient (K) Functions #####
def correct_abs_date(num_day:int):
    if num_day > 365:
        num_day = num_day - 365
    
    return num_day

def create_KC_Curve(crop: str, climate: str):
    # First, check that both crops and climate are in the K_Crops table
    if crop not in K_Crops['Crop'].unique() or climate not in K_Crops['Climate_Zone'].unique():
        raise ValueError("Specified crop or climate not found in Kc data.")
    # else:
    #    print(f"Creating Kc curve for crop: {crop} in climate: {climate}")
    
    # Retrieve data for the specified crop and climate
    crop_data = K_Crops.filter((pl.col('Crop') == crop) & (pl.col('Climate_Zone') == climate))
    # print("Crop data retrieved:")
    # print(crop_data)

    # Stage dataframe for 365 days
    Kc_df = pl.DataFrame({
        "day_num": range(1, 366),
        "Stage_id": [0.0] * 365,
        "Kc": [0.0] * 365
    })

    # Calculating start and end dates for each phase
    KC_Dates = pl.DataFrame({
        "Stage": ["Planting", "Initial", "Development", "Mid", "Late"],
        "Stage_id": [1, 2, 3, 4, 5],
        "Start_Day": [0, 0, 0, 0, 0],
        "End_Day": [0, 0, 0, 0, 0]
    })

    # print("Retrieving stages dates for crop:", crop)
    
    # Get the phanse day number as an integer value
    planting_day_num_start = abs_date_table.filter(
        pl.col('Date') == crop_data['Planting_Greenup_Date']).select('Day_Num').item()
    planting_day_num_end = planting_day_num_start
    Initial_day_start = planting_day_num_end + 1
    Initial_day_end = Initial_day_start + crop_data['Initial_days'][0]
    Development_day_start = Initial_day_end + 1
    Development_day_end = Development_day_start + crop_data['Dev_Days'][0]
    Mid_day_start = Development_day_end + 1
    Mid_day_end = Mid_day_start + crop_data['Mid_Days'][0]
    Late_day_start = Mid_day_end + 1
    Late_day_end = Late_day_start + crop_data['Late_days'][0]

    # Correct absolute day numbers if they're above 365
    Initial_cday_start = correct_abs_date(Initial_day_start)
    Initial_cday_end = correct_abs_date(Initial_day_end)
    Development_cday_start = correct_abs_date(Development_day_start)    
    Development_cday_end = correct_abs_date(Development_day_end)
    Mid_cday_start = correct_abs_date(Mid_day_start)
    Mid_cday_end = correct_abs_date(Mid_day_end)
    Late_cday_start = correct_abs_date(Late_day_start)
    Late_cday_end = correct_abs_date(Late_day_end)

    # Filling planting date
    KC_Dates = KC_Dates.with_columns(
        pl.when(pl.col('Stage') == 'Planting').then(planting_day_num_start).when(pl.col('Stage') == 'Initial').then(Initial_cday_start).when(pl.col('Stage') == 'Development').then(Development_cday_start).when(pl.col('Stage') == 'Mid').then(Mid_cday_start).otherwise(Late_cday_start).alias('Start_Day'),
        pl.when(pl.col('Stage') == 'Planting').then(planting_day_num_end).when(pl.col('Stage') == 'Initial').then(Initial_cday_end).when(pl.col('Stage') == 'Development').then(Development_cday_end).when(pl.col('Stage') == 'Mid').then(Mid_cday_end).otherwise(Late_cday_end).alias('End_Day'),
    )

    # print("Stages dates:")
    # print(KC_Dates)
    
    # Filling Kc values
    # print("Assigning Kc values to stages")

    # Initial stage
    # print("Processing initial stage")
    Kc_df = Kc_df.with_columns(
        pl.when((pl.col("day_num") >= planting_day_num_start) & (pl.col("day_num") <= Initial_day_end))
        .then(crop_data['K_ini'][0])
        .otherwise(0)
        .alias("Kc"),

        pl.when(pl.col("day_num") == planting_day_num_start).then(1).when((pl.col("day_num") > planting_day_num_start) & (pl.col("day_num") <= Initial_day_end)).then(2).otherwise(pl.col("Stage_id")).alias("Stage_id")
    )    
    
    # development stage
    # print("Processing development stage")
    dev_stage_duration = Development_day_end - Development_day_start + 1
    Kc_dev_slope = (crop_data['K_mid'][0] - crop_data['K_ini'][0]) / dev_stage_duration

    day = Development_day_start

    while day in range(Development_day_start, Development_day_end + 1):
        day_corrected = correct_abs_date(day)
               
        Kc_value = crop_data['K_ini'][0] + Kc_dev_slope * (day - Development_day_start)
        
        Kc_df = Kc_df.with_columns(
            pl.when(pl.col("day_num") == day_corrected)
            .then(Kc_value)
            .otherwise(pl.col("Kc")).alias("Kc"),
            pl.when(pl.col("day_num") == day_corrected).then(3).otherwise(pl.col("Stage_id")).alias("Stage_id")
        )     

        day += 1

    # mid stage (Kc stays the same)
    # print("Processing mid stage")
    Kc_df = Kc_df.with_columns(
        pl.when((pl.col("day_num") >= Mid_cday_start) & (pl.col("day_num") <= Mid_cday_end))
        .then(crop_data['K_mid'][0])
        .otherwise(pl.col("Kc"))
        .alias("Kc"),

        pl.when((pl.col("day_num") >= Mid_cday_start) & (pl.col("day_num") <= Mid_cday_end)).then(4).otherwise(pl.col("Stage_id")).alias("Stage_id")
    )

    # late stage
    # print("Processing late stage")
    late_stage_duration = Late_day_end - Late_day_start + 1
    Kc_late_slope = (crop_data['K_Late'][0] - crop_data['K_mid'][0]) / late_stage_duration
    day = Late_day_start
    while day in range(Late_day_start, Late_day_end + 1):
        day_corrected = correct_abs_date(day)
               
        Kc_value = crop_data['K_mid'][0] + Kc_late_slope * (day - Late_day_start)
        
        Kc_df = Kc_df.with_columns(
            pl.when(pl.col("day_num") == day_corrected)
            .then(Kc_value)
            .otherwise(pl.col("Kc")).alias("Kc"),
            pl.when(pl.col("day_num") == day_corrected).then(5).otherwise(pl.col("Stage_id")).alias("Stage_id")
        )     

        day += 1

    # Kc_df =Kc_df.join(KC_Dates.select('Stage', 'Stage_id'), on="Stage_id", how="left")
    # print("Kc curve created successfully for crop:", crop, "in climate:", climate)
    
    return Kc_df

def monthly_KC_curve(crop: str, climate: str):
    daily_Kc_curve = create_KC_Curve(crop, climate)

    # Group by month and calculate average Kc for each month
    absday_month = abs_date_table.select(
        pl.col('Day_Num'),
        pl.col('Month')
    )
    
    monthly_Kc = daily_Kc_curve.join(absday_month, left_on='day_num', right_on='Day_Num', how='left').group_by('Month').agg(
        pl.col('Kc').mean().alias('Kc')
    ).sort('Month')

    return monthly_Kc


def calculate_PET_crop_based(crop: str, climate_zone:str, monthly_temps, year: int, lat: float):
    """
    Calculate crop-adjusted Potential Evapotranspiration (PET) on a daily, monthly, and annual basis.

    This function uses crop-specific crop coefficient (Kc) curves and location-based monthly PET
    estimates to derive daily PET values for a full calendar year, accounting for crop phenology
    and climate zone.

    Parameters
    ----------
    crop : str
        Name of the crop for which to calculate PET (must exist in the K_Crops table).
    climate_zone : str
        Climate zone associated with the crop (must match entries in the K_Crops table).
    monthly_temps : list or array-like
        Monthly average temperatures for the location (length 12).
    year : int
        Year used for determining leap years and calendar structure.
    lat : float
        Latitude of the location (in decimal degrees, used for PET calculation).

    Returns
    -------
    dict
        A dictionary with the following keys:
        - 'PET_Annaul' : float
            Total annual PET (mm/year) for the specified crop and location.
        - 'PET_Monthly' : pl.DataFrame
            Monthly PET values (mm/month), grouped by calendar month.
        - 'PET_Daily' : pl.DataFrame
            Daily PET values (mm/day), including day of year, crop stage, and stage Kc.

    Notes
    -----
    - This function assumes availability of global variables or dataframes:
      `K_Crops`, `abs_date_table`, `days_in_month_table`, and the functions
      `create_KC_Curve` and `calculate_PET_location_based`.
    - Monthly PET is computed from Kc-adjusted daily values.
    - Leap years are handled automatically through day correction logic.
    """

    # Create Kc curve for the specified crop and climate zone
    Kc_curve = create_KC_Curve(crop, climate_zone)

    # Calculate PET using the Kc values from the curve and monthly temperatures
    PET0 = calculate_PET_location_based(monthly_temps, year, lat)
    PET0 = pl.DataFrame({'Month': range(1, 13), 'PET0_Month': PET0})

    PET_daily = Kc_curve
    PET_daily = PET_daily.join(abs_date_table, left_on='day_num', right_on='Day_Num', how='left')
    PET_daily = PET_daily.join(PET0, on='Month', how='left').join(days_in_month_table, on='Month', how='left')
    PET_daily = PET_daily.with_columns(
        PET_Daily = pl.col('Kc') * pl.col('PET0_Month')/pl.col("Days_in_Month")# Convert monthly PET to daily PET
    )

    PET_Monthly = PET_daily.group_by('Month').agg(
        pl.col('PET_Daily').sum().alias('PET_Monthly')
    ).sort('Month')

    PET_Annual = PET_Monthly['PET_Monthly'].sum()

    results = {
        'PET_Annaul': PET_Annual,
        'PET_Monthly': PET_Monthly,
        'PET_Daily': PET_daily
    }

    print("PET calculation completed successfully for crop:", crop, "in climate zone:", climate_zone)
    return results


def _compute_crop_based_pet_rasters(
    crop_name: str,
    landuse: Union[str, Path, np.ndarray],
    pet_base_raster_path: str,
    thermal_zone_raster_path: str,
    output_monthly_path: Optional[Union[str, Path]] = None,
    output_annual_path: Optional[Union[str, Path]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Core implementation for crop-based PET raster generation."""

    if crop_name not in K_Crops['Crop'].unique():
        raise ValueError(f"Crop '{crop_name}' not found in K_Crops table.")

    pet_base_path = Path(pet_base_raster_path)
    thermal_path = Path(thermal_zone_raster_path)

    with rasterio.open(pet_base_path) as pet_src:
        pet_base = pet_src.read().astype("float32")
        profile = pet_src.profile.copy()
        base_crs = pet_src.crs
        base_transform = pet_src.transform
        base_width = pet_src.width
        base_height = pet_src.height

    with rasterio.open(thermal_path) as thermal_src:
        if (
            thermal_src.crs != base_crs
            or thermal_src.transform != base_transform
            or thermal_src.width != base_width
            or thermal_src.height != base_height
        ):
            raise ValueError("CRS/transform/shape mismatch")

        thermal_zones = thermal_src.read(1).astype(int)

    if isinstance(landuse, (str, Path)):
        landuse_path = Path(landuse)
        with rasterio.open(landuse_path) as landuse_src:
            if (
                landuse_src.crs != base_crs
                or landuse_src.transform != base_transform
                or landuse_src.width != base_width
                or landuse_src.height != base_height
            ):
                raise ValueError("CRS/transform/shape mismatch")

            landuse_data = landuse_src.read(1).astype(int)
    else:
        landuse_data = np.asarray(landuse).astype(int)
        if landuse_data.shape != (base_height, base_width):
            raise ValueError("Land-use array must match base raster dimensions.")

    H, W = pet_base.shape[1:]
    pet_monthly = np.full_like(pet_base, np.nan, dtype="float32")
    pet_annual = np.full((H, W), np.nan, dtype="float32")

    unique_groups = list(zone_ids_by_group)
    kc_by_group = {}
    for grp in tqdm(unique_groups, desc="Precomputing Kc curves for group"):
        print(f"Precomputing Kc curve for group: {grp}")
        kc_df = monthly_KC_curve(crop_name, grp)
        kc_by_group[grp] = kc_df.sort("Month")["Kc"].to_numpy()

    for grp, kc_vec in tqdm(kc_by_group.items(), desc="Applying Kc to thremal groups"):
        valid_zones = zone_ids_by_group[grp]
        mask = np.isin(thermal_zones, valid_zones) & (landuse_data == 1)
        if not mask.any():
            continue

        monthly_crop = pet_base[:, mask] * kc_vec[:, None]
        pet_monthly[:, mask] = monthly_crop.astype("float32")
        pet_annual[mask] = np.nansum(monthly_crop, axis=0).astype("float32")

    if output_monthly_path is not None:
        monthly_profile = profile.copy()
        monthly_profile.update(count=pet_monthly.shape[0], dtype="float32", nodata=np.nan)
        with rasterio.open(output_monthly_path, "w", **monthly_profile) as dst:
            dst.write(pet_monthly)

    if output_annual_path is not None:
        annual_profile = profile.copy()
        annual_profile.update(count=1, dtype="float32", nodata=np.nan)
        with rasterio.open(output_annual_path, "w", **annual_profile) as dst:
            dst.write(pet_annual, 1)

    return pet_monthly, pet_annual


def calculate_crop_based_PET_raster_optimized(
    crop_name: str,
    landuse_raster_path: str,
    output_monthly_path: str,
    output_annual_path: str,
    pet_base_raster_path: str = "SOC_Data_Processing/uhth_pet_locationonly.tif",
    thermal_zone_raster_path: str = "SOC_Data_Processing/uhth_thermal_climates.tif"
):
    """Generate crop-specific PET rasters using optimized raster operations.

    Parameters
    ----------
    crop_name : str
        Name of the crop to calculate PET for. Must exist in the :data:`K_Crops` table.
    landuse_raster_path : str
        Path to the land-use raster aligned with the base PET and thermal zone rasters.
    output_monthly_path : str
        Path where the 12-band monthly PET raster will be written.
    output_annual_path : str
        Path where the single-band annual PET raster will be written.
    pet_base_raster_path : str, optional
        Path to the base PET raster (default ``"SOC_Data_Processing/uhth_pet_locationonly.tif"``).
    thermal_zone_raster_path : str, optional
        Path to the thermal zone raster (default ``"SOC_Data_Processing/uhth_thermal_climates.tif"``).

    Raises
    ------
    ValueError
        If the crop name is not present in :data:`K_Crops` or the rasters are misaligned.
    """

    _compute_crop_based_pet_rasters(
        crop_name=crop_name,
        landuse=landuse_raster_path,
        pet_base_raster_path=pet_base_raster_path,
        thermal_zone_raster_path=thermal_zone_raster_path,
        output_monthly_path=output_monthly_path,
        output_annual_path=output_annual_path,
    )

    print(f"PET calculation completed for crop '{crop_name}' succesfully. Storing...")




def calculate_crop_based_PET_raster_vPipeline(
    crop_name: str,
    landuse_array: np.ndarray,
    output_monthly_path: str,
    pet_base_raster_path: str = "data/soil_weather/uhth_pet_locationonly.tif",
    thermal_zone_raster_path: str = "data/soil_weather/uhth_thermal_climates.tif"
):
    """Generate crop-specific monthly PET rasters for array-based workflows.

    Parameters
    ----------
    crop_name : str
        Name of the crop to calculate PET for. Must exist in the :data:`K_Crops` table.
    landuse_array : numpy.ndarray
        Land-use mask already loaded in memory. The array must align with the base PET raster.
    output_monthly_path : str
        Path where the 12-band monthly PET raster will be written.
    pet_base_raster_path : str, optional
        Path to the base PET raster (default ``"data/soil_weather/uhth_pet_locationonly.tif"``).
    thermal_zone_raster_path : str, optional
        Path to the thermal zone raster (default ``"data/soil_weather/uhth_thermal_climates.tif"``).

    Returns
    -------
    numpy.ndarray
        The computed monthly PET raster stack (12, H, W) masked by land use and thermal zone.

    Raises
    ------
    ValueError
        If the crop name is not present in :data:`K_Crops` or the rasters are misaligned.
    """

    pet_monthly, _ = _compute_crop_based_pet_rasters(
        crop_name=crop_name,
        landuse=landuse_array,
        pet_base_raster_path=pet_base_raster_path,
        thermal_zone_raster_path=thermal_zone_raster_path,
        output_monthly_path=output_monthly_path,
    )

    print(f"PET calculation completed for crop '{crop_name}' succesfully")

    return pet_monthly
