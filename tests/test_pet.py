import numpy as np
import rasterio
from rasterio.transform import from_origin

from sbtn_leaf.PET import (
    calculate_PET_location_based,
    calculate_crop_based_PET_raster_optimized,
    calculate_crop_based_PET_raster_vPipeline,
)


def test_january_daylight_duration_does_not_raise_index_error():
    """Ensure the daylight lookup uses a valid month index for January."""

    temps = [10.0] * 12

    # A successful calculation implies no IndexError was raised for January.
    pet_values = calculate_PET_location_based(temps, 2024, 0.0)

    assert len(pet_values) == 12


def test_raster_helpers_produce_identical_results(tmp_path):
    """The raster helper should yield identical results for path and array land-use inputs."""

    width = height = 2
    transform = from_origin(0, 2, 1, 1)

    base_month = np.array([[1, 2], [3, 4]], dtype="float32")
    pet_data = np.stack([base_month + month for month in range(12)], axis=0)
    thermal_data = np.ones((height, width), dtype="int16")
    landuse_data = np.ones((height, width), dtype="int16")

    pet_profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 12,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": transform,
    }

    thermal_profile = pet_profile.copy()
    thermal_profile.update(count=1, dtype="int16")

    pet_path = tmp_path / "pet.tif"
    thermal_path = tmp_path / "thermal.tif"
    landuse_path = tmp_path / "landuse.tif"
    monthly_out_path = tmp_path / "monthly_from_path.tif"
    annual_out_path = tmp_path / "annual_from_path.tif"
    monthly_out_array_path = tmp_path / "monthly_from_array.tif"

    with rasterio.open(pet_path, "w", **pet_profile) as dst:
        dst.write(pet_data)

    with rasterio.open(thermal_path, "w", **thermal_profile) as dst:
        dst.write(thermal_data, 1)

    with rasterio.open(landuse_path, "w", **thermal_profile) as dst:
        dst.write(landuse_data, 1)

    crop_name = "Apple"

    calculate_crop_based_PET_raster_optimized(
        crop_name=crop_name,
        landuse_raster_path=str(landuse_path),
        output_monthly_path=str(monthly_out_path),
        output_annual_path=str(annual_out_path),
        pet_base_raster_path=str(pet_path),
        thermal_zone_raster_path=str(thermal_path),
    )

    array_monthly = calculate_crop_based_PET_raster_vPipeline(
        crop_name=crop_name,
        landuse_array=landuse_data,
        output_monthly_path=str(monthly_out_array_path),
        pet_base_raster_path=str(pet_path),
        thermal_zone_raster_path=str(thermal_path),
    )

    with rasterio.open(monthly_out_path) as src:
        monthly_from_path = src.read()

    with rasterio.open(monthly_out_array_path) as src:
        monthly_from_array = src.read()

    with rasterio.open(annual_out_path) as src:
        annual_from_path = src.read(1)

    np.testing.assert_allclose(monthly_from_path, monthly_from_array)
    np.testing.assert_allclose(annual_from_path, np.nansum(array_monthly, axis=0))
