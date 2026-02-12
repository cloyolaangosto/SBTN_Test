import numpy as np
import geopandas as gpd
from shapely.geometry import box
from rasterio.transform import from_origin
import rasterio

from sbtn_leaf.cropcalcs import (
    create_crop_yield_raster,
    create_crop_yield_raster_withIrrigationPracticeScaling,
)


def _write_raster(path, array, transform, *, nodata=None, crs="EPSG:4326"):
    array = np.asarray(array, dtype="float32")
    height, width = array.shape
    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": "float32",
        "transform": transform,
        "crs": crs,
        "nodata": nodata,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array, 1)


def _make_fao_gdf(
    avg_column,
    ratio_column,
    *,
    avg_value,
    ratio_value,
    sd_column="sd_yield",
    sd_value=0.0,
    crs="EPSG:4326",
):
    geom = [box(0, 0, 2, 2)]
    data = {
        avg_column: [avg_value],
        ratio_column: [ratio_value],
        sd_column: [sd_value],
    }
    return gpd.GeoDataFrame(data, geometry=geom, crs=crs)


def test_create_crop_yield_raster_base(tmp_path):
    transform = from_origin(0, 2, 1, 1)

    lu_path = tmp_path / "lu.tif"
    _write_raster(lu_path, np.ones((2, 2), dtype="float32"), transform)

    spam_path = tmp_path / "spam.tif"
    spam_array = np.array([[10, -9999], [30, 40]], dtype="float32")
    _write_raster(spam_path, spam_array, transform, nodata=-9999.0)

    fao_gdf = _make_fao_gdf("avg_yield_1423", "ratio_yield_20_toavg", avg_value=1000.0, ratio_value=0.5)

    out_path = tmp_path / "yield.tif"
    create_crop_yield_raster(str(lu_path), fao_gdf, str(spam_path), str(out_path))

    with rasterio.open(out_path) as src:
        data = src.read(1)

    expected = np.array([[5.0, 1.0], [15.0, 20.0]], dtype="float32")
    np.testing.assert_allclose(data, expected, rtol=1e-6, atol=1e-6)


def test_create_crop_yield_raster_with_irrigation_scaling(tmp_path):
    transform = from_origin(0, 2, 1, 1)

    lu_path = tmp_path / "lu.tif"
    _write_raster(lu_path, np.ones((2, 2), dtype="float32"), transform)

    spam_path = tmp_path / "spam.tif"
    spam_array = np.array([[10, -9999], [30, 40]], dtype="float32")
    _write_raster(spam_path, spam_array, transform, nodata=-9999.0)

    all_path = tmp_path / "all.tif"
    # Keep one nodata pixel in all-yields coverage while leaving a valid fallback
    # value for the SPAM nodata pixel.
    all_array = np.array([[10, 20], [30, -9999]], dtype="float32")
    _write_raster(all_path, all_array, transform, nodata=-9999.0)

    irr_path = tmp_path / "irr.tif"
    irr_array = np.array([[5, 18], [24, -9999]], dtype="float32")
    _write_raster(irr_path, irr_array, transform, nodata=-9999.0)

    rf_path = tmp_path / "rf.tif"
    rf_array = np.array([[3, 10], [12, -9999]], dtype="float32")
    _write_raster(rf_path, rf_array, transform, nodata=-9999.0)

    fao_gdf = _make_fao_gdf("avg_yield", "yld_ratio", avg_value=1000.0, ratio_value=0.5)

    out_path = tmp_path / "yield_scaled.tif"
    create_crop_yield_raster_withIrrigationPracticeScaling(
        str(lu_path),
        fao_gdf,
        str(spam_path),
        str(out_path),
        irr_yield_scaling="irr",
        all_fp=str(all_path),
        irr_fp=str(irr_path),
        rf_fp=str(rf_path),
        apply_ecoregion_fill=False,
    )

    with rasterio.open(out_path) as src:
        data = src.read(1)

    avg_wat_ratio = np.mean([5 / 10, 18 / 20, 24 / 30])
    # Direct SPAM-backed pixels should keep SPAM values in GAEZ mode.
    np.testing.assert_allclose(data[[0, 1, 1], [0, 0, 1]], [10.0, 30.0, 40.0], rtol=1e-6, atol=1e-6)

    fallback_value = 20.0 * avg_wat_ratio

    # The fallback for the SPAM nodata pixel should use all-SPAM yields scaled
    # by the average irrigation ratio, not the FAO yield ratio.
    assert fallback_value != 20.0 * 0.5
    np.testing.assert_allclose(data[0, 1], fallback_value, rtol=1e-6, atol=1e-6)
