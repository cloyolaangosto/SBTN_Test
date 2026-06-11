"""Tests for the data-coverage threshold (NaN assignment) in CF aggregation."""

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

from sbtn_leaf.map_calculations import (
    calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers,
)

# Equal-area CRS so that ``geom.area`` and pixel areas are both in square metres
# and no reprojection happens inside the calculator.
EQUAL_AREA_CRS = "EPSG:6933"
NODATA = -9999.0


def _write_partial_coverage_raster(path, covered_rows, n=10, pixel=1000.0):
    """Write an ``n`` x ``n`` raster (pixel x pixel metres) in an equal-area CRS.

    The top ``covered_rows`` rows hold a valid constant value; the rest are
    nodata. With ``n`` rows the covered fraction of a polygon spanning the whole
    raster is ``covered_rows / n``.
    """
    data = np.full((n, n), NODATA, dtype=np.float32)
    data[:covered_rows, :] = 5.0
    transform = from_origin(0, n * pixel, pixel, pixel)
    profile = {
        "driver": "GTiff",
        "height": n,
        "width": n,
        "count": 1,
        "dtype": "float32",
        "transform": transform,
        "crs": EQUAL_AREA_CRS,
        "nodata": NODATA,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data, 1)
    # Polygon covering the full raster extent -> total area = (n*pixel)**2.
    extent = box(0, 0, n * pixel, n * pixel)
    er_gdf = gpd.GeoDataFrame(
        {
            "ECO_ID": [1],
            "ECO_NAME": ["synthetic"],
            "BIOME_NAME": ["test_biome"],
            "geometry": [extent],
        },
        crs=EQUAL_AREA_CRS,
    )
    return er_gdf


def _run(path, er_gdf, threshold):
    results_df, _ = calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers(
        raster_input_filepath=str(path),
        cf_name="test_cf",
        cf_unit="unit",
        flow_name="flow",
        area_type="ecoregion",
        er_gdf=er_gdf,
        equal_area_crs=EQUAL_AREA_CRS,
        min_coverage_fraction=threshold,
        return_gdf=False,
        suppress_logging=True,
    )
    return results_df


def test_below_threshold_assigns_nan(tmp_path):
    # 30% of the polygon is covered by valid data.
    er_gdf = _write_partial_coverage_raster(tmp_path / "cov30.tif", covered_rows=3)
    results_df = _run(tmp_path / "cov30.tif", er_gdf, threshold=0.5)

    # The region is still present in the output but its stats are NaN.
    assert len(results_df) == 1
    row = results_df.iloc[0]
    assert np.isnan(row["cf"])
    assert np.isnan(row["cf_median"])
    assert np.isnan(row["cf_std"])


def test_above_threshold_keeps_value(tmp_path):
    # 30% coverage with a 20% threshold -> value retained.
    er_gdf = _write_partial_coverage_raster(tmp_path / "cov30b.tif", covered_rows=3)
    results_df = _run(tmp_path / "cov30b.tif", er_gdf, threshold=0.2)

    row = results_df.iloc[0]
    assert np.isfinite(row["cf"])
    assert row["cf"] == 5.0


def test_threshold_disabled_keeps_value(tmp_path):
    # A sliver of coverage (10%) is retained when the gate is disabled (0.0).
    er_gdf = _write_partial_coverage_raster(tmp_path / "cov10.tif", covered_rows=1)
    results_df = _run(tmp_path / "cov10.tif", er_gdf, threshold=0.0)

    row = results_df.iloc[0]
    assert np.isfinite(row["cf"])
    assert row["cf"] == 5.0


def test_threshold_boundary_inclusive(tmp_path):
    # Coverage exactly at the threshold (30% == 0.30) is retained (>= comparison).
    er_gdf = _write_partial_coverage_raster(tmp_path / "cov30c.tif", covered_rows=3)
    results_df = _run(tmp_path / "cov30c.tif", er_gdf, threshold=0.3)

    row = results_df.iloc[0]
    assert np.isfinite(row["cf"])
    assert row["cf"] == 5.0


def test_invalid_threshold_raises(tmp_path):
    er_gdf = _write_partial_coverage_raster(tmp_path / "covx.tif", covered_rows=3)
    import pytest

    with pytest.raises(ValueError):
        _run(tmp_path / "covx.tif", er_gdf, threshold=1.5)
