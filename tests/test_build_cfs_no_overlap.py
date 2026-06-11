"""Regression tests for no-overlap handling in ``build_cfs_gpkg_from_rasters``.

These guard against the duplicate-key join bug where a CF value computed for an
overlapping polygon was broadcast onto same-keyed polygons that have no overlap
with raster data pixels (observed for ecoregions, whose raw shapefile splits a
single ``ECO_ID`` across many polygon rows).
"""

import logging
import math

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

from sbtn_leaf.map_calculations import build_cfs_gpkg_from_rasters

NODATA = -9999.0
DATA_VALUE = 50.0

# Columns the ecoregion code path requires: csv enrichment + drop_cols + labels.
_ECO_DROP_COLS = [
    "NNH", "SHAPE_LENG", "SHAPE_AREA", "NNH_NAME",
    "COLOR", "COLOR_BIO", "COLOR_NNH", "LICENSE",
]


def _write_corner_raster(path):
    """10x10 EPSG:4326 raster, valid data only in the top-left 2x2 pixels.

    Top-left pixel origin is (lon=0, lat=10) with 1-degree pixels, so the data
    block covers lon[0, 2], lat[8, 10]; everything else is nodata.
    """
    arr = np.full((10, 10), NODATA, dtype="float32")
    arr[0:2, 0:2] = DATA_VALUE  # top-left corner -> lon[0,2], lat[8,10]
    transform = from_origin(0, 10, 1, 1)
    profile = {
        "driver": "GTiff",
        "height": 10,
        "width": 10,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": transform,
        "nodata": NODATA,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr, 1)


def _make_eco_gdf(rows):
    """Build a minimal ecoregion GeoDataFrame with all columns the code touches.

    ``rows`` is a list of (eco_id, geometry) tuples.
    """
    records = []
    for eco_id, geom in rows:
        rec = {
            "ECO_ID": eco_id,
            "ECO_NAME": f"eco_{eco_id}",
            "BIOME_NUM": 1,
            "BIOME_NAME": "biome_1",
            "REALM": "realm_1",
            "geometry": geom,
        }
        for col in _ECO_DROP_COLS:
            rec[col] = 0
        records.append(rec)
    return gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")


def _run(tmp_path, gdf, logger):
    raster_dir = tmp_path / "rasters"
    out_dir = tmp_path / "out"
    raster_dir.mkdir()
    out_dir.mkdir()
    _write_corner_raster(raster_dir / "soc_tomato.tif")

    _, results_df = build_cfs_gpkg_from_rasters(
        input_folder=str(raster_dir) + "/",
        output_folder=str(out_dir) + "/",
        layer_name="soc_eco",
        master_gdf=gdf,
        master_key="ECO_ID",
        result_key="er_id",
        cf_name="soc",
        cf_unit="t C/ha",
        area_type="ecoregion",
        write_gpkg=False,
        reset_gpkg=True,
        logger=logger,
    )
    return results_df


def _mean_value(results_df, eco_id):
    sel = results_df[(results_df["ECO_ID"] == eco_id) & (results_df["metric"] == "cf_mean")]
    return sel


def test_duplicate_eco_id_does_not_bleed_to_non_overlapping_polygon(tmp_path, caplog):
    """ECO_ID shared by an overlapping + a non-overlapping polygon: after the
    internal dissolve there is exactly one row per ECO_ID, the value comes from
    the overlapping part, and a fully non-overlapping ECO_ID stays NaN."""
    logger = logging.getLogger("test_cfs_dup")
    gdf = _make_eco_gdf([
        (1, box(0, 8, 2, 10)),   # overlaps the data corner
        (1, box(0, 0, 2, 2)),    # same ECO_ID, over nodata only (no overlap)
        (2, box(4, 4, 6, 6)),    # unique ECO_ID, fully over nodata (no overlap)
    ])

    with caplog.at_level(logging.WARNING, logger="test_cfs_dup"):
        results_df = _run(tmp_path, gdf, logger)

    # Duplicate key was detected and dissolved -> warning emitted.
    assert any("duplicated" in r.message for r in caplog.records)

    # Exactly one row per ECO_ID per metric (no Cartesian duplication).
    counts = results_df[results_df["metric"] == "cf_mean"]["ECO_ID"].value_counts()
    assert counts.get(1, 0) == 1
    assert counts.get(2, 0) == 1

    # Overlapping ecoregion gets the data value; non-overlapping one is NaN.
    eco1 = _mean_value(results_df, 1)["value"].iloc[0]
    eco2 = _mean_value(results_df, 2)["value"].iloc[0]
    assert math.isclose(eco1, DATA_VALUE, rel_tol=1e-3)
    assert pd.isna(eco2)


def test_unique_keys_preserve_nan_for_no_overlap(tmp_path, caplog):
    """Control: with unique ECO_IDs no dissolve runs, and non-overlapping
    ecoregions correctly come back as NaN."""
    logger = logging.getLogger("test_cfs_unique")
    gdf = _make_eco_gdf([
        (1, box(0, 8, 2, 10)),   # overlaps the data corner
        (2, box(4, 4, 6, 6)),    # no overlap
        (3, box(7, 1, 9, 3)),    # no overlap
    ])

    with caplog.at_level(logging.WARNING, logger="test_cfs_unique"):
        results_df = _run(tmp_path, gdf, logger)

    # No duplicate keys -> no dissolve warning.
    assert not any("duplicated" in r.message for r in caplog.records)

    assert math.isclose(_mean_value(results_df, 1)["value"].iloc[0], DATA_VALUE, rel_tol=1e-3)
    assert pd.isna(_mean_value(results_df, 2)["value"].iloc[0])
    assert pd.isna(_mean_value(results_df, 3)["value"].iloc[0])
