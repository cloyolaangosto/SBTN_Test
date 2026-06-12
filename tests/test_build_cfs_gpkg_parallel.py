"""Tests for build_cfs_gpkg_from_rasters parallel execution and tqdm/logging safety."""

import logging
import threading
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

import sbtn_leaf.map_calculations as mc

CRS = "EPSG:6933"
PIXEL = 100.0  # meters


def _write_raster(path: Path, left_value: float, right_value: float) -> None:
    """10x10 raster spanning x 0..1000, y 0..1000; left half=left_value, right half=right_value."""
    data = np.full((10, 10), right_value, dtype=np.float32)
    data[:, :5] = left_value
    profile = {
        "driver": "GTiff",
        "height": 10,
        "width": 10,
        "count": 1,
        "dtype": "float32",
        "transform": from_origin(0, 1000, PIXEL, PIXEL),
        "crs": CRS,
        "nodata": -9999.0,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data, 1)


def _make_inputs(tmp_path: Path, n_rasters: int = 4):
    input_folder = tmp_path / "rasters"
    input_folder.mkdir()
    for i in range(n_rasters):
        _write_raster(input_folder / f"flow_{i}.tif", left_value=10.0 * (i + 1), right_value=20.0 * (i + 1))

    master_gdf = gpd.GeoDataFrame(
        {"ADM0_NAME": ["Left", "Right"]},
        geometry=[box(0, 0, 500, 1000), box(500, 0, 1000, 1000)],
        crs=CRS,
    )
    return input_folder, master_gdf


def _run(input_folder: Path, master_gdf: gpd.GeoDataFrame, output_folder: Path, *, max_workers: int, logger) -> pd.DataFrame:
    output_folder.mkdir(exist_ok=True)
    _, results_df = mc.build_cfs_gpkg_from_rasters(
        input_folder=str(input_folder),
        output_folder=str(output_folder) + "/",
        layer_name="test_layer",
        master_gdf=master_gdf,
        master_key="ADM0_NAME",
        result_key="country",
        cf_name="test_cf",
        cf_unit="unit",
        area_type="country",
        reset_gpkg=True,
        write_gpkg=False,
        logger=logger,
        max_workers=max_workers,
    )
    return results_df.sort_values(["ADM0_NAME", "flow_name", "metric"]).reset_index(drop=True)


def test_parallel_matches_sequential_and_uses_workers(tmp_path, monkeypatch, capsys):
    input_folder, master_gdf = _make_inputs(tmp_path)
    logger = logging.getLogger("test_build_cfs_gpkg")
    logger.setLevel(logging.INFO)

    sequential = _run(input_folder, master_gdf, tmp_path / "seq", max_workers=1, logger=logger)

    # Wrap the calculator to record which threads run it and force two calls
    # to overlap, proving the executor actually runs rasters concurrently.
    original = mc.calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers
    thread_names: set[str] = set()
    barrier = threading.Barrier(2)

    def tracking_calculator(*args, **kwargs):
        thread_names.add(threading.current_thread().name)
        try:
            barrier.wait(timeout=5)
        except threading.BrokenBarrierError:
            pass
        return original(*args, **kwargs)

    monkeypatch.setattr(
        mc,
        "calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers",
        tracking_calculator,
    )

    capsys.readouterr()  # drop output from the sequential run
    parallel = _run(input_folder, master_gdf, tmp_path / "par", max_workers=4, logger=logger)

    pd.testing.assert_frame_equal(sequential, parallel)
    assert len(thread_names) >= 2, f"expected >=2 worker threads, got {thread_names}"

    # Worker-thread log records are routed through tqdm.write by
    # logging_redirect_tqdm; any handler failure (e.g. refreshing a
    # half-constructed bar) would be printed to stderr by logging.
    err = capsys.readouterr().err
    assert "Logging error" not in err
    assert "AttributeError" not in err


def test_results_values_are_area_weighted(tmp_path):
    input_folder, master_gdf = _make_inputs(tmp_path, n_rasters=2)
    results = _run(input_folder, master_gdf, tmp_path / "out", max_workers=2, logger=None)

    # Each region is covered by a single constant value, so mean == median and std == 0.
    for i in range(2):
        flow = results[results["flow_name"] == f"flow_{i}"]
        left_mean = flow[(flow["ADM0_NAME"] == "Left") & (flow["metric"] == "cf_mean")]["value"].item()
        right_mean = flow[(flow["ADM0_NAME"] == "Right") & (flow["metric"] == "cf_mean")]["value"].item()
        assert left_mean == 10.0 * (i + 1)
        assert right_mean == 20.0 * (i + 1)
        assert (flow[flow["metric"] == "cf_std"]["value"] == 0).all()
