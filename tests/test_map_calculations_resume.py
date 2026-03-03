"""Tests for checkpointed/resumable raster-to-tabular builds."""

from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from sbtn_leaf import map_calculations as mc


def _build_master_gdf() -> gpd.GeoDataFrame:
    rows = [
        {
            "ECO_ID": 1,
            "ECO_NAME": "Region A",
            "BIOME_NUM": 10,
            "BIOME_NAME": "Biome",
            "REALM": "Realm",
            "NNH": "",
            "SHAPE_LENG": 0.0,
            "SHAPE_AREA": 0.0,
            "NNH_NAME": "",
            "COLOR": "",
            "COLOR_BIO": "",
            "COLOR_NNH": "",
            "LICENSE": "",
            "geometry": Point(0, 0),
        },
        {
            "ECO_ID": 2,
            "ECO_NAME": "Region B",
            "BIOME_NUM": 10,
            "BIOME_NAME": "Biome",
            "REALM": "Realm",
            "NNH": "",
            "SHAPE_LENG": 0.0,
            "SHAPE_AREA": 0.0,
            "NNH_NAME": "",
            "COLOR": "",
            "COLOR_BIO": "",
            "COLOR_NNH": "",
            "LICENSE": "",
            "geometry": Point(1, 1),
        },
    ]
    return gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:6933")


def _flow_output(flow_name: str):
    df = pd.DataFrame(
        {
            "ECO_ID": [1, 2],
            "cf": [1.0, 2.0],
            "cf_median": [1.1, 2.1],
            "cf_std": [0.1, 0.2],
        }
    )
    gdf = gpd.GeoDataFrame({"ECO_ID": [1, 2], "geometry": [Point(0, 0), Point(1, 1)]}, geometry="geometry", crs="EPSG:6933")
    return df, gdf


def test_build_cfs_gpkg_resume_checkpoint(tmp_path, monkeypatch):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    for name in ("flow_a.tif", "flow_b.tif", "flow_c.tif"):
        (input_dir / name).write_text("placeholder")

    state = {"flow_b_failures": 0, "calls": []}

    def flaky_calc(*, raster_input_filepath, flow_name, **kwargs):
        state["calls"].append(Path(raster_input_filepath).name)
        if flow_name == "flow_b" and state["flow_b_failures"] == 0:
            state["flow_b_failures"] += 1
            raise RuntimeError("simulated failure")
        return _flow_output(flow_name)

    monkeypatch.setattr(mc, "calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers", flaky_calc)

    with pytest.raises(RuntimeError):
        mc.build_cfs_gpkg_from_rasters(
            input_folder=str(input_dir),
            output_folder=str(output_dir) + "/",
            layer_name="test_layer",
            master_gdf=_build_master_gdf(),
            master_key="ECO_ID",
            result_key="ECO_ID",
            cf_name="impact",
            cf_unit="unit",
            area_type="ecoregion",
            write_gpkg=False,
            fail_fast=True,
            resume=True,
        )

    csv_path = output_dir / "impact_ecoregion.csv"
    checkpoint_path = output_dir / "impact_ecoregion.checkpoint.json"
    assert csv_path.exists()
    assert checkpoint_path.exists()

    first_run = pd.read_csv(csv_path)
    assert set(first_run["flow_name"].unique()) == {"flow_a"}

    state["calls"] = []

    mc.build_cfs_gpkg_from_rasters(
        input_folder=str(input_dir),
        output_folder=str(output_dir) + "/",
        layer_name="test_layer",
        master_gdf=_build_master_gdf(),
        master_key="ECO_ID",
        result_key="ECO_ID",
        cf_name="impact",
        cf_unit="unit",
        area_type="ecoregion",
        write_gpkg=False,
        fail_fast=True,
        resume=True,
    )

    assert "flow_a.tif" not in state["calls"]
    assert state["calls"] == ["flow_b.tif", "flow_c.tif"]

    rerun = pd.read_csv(csv_path)
    assert set(rerun["flow_name"].unique()) == {"flow_a", "flow_b", "flow_c"}


def test_build_cfs_gpkg_reset_false_preserves_existing_layers(tmp_path, monkeypatch):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    for name in ("flow_a.tif", "flow_b.tif"):
        (input_dir / name).write_text("placeholder")

    gpkg_path = output_dir / "impact_ecoregion.gpkg"
    gpkg_path.write_text("placeholder")

    def stable_calc(*, flow_name, **kwargs):
        return _flow_output(flow_name)

    calls = []

    def fake_write_df(df, path, layer, driver, append, **kwargs):
        calls.append({"layer": layer, "append": append, "cols": list(df.columns)})

    def fake_list_layers(path):
        return [("geometry_layer", "Point"), ("test_layer", "None"), ("test_layer_metadata", "None")]

    def fake_read_df(path, layer=None, max_features=None, read_geometry=None, **kwargs):
        if layer == "test_layer":
            return pd.DataFrame(columns=["ECO_ID", "flow_name", "cf", "cf_median", "cf_std", "_source_file"])
        if layer == "test_layer_metadata":
            return pd.DataFrame(columns=["flow_name", "impact_category", "unit", "source_file"])
        return pd.DataFrame()

    monkeypatch.setattr(mc, "calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers", stable_calc)
    monkeypatch.setattr(mc, "write_df", fake_write_df)
    monkeypatch.setattr(mc, "ogr_list_layers", fake_list_layers)
    monkeypatch.setattr(mc, "read_df", fake_read_df)

    mc.build_cfs_gpkg_from_rasters(
        input_folder=str(input_dir),
        output_folder=str(output_dir) + "/",
        layer_name="test_layer",
        master_gdf=_build_master_gdf(),
        master_key="ECO_ID",
        result_key="ECO_ID",
        cf_name="impact",
        cf_unit="unit",
        area_type="ecoregion",
        write_gpkg=True,
        reset_gpkg=False,
        resume=False,
    )

    geometry_writes = [c for c in calls if c["layer"] == "geometry_layer"]
    assert not geometry_writes, "geometry_layer should not be overwritten when reset_gpkg=False"

    value_writes = [c for c in calls if c["layer"] == "test_layer"]
    assert value_writes
    assert value_writes[0]["append"] is True
