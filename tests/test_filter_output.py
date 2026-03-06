"""Tests for _clip_soc_output (RothC_Raster) and filter_output_raster (map_calculations)."""

import numpy as np
import rasterio
from rasterio.transform import from_origin

from sbtn_leaf.RothC_Raster import _clip_soc_output
from sbtn_leaf.map_calculations import filter_output_raster


def _write_multiband_raster(path, data, transform, *, crs="EPSG:4326"):
    """Write a (bands, height, width) float32 array as a GeoTIFF."""
    data = np.asarray(data, dtype="float32")
    n_bands, height, width = data.shape
    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": n_bands,
        "dtype": "float32",
        "transform": transform,
        "crs": crs,
        "nodata": None,
    }
    with rasterio.open(path, "w", **profile) as dst:
        for i in range(n_bands):
            dst.write(data[i], i + 1)


# ---------------------------------------------------------------------------
# Tests for _clip_soc_output
# ---------------------------------------------------------------------------

def test_clip_soc_output_caps_at_max_gain():
    soc0 = np.array([[10.0, 20.0]], dtype="float32")
    # 3 years of SOC data (year 0, 1, 2)
    soc_annual = np.array([
        [[10.0, 20.0]],  # year 0
        [[18.0, 30.0]],  # year 1: gain of 8 and 10
        [[30.0, 50.0]],  # year 2: gain of 20 and 30
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0, max_annual_gain=5.0,
                              global_percentile_cap=None)

    # Year 0 unchanged
    np.testing.assert_allclose(result[0], soc_annual[0])
    # Year 1: cap = soc0 + 5*1 = [15, 25]
    np.testing.assert_allclose(result[1, 0, 0], 15.0, rtol=1e-6)
    np.testing.assert_allclose(result[1, 0, 1], 25.0, rtol=1e-6)
    # Year 2: cap = soc0 + 5*2 = [20, 30]
    np.testing.assert_allclose(result[2, 0, 0], 20.0, rtol=1e-6)
    np.testing.assert_allclose(result[2, 0, 1], 30.0, rtol=1e-6)


def test_clip_soc_output_preserves_nan():
    soc0 = np.array([[10.0, np.nan]], dtype="float32")
    soc_annual = np.array([
        [[10.0, np.nan]],
        [[100.0, np.nan]],
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0, max_annual_gain=5.0,
                              global_percentile_cap=None)

    np.testing.assert_allclose(result[1, 0, 0], 15.0, rtol=1e-6)
    assert np.isnan(result[1, 0, 1])


def test_clip_soc_output_no_clip_when_within_bounds():
    soc0 = np.array([[10.0]], dtype="float32")
    soc_annual = np.array([
        [[10.0]],
        [[12.0]],  # gain of 2, within 5/year
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0, max_annual_gain=5.0,
                              global_percentile_cap=None)
    np.testing.assert_allclose(result[1, 0, 0], 12.0, rtol=1e-6)


# ---------------------------------------------------------------------------
# Tests for filter_output_raster
# ---------------------------------------------------------------------------

def test_filter_output_raster_abs_bounds(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    data = np.array([[[-5.0, 10.0, 100.0]]], dtype="float32")
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, data, transform)

    filter_output_raster(inp, out, abs_min=0.0, abs_max=50.0)

    with rasterio.open(out) as src:
        result = src.read()
    np.testing.assert_allclose(result[0, 0, 0], 0.0, rtol=1e-6)
    np.testing.assert_allclose(result[0, 0, 1], 10.0, rtol=1e-6)
    np.testing.assert_allclose(result[0, 0, 2], 50.0, rtol=1e-6)


def test_filter_output_raster_percentile(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    rng = np.random.default_rng(42)
    vals = rng.normal(10.0, 2.0, size=(1, 10, 10)).astype("float32")
    vals[0, 0, 0] = 1000.0  # extreme outlier
    vals[0, 9, 9] = -500.0  # extreme low
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, vals, transform)

    filter_output_raster(inp, out, percentile_bounds=(1.0, 99.0))

    with rasterio.open(out) as src:
        result = src.read()
    assert result[0, 0, 0] < 1000.0
    assert result[0, 9, 9] > -500.0


def test_filter_output_raster_max_annual_gain(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    # 3 bands (years): baseline=10, year1=20 (gain=10), year2=50 (gain=40)
    data = np.array([
        [[10.0, 10.0]],
        [[20.0, 20.0]],
        [[50.0, 50.0]],
    ], dtype="float32")
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, data, transform)

    filter_output_raster(inp, out, max_annual_gain=5.0, baseline_band=0)

    with rasterio.open(out) as src:
        result = src.read()
    # Band 0 unchanged
    np.testing.assert_allclose(result[0], 10.0, rtol=1e-6)
    # Band 1: cap = 10 + 5*1 = 15
    np.testing.assert_allclose(result[1], 15.0, rtol=1e-6)
    # Band 2: cap = 10 + 5*2 = 20
    np.testing.assert_allclose(result[2], 20.0, rtol=1e-6)


def test_filter_output_raster_overwrite_requires_flag(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    data = np.array([[[1.0]]], dtype="float32")
    inp = tmp_path / "input.tif"
    _write_multiband_raster(inp, data, transform)

    try:
        filter_output_raster(inp, None, abs_min=0.0)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

    # With overwrite=True it should work
    filter_output_raster(inp, None, abs_min=0.0, overwrite=True)
