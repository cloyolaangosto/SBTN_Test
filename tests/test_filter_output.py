"""Tests for _clip_soc_output (RothC_Raster), filter_output_raster, and filter_raster_outliers (map_calculations)."""

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from sbtn_leaf.RothC_Raster import _clip_soc_output
from sbtn_leaf.map_calculations import (
    filter_output_raster,
    filter_raster_outliers,
    _apply_log_winsor,
    _apply_local_zscore,
)


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


def test_clip_soc_output_abs_cap_disabled():
    """When max_soc_tc_ha=None, values above 500 should pass through."""
    soc0 = np.array([[40.0]], dtype="float32")
    soc_annual = np.array([
        [[40.0]],
        [[600.0]],  # exceeds default 500 cap
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0,
                              max_soc_tc_ha=None,
                              max_annual_gain=None,
                              global_percentile_cap=None)
    # With all filters disabled, value should be unchanged
    np.testing.assert_allclose(result[1, 0, 0], 600.0, rtol=1e-6)


def test_clip_soc_output_delta_cap_disabled():
    """When max_annual_gain=None, large year-over-year gains pass through."""
    soc0 = np.array([[10.0]], dtype="float32")
    soc_annual = np.array([
        [[10.0]],
        [[100.0]],  # gain of 90 in 1 year
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0,
                              max_soc_tc_ha=None,
                              max_annual_gain=None,
                              global_percentile_cap=None)
    np.testing.assert_allclose(result[1, 0, 0], 100.0, rtol=1e-6)


def test_clip_soc_output_all_disabled():
    """When all filters are None, output equals input for finite values."""
    soc0 = np.array([[10.0, 20.0]], dtype="float32")
    soc_annual = np.array([
        [[10.0, 20.0]],
        [[600.0, 800.0]],
        [[1000.0, 2000.0]],
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0,
                              max_soc_tc_ha=None,
                              max_annual_gain=None,
                              global_percentile_cap=None)
    np.testing.assert_array_equal(result, soc_annual)


def test_clip_soc_output_only_abs_cap():
    """Only abs cap active; delta and percentile disabled."""
    soc0 = np.array([[10.0]], dtype="float32")
    soc_annual = np.array([
        [[10.0]],
        [[600.0]],  # exceeds 500 cap, but huge gain should pass if delta disabled
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0,
                              max_soc_tc_ha=500.0,
                              max_annual_gain=None,
                              global_percentile_cap=None)
    # Clipped by abs cap only
    np.testing.assert_allclose(result[1, 0, 0], 500.0, rtol=1e-6)


def test_clip_soc_output_only_delta_cap():
    """Only delta cap active; abs and percentile disabled."""
    soc0 = np.array([[10.0]], dtype="float32")
    soc_annual = np.array([
        [[10.0]],
        [[600.0]],  # gain of 590 in 1 year
    ], dtype="float32")

    result = _clip_soc_output(soc_annual, soc0,
                              max_soc_tc_ha=None,
                              max_annual_gain=5.0,
                              global_percentile_cap=None)
    # cap = 10 + 5*1 = 15
    np.testing.assert_allclose(result[1, 0, 0], 15.0, rtol=1e-6)


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


# ---------------------------------------------------------------------------
# Tests for _apply_log_winsor
# ---------------------------------------------------------------------------

def test_log_winsor_clips_extreme_high():
    arr = np.array([[1.0, 2.0, 3.0, 4.0, 1000.0]], dtype="float32")
    result = _apply_log_winsor(arr, (1.0, 99.0))
    assert result[0, -1] < 1000.0
    assert result[0, 2] == pytest.approx(arr[0, 2], rel=1e-5)


def test_log_winsor_clips_extreme_low():
    arr = np.array([[0.001, 5.0, 5.0, 5.0, 5.0]], dtype="float32")
    result = _apply_log_winsor(arr, (5.0, 99.0))
    assert result[0, 0] > 0.001


def test_log_winsor_preserves_nan():
    arr = np.array([[1.0, np.nan, 3.0]], dtype="float32")
    result = _apply_log_winsor(arr, (1.0, 99.0))
    assert np.isnan(result[0, 1])


def test_log_winsor_rejects_invalid_bounds():
    arr = np.ones((1, 4), dtype="float32")
    with pytest.raises(ValueError, match="log_winsor_bounds"):
        _apply_log_winsor(arr, (99.0, 1.0))


def test_log_winsor_rejects_values_le_minus1():
    arr = np.array([[-2.0, 1.0, 2.0]], dtype="float32")
    with pytest.raises(ValueError, match="> -1"):
        _apply_log_winsor(arr, (1.0, 99.0))


def test_log_winsor_all_nan_passthrough():
    arr = np.full((3, 3), np.nan, dtype="float32")
    result = _apply_log_winsor(arr, (1.0, 99.0))
    assert np.all(np.isnan(result))


# ---------------------------------------------------------------------------
# Tests for _apply_local_zscore
# ---------------------------------------------------------------------------

def test_local_zscore_clips_spike():
    rng = np.random.default_rng(0)
    arr = rng.normal(10.0, 1.0, (30, 30)).astype("float32")
    arr[15, 15] = 500.0
    result = _apply_local_zscore(arr, local_window=9, local_k=2.5, local_min_neighbors=5)
    assert result[15, 15] < 500.0


def test_local_zscore_leaves_normal_values_unchanged():
    arr = np.full((20, 20), 5.0, dtype="float32")
    result = _apply_local_zscore(arr, local_window=5, local_k=2.5, local_min_neighbors=3)
    np.testing.assert_allclose(result, arr, rtol=1e-5)


def test_local_zscore_preserves_nan():
    arr = np.full((10, 10), 5.0, dtype="float32")
    arr[3, 3] = np.nan
    result = _apply_local_zscore(arr, local_window=5, local_k=2.5, local_min_neighbors=3)
    assert np.isnan(result[3, 3])


def test_local_zscore_rejects_even_window():
    arr = np.ones((5, 5), dtype="float32")
    with pytest.raises(ValueError, match="odd"):
        _apply_local_zscore(arr, local_window=4, local_k=2.5, local_min_neighbors=3)


def test_local_zscore_all_nan_passthrough():
    arr = np.full((5, 5), np.nan, dtype="float32")
    result = _apply_local_zscore(arr, local_window=3, local_k=2.5, local_min_neighbors=2)
    assert np.all(np.isnan(result))


# ---------------------------------------------------------------------------
# Tests for filter_raster_outliers (file-level)
# ---------------------------------------------------------------------------

def test_filter_raster_outliers_log_winsor(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    rng = np.random.default_rng(7)
    vals = rng.lognormal(mean=2.0, sigma=0.5, size=(1, 10, 10)).astype("float32")
    vals[0, 0, 0] = 5000.0  # outlier
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, vals, transform)

    filter_raster_outliers(inp, out, log_winsor_bounds=(1.0, 99.0))

    with rasterio.open(out) as src:
        result = src.read()
    assert result[0, 0, 0] < 5000.0


def test_filter_raster_outliers_local_zscore(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    rng = np.random.default_rng(8)
    vals = rng.normal(10.0, 1.0, (1, 30, 30)).astype("float32")
    vals[0, 15, 15] = 500.0
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, vals, transform)

    filter_raster_outliers(
        inp, out,
        apply_local_zscore=True,
        local_window=9,
        local_k=2.5,
        local_min_neighbors=5,
    )

    with rasterio.open(out) as src:
        result = src.read()
    assert result[0, 15, 15] < 500.0


def test_filter_raster_outliers_both_filters(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    rng = np.random.default_rng(9)
    vals = rng.lognormal(mean=2.0, sigma=0.5, size=(1, 30, 30)).astype("float32")
    vals[0, 0, 0] = 5000.0
    vals[0, 15, 15] = 4000.0
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, vals, transform)

    filter_raster_outliers(
        inp, out,
        log_winsor_bounds=(0.5, 99.5),
        apply_local_zscore=True,
        local_window=9,
        local_k=2.5,
        local_min_neighbors=5,
    )

    with rasterio.open(out) as src:
        result = src.read()
    assert result[0, 0, 0] < 5000.0
    assert result[0, 15, 15] < 4000.0


def test_filter_raster_outliers_global_percentile_cap(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    vals = np.arange(100, dtype="float32").reshape(1, 10, 10)
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, vals, transform)

    filter_raster_outliers(inp, out, global_percentile_cap=90.0)

    with rasterio.open(out) as src:
        result = src.read()
    cap = np.nanpercentile(vals, 90.0)
    assert result.max() <= cap + 1e-4


def test_filter_raster_outliers_overwrite_requires_flag(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    data = np.array([[[1.0]]], dtype="float32")
    inp = tmp_path / "input.tif"
    _write_multiband_raster(inp, data, transform)

    with pytest.raises(ValueError):
        filter_raster_outliers(inp, None, log_winsor_bounds=(1.0, 99.0))

    filter_raster_outliers(inp, None, log_winsor_bounds=(1.0, 99.0), overwrite=True)


def test_filter_raster_outliers_no_filter_passthrough(tmp_path):
    transform = from_origin(0, 1, 1, 1)
    data = np.array([[[1.0, 2.0, 3.0]]], dtype="float32")
    inp = tmp_path / "input.tif"
    out = tmp_path / "output.tif"
    _write_multiband_raster(inp, data, transform)

    filter_raster_outliers(inp, out)

    with rasterio.open(out) as src:
        result = src.read()
    np.testing.assert_allclose(result, data, rtol=1e-5)
