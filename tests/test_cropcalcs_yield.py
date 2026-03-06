import numpy as np
import geopandas as gpd
from shapely.geometry import box
from rasterio.transform import from_origin
import rasterio
import sbtn_leaf.cropcalcs as cropcalcs

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

    fao_gdf = _make_fao_gdf("avg_yield", "yld_ratio", avg_value=1000.0, ratio_value=0.5)

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


def test_fill_with_ecoregions_uses_global_yield_fallback_scale(monkeypatch):
    result = np.array([[np.nan, np.nan]], dtype=float)
    lu_mask = np.array([[True, True]])

    def fake_ecoregion_stats(_result, _croplu):
        zone_array = np.array([[-1, -1]], dtype=int)
        return {}, {}, zone_array, {}

    monkeypatch.setattr(cropcalcs, "calculate_average_yield_by_ecoregion_and_biome", fake_ecoregion_stats)

    filled, _ = cropcalcs._fill_with_ecoregions(
        result.copy(),
        "dummy.tif",
        lu_mask,
        global_fao_yield_fallback=3.75,
        enable_ecoregion_fill=True,
        enable_nearest_fill=False,
    )

    np.testing.assert_allclose(filled[lu_mask], 3.75, rtol=1e-6, atol=1e-6)
    assert not np.any(np.isclose(filled[lu_mask], 1.0))


def test_fill_with_ecoregions_prioritizes_ecoregion_then_biome_before_global(monkeypatch):
    result = np.array([[5.0, np.nan, np.nan]], dtype=float)
    lu_mask = np.array([[True, True, True]])

    def fake_ecoregion_stats(_result, _croplu):
        ecoregion_avg = {0: 5.0}
        biome_avg = {"BiomeB": 2.5}
        zone_array = np.array([[0, 1, 2]], dtype=int)
        biome_name_map = {1: "BiomeB"}
        return ecoregion_avg, biome_avg, zone_array, biome_name_map

    monkeypatch.setattr(cropcalcs, "calculate_average_yield_by_ecoregion_and_biome", fake_ecoregion_stats)

    filled, _ = cropcalcs._fill_with_ecoregions(
        result.copy(),
        "dummy.tif",
        lu_mask,
        global_fao_yield_fallback=7.0,
        enable_ecoregion_fill=True,
        enable_nearest_fill=False,
    )

    np.testing.assert_allclose(filled[0, 1], 2.5, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(filled[0, 2], 7.0, rtol=1e-6, atol=1e-6)



def test_pre_filter_local_zscore_removes_isolated_spike():
    arr = np.ones((9, 9), dtype="float32")
    arr[4, 4] = 100.0
    arr[0, 0] = np.nan

    zone_array = np.ones((9, 9), dtype=int)
    zone_array[0, 8] = 0

    fao_gdf = gpd.GeoDataFrame({"zone_id": [1], "avg_yield": [10.0]}, geometry=[box(0, 0, 1, 1)], crs="EPSG:4326")

    (filtered,) = cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="local_zscore",
        percentile_bounds=(1.0, 99.0),
        k_sd=1000.0,
        local_window=5,
        local_k=1.0,
        local_min_neighbors=8,
    )

    assert np.isnan(filtered[0, 0])
    assert np.isnan(filtered[0, 8])
    assert filtered[4, 4] < 30.0


def test_pre_filter_local_zscore_preserves_local_gradient():
    arr = np.tile(np.linspace(1, 20, 21, dtype="float32"), (21, 1))
    zone_array = np.ones_like(arr, dtype=int)
    fao_gdf = gpd.GeoDataFrame({"zone_id": [1], "avg_yield": [10.0]}, geometry=[box(0, 0, 1, 1)], crs="EPSG:4326")

    (filtered,) = cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="local_zscore",
        percentile_bounds=(1.0, 99.0),
        k_sd=1000.0,
        local_window=3,
        local_k=3.0,
        local_min_neighbors=3,
    )

    np.testing.assert_allclose(filtered, arr, rtol=1e-6, atol=1e-6)


def test_pre_filter_local_zscore_runtime_sanity():
    rng = np.random.default_rng(42)
    arr = rng.normal(10.0, 2.0, size=(512, 512)).astype("float32")
    arr[rng.random(arr.shape) < 0.1] = np.nan
    zone_array = np.ones(arr.shape, dtype=int)
    fao_gdf = gpd.GeoDataFrame({"zone_id": [1], "avg_yield": [10.0]}, geometry=[box(0, 0, 1, 1)], crs="EPSG:4326")

    import time
    t0 = time.perf_counter()
    cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="local_zscore",
        percentile_bounds=(1.0, 99.0),
        k_sd=1000.0,
        local_window=5,
        local_k=2.5,
        local_min_neighbors=8,
    )
    elapsed = time.perf_counter() - t0

    assert elapsed < 5.0


def test_apply_uncertainty_to_monthly_residues_preserves_sparse_zero_months():
    monthly = np.zeros((12, 2, 2), dtype="float32")
    monthly[2] = np.array([[1.0, 0.0], [2.0, 0.0]], dtype="float32")
    monthly[7] = np.array([[0.5, 0.0], [1.5, 0.0]], dtype="float32")

    lu_mask = np.array([[True, False], [True, False]])
    fao_avg = np.array([[3.0, 0.0], [4.0, 0.0]], dtype="float32")
    fao_sd = np.array([[0.3, 0.0], [0.4, 0.0]], dtype="float32")

    randomized = cropcalcs._apply_uncertainty_to_monthly_residues(
        monthly_residues=monthly,
        fao_avg_yields_array=fao_avg,
        fao_sd_yields_array=fao_sd,
        lu_mask=lu_mask,
        random_runs=20,
        rng=np.random.default_rng(0),
    )

    assert randomized.shape == monthly.shape
    np.testing.assert_allclose(randomized[monthly == 0.0], 0.0, rtol=0.0, atol=0.0)


# ---------------------------------------------------------------------------
# Tests for fao_max_ratio yield cap
# ---------------------------------------------------------------------------

def test_fao_max_ratio_clips_high_yields_sd_strategy():
    """Yields exceeding fao_max_ratio * FAO avg should be clipped."""
    arr = np.array([[5.0, 50.0], [100.0, 200.0]], dtype="float32")
    zone_array = np.ones((2, 2), dtype=int)
    fao_gdf = gpd.GeoDataFrame(
        {"zone_id": [1], "avg_yield": [10.0]},
        geometry=[box(0, 0, 1, 1)],
        crs="EPSG:4326",
    )

    (filtered,) = cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="sd",
        k_sd=100.0,  # very wide — would not clip on its own
        fao_max_ratio=3.0,  # cap at 30.0
    )

    # All pixels should be <= 3.0 * 10.0 = 30.0
    assert np.nanmax(filtered) <= 30.0 + 1e-6
    np.testing.assert_allclose(filtered[0, 0], 5.0, rtol=1e-6)
    np.testing.assert_allclose(filtered[0, 1], 30.0, rtol=1e-6)
    np.testing.assert_allclose(filtered[1, 0], 30.0, rtol=1e-6)
    np.testing.assert_allclose(filtered[1, 1], 30.0, rtol=1e-6)


def test_fao_max_ratio_none_disables_cap():
    """When fao_max_ratio is None, no FAO-based cap is applied."""
    arr = np.array([[5.0, 200.0]], dtype="float32")
    zone_array = np.ones((1, 2), dtype=int)
    fao_gdf = gpd.GeoDataFrame(
        {"zone_id": [1], "avg_yield": [10.0]},
        geometry=[box(0, 0, 1, 1)],
        crs="EPSG:4326",
    )

    (filtered,) = cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="sd",
        k_sd=100.0,
        fao_max_ratio=None,  # disabled
    )

    # 200.0 should pass through (only sd strategy with huge k_sd)
    np.testing.assert_allclose(filtered[0, 1], 200.0, rtol=1e-6)


def test_fao_max_ratio_with_local_zscore_strategy():
    """fao_max_ratio should also cap values when using local_zscore strategy."""
    arr = np.ones((5, 5), dtype="float32") * 50.0
    zone_array = np.ones((5, 5), dtype=int)
    fao_gdf = gpd.GeoDataFrame(
        {"zone_id": [1], "avg_yield": [10.0]},
        geometry=[box(0, 0, 1, 1)],
        crs="EPSG:4326",
    )

    (filtered,) = cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="local_zscore",
        fao_max_ratio=2.0,  # cap at 20.0
        local_window=3,
        local_k=2.5,
        local_min_neighbors=4,
    )

    # All values should be capped at 2.0 * 10.0 = 20.0
    finite = filtered[np.isfinite(filtered)]
    assert np.all(finite <= 20.0 + 1e-6)


# ---------------------------------------------------------------------------
# Tests for global_percentile_cap
# ---------------------------------------------------------------------------

def test_global_percentile_cap_clips_extreme_values():
    """global_percentile_cap should clip values above the computed percentile."""
    # 100 pixels: 99 at value 10, 1 at value 1000
    arr = np.full((10, 10), 10.0, dtype="float32")
    arr[0, 0] = 1000.0
    zone_array = np.ones((10, 10), dtype=int)
    fao_gdf = gpd.GeoDataFrame(
        {"zone_id": [1], "avg_yield": [10.0]},
        geometry=[box(0, 0, 1, 1)],
        crs="EPSG:4326",
    )

    (filtered,) = cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="sd",
        k_sd=100.0,  # very wide, won't clip on its own
        fao_max_ratio=None,  # disable FAO cap
        global_percentile_cap=99.0,  # should clip the single outlier
    )

    # The 99th percentile of 99x10 + 1x1000 should bring the outlier down
    assert np.nanmax(filtered) < 1000.0


def test_global_percentile_cap_none_disables():
    """When global_percentile_cap is None, no global clipping occurs."""
    arr = np.full((10, 10), 10.0, dtype="float32")
    arr[0, 0] = 1000.0
    zone_array = np.ones((10, 10), dtype=int)
    fao_gdf = gpd.GeoDataFrame(
        {"zone_id": [1], "avg_yield": [10.0]},
        geometry=[box(0, 0, 1, 1)],
        crs="EPSG:4326",
    )

    (filtered,) = cropcalcs._pre_filter_yields_rasters(
        yld_arrays=(arr,),
        fao_gdf=fao_gdf,
        zone_array=zone_array,
        fao_avg_yield_name="avg_yield",
        filter_outlier_strategy="sd",
        k_sd=100.0,
        fao_max_ratio=None,
        global_percentile_cap=None,
    )

    # Outlier should pass through untouched
    np.testing.assert_allclose(filtered[0, 0], 1000.0, rtol=1e-6)
