import numpy as np
import numpy.testing as npt
import xarray as xr
import geopandas as gpd
import polars as pl
import rasterio
from rasterio.transform import from_origin

import sbtn_leaf.RothC_Raster as rr

_ORIGINAL_READ_FILE = gpd.read_file


def _safe_read_file(path, *args, **kwargs):
    if isinstance(path, str) and path.startswith("../data/"):
        return gpd.GeoDataFrame()
    return _ORIGINAL_READ_FILE(path, *args, **kwargs)


gpd.read_file = _safe_read_file

_ORIGINAL_READ_EXCEL = pl.read_excel


def _safe_read_excel(source, *args, **kwargs):
    if isinstance(source, str) and source.endswith("forest_residues_IPCC.xlsx"):
        return pl.DataFrame(
            {
                "IPCC Climate": ["Temperate"],
                "BD_mean": [1.0],
                "NE_mean": [1.0],
                "BD_TP": [20],
                "NE_TP": [20],
            }
        )
    if isinstance(source, str) and source.endswith("grassland_residues_IPCC.xlsx"):
        return pl.DataFrame(
            {
                "FAO_ID": [1],
                "Residue": [0.0],
                "ResXErr": [0.0],
            }
        )
    if isinstance(source, str) and source.endswith("Animals_Dung_IPCC.xlsx"):
        return pl.DataFrame(
            {
                "region": ["Western Europe"],
                "cattle_other": [0.0],
                "cattle_dairy": [0.0],
                "goat": [0.0],
                "sheep": [0.0],
            }
        )
    return _ORIGINAL_READ_EXCEL(source, *args, **kwargs)


pl.read_excel = _safe_read_excel

def _write_raster(path, data, *, nodata=None, dtype="float32"):
    if data.ndim == 2:
        count = 1
    elif data.ndim == 3:
        count = data.shape[0]
    else:
        raise ValueError("Raster data must be 2-D or 3-D")

    height, width = data.shape[-2:]
    transform = from_origin(0, height, 1, 1)

    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": count,
        "dtype": dtype,
        "crs": "EPSG:4326",
        "transform": transform,
    }
    if nodata is not None:
        profile["nodata"] = nodata

    data_to_write = data.astype(dtype)

    with rasterio.open(path, "w", **profile) as dst:
        if count == 1:
            dst.write(data_to_write, 1)
        else:
            dst.write(data_to_write)

from sbtn_leaf.RothC_Raster import (
    _load_forest_data,
    raster_rothc_annual_results,
    run_RothC_forest,
)
import sbtn_leaf.cropcalcs as cropcalcs


def test_reduced_tillage_trm_uses_full_soc(monkeypatch):
    n_years = 1
    y = x = 2
    months = 12

    clay = np.full((y, x), 30.0, dtype=float)
    soc0 = np.full((y, x), 50.0, dtype=float)
    tmp = np.full((months, y, x), 15.0, dtype=float)
    rain = np.full((months, y, x), 80.0, dtype=float)
    evap = np.full((months, y, x), 20.0, dtype=float)
    pc = np.ones((months, y, x), dtype=int)
    sand = np.full((y, x), 40.0, dtype=float)

    call_counter = {"count": 0}

    def fake_trm(sand_arr, soc_arr):
        call_counter["count"] += 1
        assert sand_arr.shape == (y, x)
        assert soc_arr.shape == (y, x)
        # return neutral modifiers so dynamics remain stable
        ones = np.ones_like(sand_arr, dtype=float)
        return ones, ones, ones, ones

    monkeypatch.setattr(
        "sbtn_leaf.RothC_Raster.RMF_TRM",
        fake_trm,
    )

    soc_annual, co2_annual = raster_rothc_annual_results(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        commodity_type="annual_crop",
        sand=sand,
        red_till=True,
    )

    # TRM should be applied once per monthly timestep (n_years * 12)
    assert call_counter["count"] == n_years * months

    # Output shapes should remain consistent with expectations
    assert soc_annual.shape == (n_years + 1, y, x)
    assert co2_annual.shape == (n_years + 1, y, x)


def test_reduced_tillage_trm_accepts_time_varying_sand(monkeypatch):
    n_years = 1
    y = x = 2
    months = 12

    clay = np.full((y, x), 30.0, dtype=float)
    soc0 = np.full((y, x), 50.0, dtype=float)
    tmp = np.full((months, y, x), 15.0, dtype=float)
    rain = np.full((months, y, x), 80.0, dtype=float)
    evap = np.full((months, y, x), 20.0, dtype=float)
    pc = np.ones((months, y, x), dtype=int)
    sand = np.full((months, y, x), 40.0, dtype=float)

    captured_shapes = []

    def fake_trm(sand_arr, soc_arr):
        captured_shapes.append((sand_arr.shape, soc_arr.shape))
        ones = np.ones_like(sand_arr, dtype=float)
        return ones, ones, ones, ones

    monkeypatch.setattr(
        "sbtn_leaf.RothC_Raster.RMF_TRM",
        fake_trm,
    )

    raster_rothc_annual_results(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        commodity_type="annual_crop",
        sand=sand,
        red_till=True,
    )

    # Every call should receive matching 2-D slices for sand and the full SOC state
    assert len(captured_shapes) == n_years * months
    assert all(s == ((y, x), (y, x)) for s in captured_shapes)


def test_unified_wrapper_accepts_shared_args(monkeypatch):
    captured_kwargs = []

    def fake_runner(**kwargs):
        captured_kwargs.append(kwargs)
        return "soc", "co2"

    monkeypatch.setattr(rr, "_raster_rothc_annual_results", fake_runner)

    months = 12
    clay = np.full((1, 1), 25.0, dtype=float)
    soc0 = np.full((1, 1), 40.0, dtype=float)
    tmp = np.full((months, 1, 1), 15.0, dtype=float)
    rain = np.full((months, 1, 1), 80.0, dtype=float)
    evap = np.full((months, 1, 1), 20.0, dtype=float)
    pc = np.ones((months, 1, 1), dtype=float)
    irr = np.zeros_like(tmp)
    c_inp = np.ones_like(tmp)
    fym = np.ones_like(tmp)
    sand = np.full((1, 1), 30.0, dtype=float)

    shared_kwargs = dict(
        n_years=1,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        irr=irr,
        c_inp=c_inp,
        fym=fym,
        crop_name="maize",
        spam_crop_raster="spam",
        practices_string_id="practice",
        irr_yield_scaling="scale",
        spam_all_fp="all",
        spam_irr_fp="irr_fp",
        spam_rf_fp="rf_fp",
        commodity_lu_fp="lu.tif",
        commodity_type="annual_crop",
        residue_runs=2,
    )

    baseline_soc, baseline_co2 = raster_rothc_annual_results(**shared_kwargs)
    assert (baseline_soc, baseline_co2) == ("soc", "co2")
    assert captured_kwargs[-1]["trm_handler"] is None
    assert captured_kwargs[-1]["sand"] is None

    captured_kwargs.clear()
    reduced_soc, reduced_co2 = raster_rothc_annual_results(
        **shared_kwargs, red_till=True, sand=sand
    )
    assert (reduced_soc, reduced_co2) == ("soc", "co2")
    assert captured_kwargs[-1]["trm_handler"] is rr.RMF_TRM
    assert captured_kwargs[-1]["sand"] is sand


def test_raster_rothc_baseline_regression(monkeypatch):
    monkeypatch.setattr("sbtn_leaf.RothC_Raster.trange", lambda n, **_: range(n))

    n_years = 2
    y = x = 2
    months = 12

    clay = np.full((y, x), 25.0, dtype=float)
    soc0 = np.full((y, x), 40.0, dtype=float)
    tmp = np.full((months, y, x), 15.0, dtype=float)
    rain = np.full((months, y, x), 80.0, dtype=float)
    evap = np.full((months, y, x), 20.0, dtype=float)
    pc = np.ones((months, y, x), dtype=float)
    irr = np.zeros_like(tmp)

    expected_soc = np.array(
        [
            [[40.0, 40.0], [40.0, 40.0]],
            [[37.534622, 37.534622], [37.534622, 37.534622]],
            [[35.909515, 35.909515], [35.909515, 35.909515]],
        ],
        dtype=np.float32,
    )
    expected_co2 = np.array(
        [
            [[0.0, 0.0], [0.0, 0.0]],
            [[2.465377, 2.465377], [2.465377, 2.465377]],
            [[1.6251065, 1.6251065], [1.6251065, 1.6251065]],
        ],
        dtype=np.float32,
    )

    soc, co2 = raster_rothc_annual_results(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        irr=irr,
        commodity_type="annual_crop",
    )

    npt.assert_allclose(soc, expected_soc)
    npt.assert_allclose(co2, expected_co2)


def test_raster_rothc_reduced_tillage_regression(monkeypatch):
    monkeypatch.setattr("sbtn_leaf.RothC_Raster.trange", lambda n, **_: range(n))

    n_years = 2
    y = x = 2
    months = 12

    clay = np.full((y, x), 25.0, dtype=float)
    soc0 = np.full((y, x), 40.0, dtype=float)
    tmp = np.full((months, y, x), 15.0, dtype=float)
    rain = np.full((months, y, x), 80.0, dtype=float)
    evap = np.full((months, y, x), 20.0, dtype=float)
    pc = np.ones((months, y, x), dtype=float)
    irr = np.zeros_like(tmp)
    sand = np.full((y, x), 30.0, dtype=float)

    expected_soc = np.array(
        [
            [[40.0, 40.0], [40.0, 40.0]],
            [[37.60099, 37.60099], [37.60099, 37.60099]],
            [[36.023407, 36.023407], [36.023407, 36.023407]],
        ],
        dtype=np.float32,
    )
    expected_co2 = np.array(
        [
            [[0.0, 0.0], [0.0, 0.0]],
            [[2.3990104, 2.3990104], [2.3990104, 2.3990104]],
            [[1.5775822, 1.5775822], [1.5775822, 1.5775822]],
        ],
        dtype=np.float32,
    )

    soc, co2 = raster_rothc_annual_results(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        irr=irr,
        commodity_type="annual_crop",
        sand=sand,
        red_till=True,
    )

    npt.assert_allclose(soc, expected_soc)
    npt.assert_allclose(co2, expected_co2)


def test_run_rothc_forest_handles_single_band_age(monkeypatch, tmp_path):
    monkeypatch.setattr("sbtn_leaf.RothC_Raster.trange", lambda n, **_: range(n))

    n_years = 2
    y = x = 2
    months = 12

    coords = {
        "time": np.arange(months),
        "y": np.arange(y),
        "x": np.arange(x),
    }

    tmp = xr.DataArray(np.full((months, y, x), 12.0, dtype=float), dims=("time", "y", "x"), coords=coords)
    rain = xr.full_like(tmp, 80.0)
    evap = xr.full_like(tmp, 20.0)
    pc = xr.full_like(tmp, 1.0)

    soc0 = xr.DataArray(np.full((y, x), 50.0, dtype=float), dims=("y", "x"))
    clay = xr.full_like(soc0, 30.0)
    iom = xr.full_like(soc0, 5.0)
    sand = xr.full_like(soc0, 40.0)
    lu = xr.full_like(soc0, 1.0)
    evap_da = xr.DataArray(np.full((months, y, x), 20.0, dtype=float), dims=("time", "y", "x"), coords=coords)
    pc_da = xr.full_like(tmp, 1.0)
    age = xr.DataArray(np.full((1, y, x), 10.0, dtype=float), dims=("band", "y", "x"))

    monkeypatch.setattr(
        "sbtn_leaf.RothC_Raster._load_environmental_data",
        lambda *_, **__: (tmp, rain, soc0, iom, clay, sand),
    )
    monkeypatch.setattr(
        "sbtn_leaf.RothC_Raster._load_forest_data",
        lambda *_: (lu, evap_da, pc_da, age),
    )
    monkeypatch.setattr("sbtn_leaf.RothC_Raster.save_annual_results", lambda *_, **__: None)

    def fake_forest_litter(age_arr, *_, **__):
        y_dim, x_dim = age_arr.shape
        return np.full((1, y_dim, x_dim), 0.5, dtype=float)

    monkeypatch.setattr(
        "sbtn_leaf.RothC_Raster.cropcalcs.get_forest_litter_monthlyrate_fromda",
        fake_forest_litter,
    )

    soc = run_RothC_forest(
        forest_type="BRDC",
        weather_type="Temperate",
        n_years=n_years,
        save_folder=str(tmp_path),
        data_description="test",
        lu_fp="lu.tif",
        evap_fp="evap.tif",
        age_fp="age.tif",
    )

    assert soc.shape == (n_years + 1, y, x)
    assert np.isfinite(soc).all()


def test_load_forest_data_preserves_age_nans(tmp_path):
    y = x = 2

    lu_data = np.array([[1, 0], [1, 1]], dtype=np.float32)
    evap_data = np.full((12, y, x), 2.0, dtype=np.float32)
    age_data = np.array([[10.0, 20.0], [-9999.0, 30.0]], dtype=np.float32)

    lu_path = tmp_path / "lu.tif"
    evap_path = tmp_path / "evap.tif"
    age_path = tmp_path / "age.tif"

    _write_raster(lu_path, lu_data)
    _write_raster(evap_path, evap_data)
    _write_raster(age_path, age_data, nodata=-9999.0)

    _, _, _, age = _load_forest_data(str(lu_path), str(evap_path), str(age_path))

    age_arr = age.values
    assert age_arr.shape == (y, x)
    assert age_arr[0, 0] == 10.0
    assert age_arr[1, 1] == 30.0
    assert np.isnan(age_arr[0, 1])
    assert np.isnan(age_arr[1, 0])


def test_run_rothc_forest_passes_nan_age_through_litter(monkeypatch, tmp_path):
    monkeypatch.setattr("sbtn_leaf.RothC_Raster.trange", lambda n, **_: range(n))

    y = x = 2
    months = 12
    n_years = 1

    lu_data = np.array([[1, 0], [1, 1]], dtype=np.float32)
    evap_data = np.full((12, y, x), 2.0, dtype=np.float32)
    age_data = np.array([[10.0, 20.0], [-9999.0, 30.0]], dtype=np.float32)

    lu_path = tmp_path / "lu.tif"
    evap_path = tmp_path / "evap.tif"
    age_path = tmp_path / "age.tif"

    _write_raster(lu_path, lu_data)
    _write_raster(evap_path, evap_data)
    _write_raster(age_path, age_data, nodata=-9999.0)

    coords = {"time": np.arange(months), "y": np.arange(y), "x": np.arange(x)}
    tmp = xr.DataArray(np.full((months, y, x), 12.0, dtype=float), dims=("time", "y", "x"), coords=coords)
    rain = xr.full_like(tmp, 80.0)
    evap_env = xr.full_like(tmp, 20.0)
    pc = xr.full_like(tmp, 1.0)
    soc0 = xr.DataArray(np.full((y, x), 50.0, dtype=float), dims=("y", "x"))
    clay = xr.full_like(soc0, 30.0)
    iom = xr.full_like(soc0, 5.0)
    sand = xr.full_like(soc0, 40.0)

    monkeypatch.setattr(
        "sbtn_leaf.RothC_Raster._load_environmental_data",
        lambda *_, **__: (tmp, rain, soc0, iom, clay, sand),
    )
    monkeypatch.setattr("sbtn_leaf.RothC_Raster.save_annual_results", lambda *_, **__: None)

    captured_inputs = []
    captured_outputs = []

    original_litter = cropcalcs.get_forest_litter_monthlyrate_fromda

    def capture_litter(age_arr, *args, **kwargs):
        captured_inputs.append(np.copy(age_arr))
        result = original_litter(age_arr, *args, **kwargs)
        captured_outputs.append(np.copy(result))
        return result

    monkeypatch.setattr(
        "sbtn_leaf.RothC_Raster.cropcalcs.get_forest_litter_monthlyrate_fromda",
        capture_litter,
    )

    run_RothC_forest(
        forest_type="BRDC",
        weather_type="Temperate",
        n_years=n_years,
        save_folder=str(tmp_path),
        data_description="test",
        lu_fp=str(lu_path),
        evap_fp=str(evap_path),
        age_fp=str(age_path),
    )

    assert captured_inputs
    for arr in captured_inputs:
        assert np.isnan(arr[0, 1])
        assert np.isnan(arr[1, 0])
        assert arr[0, 0] == 10.0
        assert arr[1, 1] == 30.0

    assert captured_outputs
    for arr in captured_outputs:
        assert np.isnan(arr[0, 1])
        assert np.isnan(arr[1, 0])


# ---------------------------------------------------------------------------
# Tests for _clip_soc_output post-filter
# ---------------------------------------------------------------------------

from sbtn_leaf.RothC_Raster import _clip_soc_output


def test_clip_soc_output_absolute_cap():
    """Values exceeding max_soc_tc_ha should be clipped."""
    soc0 = np.array([[40.0, 50.0]], dtype="float32")
    soc_annual = np.array([
        [[40.0, 50.0]],
        [[100.0, 600.0]],  # 600 exceeds default 500 cap
    ], dtype="float32")

    result = _clip_soc_output(
        soc_annual, soc0, max_annual_gain=1000.0,  # high to not trigger delta cap
        max_soc_tc_ha=500.0, global_percentile_cap=None,
    )

    assert result[1, 0, 0] == 100.0  # unchanged, below cap
    assert result[1, 0, 1] == 500.0  # clipped to cap


def test_clip_soc_output_annual_delta_cap():
    """Year-over-year SOC gain should be capped at soc0 + max_annual_gain * year."""
    soc0 = np.array([[40.0]], dtype="float32")
    soc_annual = np.array([
        [[40.0]],
        [[60.0]],  # gain of 20 from soc0
        [[90.0]],  # gain of 50 from soc0
    ], dtype="float32")

    result = _clip_soc_output(
        soc_annual, soc0, max_annual_gain=5.0,
        max_soc_tc_ha=1000.0, global_percentile_cap=None,
    )

    # year 1: cap = 40 + 5*1 = 45
    assert result[1, 0, 0] == 45.0
    # year 2: cap = 40 + 5*2 = 50
    assert result[2, 0, 0] == 50.0


def test_clip_soc_output_global_percentile():
    """Global percentile cap should clip high values across final year."""
    soc0 = np.full((1, 100), 40.0, dtype="float32")
    year1 = np.full((1, 100), 50.0, dtype="float32")
    year1[0, 99] = 200.0  # single outlier

    soc_annual = np.stack([soc0, year1])

    result = _clip_soc_output(
        soc_annual, soc0, max_annual_gain=1000.0,
        max_soc_tc_ha=1000.0, global_percentile_cap=99.0,
    )

    # The outlier at 200.0 should be clipped down
    assert result[1, 0, 99] < 200.0


def test_clip_soc_output_preserves_nans():
    """NaN pixels should be left unchanged."""
    soc0 = np.array([[40.0, np.nan]], dtype="float32")
    soc_annual = np.array([
        [[40.0, np.nan]],
        [[60.0, np.nan]],
    ], dtype="float32")

    result = _clip_soc_output(
        soc_annual, soc0, max_annual_gain=5.0,
        max_soc_tc_ha=500.0, global_percentile_cap=None,
    )

    assert np.isnan(result[0, 0, 1])
    assert np.isnan(result[1, 0, 1])
    assert result[1, 0, 0] == 45.0  # 40 + 5*1
