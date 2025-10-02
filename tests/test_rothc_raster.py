import numpy as np
import numpy.testing as npt

from sbtn_leaf.RothC_Raster import (
    raster_rothc_annual_results_1yrloop,
    raster_rothc_ReducedTillage_annual_results_1yrloop,
)


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

    soc_annual, co2_annual = raster_rothc_ReducedTillage_annual_results_1yrloop(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        sand=sand,
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

    raster_rothc_ReducedTillage_annual_results_1yrloop(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        sand=sand,
    )

    # Every call should receive matching 2-D slices for sand and the full SOC state
    assert len(captured_shapes) == n_years * months
    assert all(s == ((y, x), (y, x)) for s in captured_shapes)


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

    soc, co2 = raster_rothc_annual_results_1yrloop(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        irr=irr,
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

    soc, co2 = raster_rothc_ReducedTillage_annual_results_1yrloop(
        n_years=n_years,
        clay=clay,
        soc0=soc0,
        tmp=tmp,
        rain=rain,
        evap=evap,
        pc=pc,
        sand=sand,
        irr=irr,
    )

    npt.assert_allclose(soc, expected_soc)
    npt.assert_allclose(co2, expected_co2)
