import numpy as np

from sbtn_leaf.RothC_Core import RMF_TRM, run_simulation


def test_run_simulation_defaults_no_optional_inputs():
    n_years = 1
    months = n_years * 12
    soc_annual, co2_annual = run_simulation(
        clay=30.0,
        depth=23.0,
        iom=5.0,
        soc0=50.0,
        tmp=np.full(months, 10.0),
        rain=np.full(months, 50.0),
        evap=np.full(months, 20.0),
        pc=np.zeros(months, dtype=int),
        dpm_rpm=1.44,
        n_years=n_years,
        do_equilibrium=True,
    )

    assert soc_annual.shape == (n_years,)
    assert co2_annual.shape == (n_years,)
    assert np.isfinite(soc_annual).all()
    assert np.isfinite(co2_annual).all()


def test_rmf_trm_preserves_modifier_values():
    sand = np.array([[10.0, 36.0], [38.0, 20.0]])
    soc = np.array([[10.0, 80.0], [70.0, 100.0]])

    trm_dpm, trm_rpm, trm_bio, trm_hum = RMF_TRM(sand, soc)

    expected_dpm = np.array([[0.72, 1.54], [0.72, 0.72]])
    expected_rpm = np.array([[0.97, 2.15], [0.97, 0.97]])
    expected_bio = np.array([[0.99, 2.38], [0.99, 0.99]])
    expected_hum = np.array([[0.94, 2.93], [0.94, 0.94]])

    np.testing.assert_allclose(trm_dpm, expected_dpm)
    np.testing.assert_allclose(trm_rpm, expected_rpm)
    np.testing.assert_allclose(trm_bio, expected_bio)
    np.testing.assert_allclose(trm_hum, expected_hum)
