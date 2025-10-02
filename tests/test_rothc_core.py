import numpy as np

from sbtn_leaf.RothC_Core import run_simulation


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
