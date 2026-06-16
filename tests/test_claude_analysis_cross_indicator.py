"""Tests for the cross-indicator aggregation analysis (``sbtn_leaf.claude_analysis``)."""

import numpy as np
import pytest

from sbtn_leaf.claude_analysis import indicator_aggregation as eng
from sbtn_leaf.claude_analysis import cross_indicator as xi
from sbtn_leaf.claude_analysis.indicators import (
    INDICATORS,
    SOC,
    SOIL_EROSION,
    canonical_flow_soc,
)


# --------------------------------------------------------------------------- #
# SOC flow canonicalisation (pure, no IO)
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("Wheat_rf_roff_2030y_SOC", "Wheat|rf|roff|ct"),
        ("Maize_rf_roff_rt_2030y_SOC", "Wheat|rf|roff|rt".replace("Wheat", "Maize")),
        ("Soybean_rf_2030y_SOC", "Soybeans|rf|na|ct"),  # singular -> plural alias
        ("Oil palm_rf_2030y_SOC", "Oil_palm|rf|na|ct"),  # space -> underscore
        ("Coffee_rf_2030y_SOC", "Coffee|rf|na|ct"),
        ("Cotton_rf_rt_2030y_SOC", "Cotton|rf|na|rt"),
        ("BRDC_Tropical_2030y_SOC", "BRDC_Tropical"),
        ("BRDC_Boreal dry_2030y_SOC", "BRDC_Boreal_dry"),
        ("NEEV_Tropical_2030y_SOC", "NEEV_Tropical"),
        ("natural_grassland_cattle_avg_2030y_SOC", "Grassland"),
        ("natural_grassland_sheep_2030y_SOC", "Grassland|sheep"),
    ],
)
def test_canonical_flow_soc(raw, expected):
    assert canonical_flow_soc(raw) == expected


def test_soc_focal_keys_match_soil_erosion_keys():
    # SOC focal keys are shared with the soil-erosion convention so the same
    # commodity lines up across indicators in the cross-indicator view.
    assert set(SOC.focal_flows) >= {"Wheat|rf|roff|ct", "Soybeans|rf|na|ct", "Grassland"}
    assert set(SOC.focal_flows) == set(SOIL_EROSION.focal_flows)


# --------------------------------------------------------------------------- #
# Harmonised loading (touch the real CSVs once per indicator, module scope)
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module", params=list(INDICATORS.values()), ids=list(INDICATORS))
def cfg(request):
    return request.param


@pytest.fixture(scope="module")
def harmonized(cfg):
    return cfg.load_harmonized(drop_na=False)


def test_schema_and_levels(cfg, harmonized):
    assert list(harmonized.columns) == eng.STD_COLS
    assert set(harmonized["level"].unique()) == set(eng.LEVELS)
    # ecoregion rows carry biome / realm; political levels do not
    eco = harmonized[harmonized["level"] == "ecoregion"]
    assert eco["biome"].notna().any()
    assert eco["realm"].notna().any()


def test_focal_flows_present_at_all_levels(cfg, harmonized):
    cov = eng.validate_flow_coverage(harmonized, focal_flows=list(cfg.focal_flows))
    for focal in cfg.focal_flows:
        assert focal in cov.index, f"{focal} absent from {cfg.key}"
        assert bool(cov.loc[focal, "in_all"]), f"{focal} missing from a level in {cfg.key}"


def test_coverage_increases_for_land_use_indicators(cfg, harmonized):
    cov = eng.coverage_table(harmonized)
    assert list(cov["level"]) == list(eng.LEVELS)
    if cfg.key == "acidification":
        # CFs are defined everywhere -> full coverage at every level
        assert (cov["coverage_pct"] > 99).all()
    else:
        # land-use indicators resolve progressively more land
        assert cov["coverage_pct"].is_monotonic_increasing


def test_summary_and_variance_bounds(cfg):
    df = cfg.load_harmonized(drop_na=True)
    summ = eng.summary_table(cfg, df, flows=list(cfg.focal_flows))
    assert {"mean_leaf", "median_leaf", "std_leaf", "cv", "indicator", "flow_label"}.issubset(summ.columns)
    assert (summ.groupby("flow", observed=True)["level"].nunique() == len(eng.LEVELS)).all()
    assert np.isfinite(summ["mean_leaf"]).all()

    var = eng.hierarchy_variance_table(cfg, df, flows=list(cfg.focal_flows))
    for col in ("subcty_within_country_frac", "eco_between_biome_eta2", "eco_between_realm_eta2"):
        vals = var[col].dropna()
        assert ((vals >= -1e-9) & (vals <= 1 + 1e-9)).all()


# --------------------------------------------------------------------------- #
# Cross-indicator tables
# --------------------------------------------------------------------------- #


def test_significance_table_one_row_per_indicator():
    sig = xi.significance_table()
    assert list(sig["indicator"]) == list(INDICATORS)
    for col in ("eco_biome_eta2", "eco_realm_eta2", "subcty_within_country_frac"):
        assert ((sig[col] >= -1e-9) & (sig[col] <= 1 + 1e-9)).all()


def test_level_reframing_country_is_unity():
    ref = xi.level_reframing_table()
    # everything is normalised to the country mean / median
    assert np.allclose(ref["rel_mean_country"], 1.0)
    assert np.allclose(ref["rel_median_country"], 1.0)
    assert {f"rel_mean_{lvl}" for lvl in eng.LEVELS}.issubset(ref.columns)


def test_dispersion_table_shape():
    disp = xi.dispersion_table()
    assert len(disp) == len(INDICATORS)
    assert {f"cv_{lvl}" for lvl in eng.LEVELS}.issubset(disp.columns)
    assert (disp[[f"cv_{lvl}" for lvl in eng.LEVELS]].to_numpy() > 0).all()


# --------------------------------------------------------------------------- #
# Plot smoke tests (Agg backend, no display)
# --------------------------------------------------------------------------- #


def test_plots_return_fig():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = SOC.load_harmonized(drop_na=True)
    flow = next(iter(SOC.focal_flows))
    for fn in (eng.plot_cross_level_box, eng.plot_within_std_box, eng.plot_distribution_overlay, eng.plot_biome_box):
        fig, ax = fn(SOC, df, flow)
        assert fig is not None and ax is not None
        plt.close(fig)
    fig, _ = eng.plot_sensitivity_heatmap(SOC, df)
    plt.close(fig)

    for fn in (xi.plot_significance_comparison, xi.plot_level_reframing, xi.plot_dispersion_by_level):
        fig, ax = fn()
        assert fig is not None
        plt.close(fig)
