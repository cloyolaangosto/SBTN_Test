"""Tests for the soil-erosion aggregation-comparison module."""

import math

import numpy as np
import pandas as pd
import pytest

from sbtn_leaf.claude_analysis import se_aggregation_analysis as se


# --------------------------------------------------------------------------- #
# Flow-name canonicalisation
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "raw, expected",
    [
        # short code (country / subcountry) vs descriptive (ecoregion)
        ("Wheat_rf_roff", "Wheat|rf|roff|ct"),
        ("se_rate_25km_clipped_Rainfed_Wheat_residues_removed_from_the_field", "Wheat|rf|roff|ct"),
        ("Maize_irr_ron", "Maize|irr|ron|ct"),
        ("se_rate_25km_clipped_Irrigated_Maize_residues_left_on_the_field", "Maize|irr|ron|ct"),
        # reduced tillage keeps short codes in both schemes
        ("Wheat_rf_roff_rt", "Wheat|rf|roff|rt"),
        ("se_rate_25km_clipped_Wheat_rf_roff_rt", "Wheat|rf|roff|rt"),
        # simple crops (no residue/tillage tokens)
        ("Soybeans_rf", "Soybeans|rf|na|ct"),
        ("se_rate_25km_clipped_Rainfed_Soybeans", "Soybeans|rf|na|ct"),
        ("Oil_palm_rf", "Oil_palm|rf|na|ct"),
        ("se_rate_25km_clipped_Rainfed_Oil_palm", "Oil_palm|rf|na|ct"),
        # lowercase cotton conventional vs reduced-till Cotton
        ("cotton_rf", "Cotton|rf|na|ct"),
        ("se_rate_25km_clipped_Rainfed_cotton", "Cotton|rf|na|ct"),
        ("Cotton_rf_rt", "Cotton|rf|na|rt"),
        # land cover
        ("BRDC_Tropical", "BRDC_Tropical"),
        ("se_rate_25km_clipped_Broadleaf_Deciduous_Tropical", "BRDC_Tropical"),
        ("NEEV_Boreal_dry", "NEEV_Boreal_dry"),
        ("se_rate_25km_clipped_Needleleaf_Evergreen_Boreal_dry", "NEEV_Boreal_dry"),
        ("Grassland", "Grassland"),
        ("Urban", "Urban"),
    ],
)
def test_canonical_flow(raw, expected):
    assert se.canonical_flow(raw) == expected


def test_flow_label_is_readable():
    assert se.flow_label("Wheat|rf|roff|ct") == "Wheat (rainfed)"
    # a non-focal crop falls back to the generated label
    assert "reduced till." in se.flow_label("Maize|irr|ron|rt")
    assert se.flow_label("NEEV_Tropical").startswith("Needleleaf-everg.")


# --------------------------------------------------------------------------- #
# Variance decomposition (pure maths, no IO)
# --------------------------------------------------------------------------- #


def test_variance_decomposition_known_case():
    res = se.variance_decomposition([1, 2, 3, 4], ["A", "A", "B", "B"])
    assert res["ss_total"] == pytest.approx(5.0)
    assert res["ss_between"] == pytest.approx(4.0)
    assert res["ss_within"] == pytest.approx(1.0)
    assert res["eta2_between"] == pytest.approx(0.8)
    assert res["frac_within"] == pytest.approx(0.2)


def test_variance_decomposition_degenerate():
    res = se.variance_decomposition([1.0, 2.0, 3.0], ["A", "A", "A"])  # single group
    assert res["n_groups"] == 1
    assert math.isnan(res["eta2_between"])


# --------------------------------------------------------------------------- #
# Loading & harmonisation (touch the real CSVs once, module scope)
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def raw():
    return se.load_harmonized(drop_na=False)


def test_load_level_country_schema():
    df = se.load_level("country", drop_na=True)
    assert list(df.columns) == se._STD_COLS
    assert (df["level"] == "country").all()
    assert df["leaf"].notna().all()
    assert df["country_name"].notna().all()


def test_load_level_ecoregion_pivot():
    df = se.load_level("ecoregion", drop_na=True)
    # long -> wide pivot must yield the three metric columns
    for col in ("leaf", "leaf_median", "leaf_std"):
        assert col in df.columns
    # ecological grouping fields are populated only here
    assert df["biome"].notna().any()
    assert df["realm"].notna().any()


def test_all_levels_present(raw):
    assert set(raw["level"].unique()) == set(se.LEVELS)


def test_flow_coverage_complete(raw):
    cov = se.validate_flow_coverage(raw)
    # the canonicaliser yields a complete 106-flow match across all three levels
    assert int(cov["in_all"].sum()) == 106
    for focal in se.FOCAL_FLOWS:
        assert cov.loc[focal, "in_all"], f"focal flow {focal} missing from a level"


def test_no_unparsed_flows(raw):
    assert not raw["flow"].astype(str).str.startswith("UNPARSED").any()


def test_coverage_table_orders_and_increases(raw):
    cov = se.coverage_table(raw)
    assert list(cov["level"]) == list(se.LEVELS)
    # finer levels resolve progressively more land (coverage increases)
    assert cov["coverage_pct"].is_monotonic_increasing


# --------------------------------------------------------------------------- #
# Summary statistics
# --------------------------------------------------------------------------- #


def test_summary_table_focal(raw):
    df = raw[raw["leaf"].notna()]
    summ = se.summary_table(df, flows=list(se.FOCAL_FLOWS))
    expected = {
        "flow",
        "level",
        "n_regions",
        "mean_leaf",
        "median_leaf",
        "std_leaf",
        "p10",
        "p90",
        "max_leaf",
        "mean_within_std",
        "cv",
        "flow_label",
    }
    assert expected.issubset(summ.columns)
    # each focal flow has a row per level
    counts = summ.groupby("flow", observed=True)["level"].nunique()
    assert (counts == len(se.LEVELS)).all()
    assert np.isfinite(summ["mean_leaf"]).all()


def test_hierarchy_variance_bounds(raw):
    df = raw[raw["leaf"].notna()]
    var = se.hierarchy_variance_table(df, flows=list(se.FOCAL_FLOWS))
    for col in ("subcty_within_country_frac", "eco_between_biome_eta2", "eco_between_realm_eta2"):
        vals = var[col].dropna()
        assert ((vals >= -1e-9) & (vals <= 1 + 1e-9)).all()


# --------------------------------------------------------------------------- #
# Plot smoke tests (Agg backend, no display)
# --------------------------------------------------------------------------- #


def test_plots_return_fig(raw):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = raw[raw["leaf"].notna()]
    flow = "Wheat|rf|roff|ct"
    for fn in (se.plot_cross_level_box, se.plot_within_std_box, se.plot_distribution_overlay, se.plot_biome_box):
        fig, ax = fn(df, flow)
        assert fig is not None and ax is not None
        plt.close(fig)
    fig, ax = se.plot_commodity_ranking(df)
    plt.close(fig)
    fig, ax = se.plot_sensitivity_heatmap(df)
    plt.close(fig)


def test_maps_skip_without_geometry(raw):
    # No boundary geometry is present in CI; the map helper must degrade to None.
    df = raw[raw["leaf"].notna()]
    with pytest.warns(UserWarning):
        result = se.plot_choropleth_panels(df, "Wheat|rf|roff|ct")
    assert result is None
