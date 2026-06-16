"""Tests for the manuscript-support statistics & figures."""

import numpy as np

from sbtn_leaf.claude_analysis import indicator_aggregation as eng
from sbtn_leaf.claude_analysis import manuscript_support as ms
from sbtn_leaf.claude_analysis.indicators import INDICATORS


# --------------------------------------------------------------------------- #
# 1. Biome significance
# --------------------------------------------------------------------------- #


def test_biome_significance_table_columns_and_bounds():
    tbl = ms.biome_significance_table()
    expected = {
        "indicator",
        "flow",
        "flow_label",
        "k_biomes",
        "n_ecoregions",
        "kruskal_H",
        "kruskal_p",
        "anova_F_log",
        "anova_p_log",
        "eta2_biome",
        "eta2_realm",
        "signif",
    }
    assert expected.issubset(tbl.columns)
    for col in ("eta2_biome", "eta2_realm"):
        vals = tbl[col].dropna()
        assert ((vals >= -1e-9) & (vals <= 1 + 1e-9)).all()
    # at least one row per indicator
    assert set(tbl["indicator"].unique()) == set(INDICATORS)


def test_biome_significance_representative_flow_is_significant():
    # The manuscript claim: biomes differ significantly. The representative flow of
    # every indicator should reject equal biome distributions overwhelmingly.
    tbl = ms.biome_significance_table()
    for key, flow in ms.REP_FLOW.items():
        row = tbl[(tbl["indicator"] == key) & (tbl["flow"] == flow)]
        assert len(row) == 1
        assert float(row["kruskal_p"].iloc[0]) < 1e-3
        assert row["signif"].iloc[0] == "***"


# --------------------------------------------------------------------------- #
# 2. Within-region SD by level
# --------------------------------------------------------------------------- #


def test_within_region_sd_subcountry_is_lowest():
    tbl = ms.within_region_sd_table()
    assert len(tbl) == len(INDICATORS)
    # country is the normalisation baseline
    assert np.allclose(tbl["within_sd_rel_country"], 1.0)
    # sub-country is the most internally homogeneous level for every indicator
    assert (tbl["within_sd_rel_subcountry"] < 1.0).all()
    assert (tbl["within_sd_rel_subcountry"] < tbl["within_sd_rel_ecoregion"]).all()


# --------------------------------------------------------------------------- #
# 3. Multi-indicator SOC <-> erosion correlation
# --------------------------------------------------------------------------- #


def test_multi_indicator_correlation_shape_and_bounds():
    corr = ms.multi_indicator_correlation_table()
    shared = ms.shared_commodities()
    assert set(corr["flow"].unique()) == set(shared)
    assert set(corr["level"].unique()) == set(eng.LEVELS)
    rho = corr["spearman_rho"].dropna()
    assert ((rho >= -1) & (rho <= 1)).all()
    assert (corr["n_regions"] >= 0).all()


def test_multi_indicator_ecoregion_correlation_significant_for_wheat():
    corr = ms.multi_indicator_correlation_table(flows=["Wheat|rf|roff|ct"], levels=["ecoregion"])
    row = corr.iloc[0]
    assert row["n_regions"] > 100
    assert np.isfinite(row["spearman_rho"])
    assert row["p"] < 1e-3


# --------------------------------------------------------------------------- #
# Plot smoke tests (Agg backend, no display)
# --------------------------------------------------------------------------- #


def test_plots_return_fig():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for fn in (ms.plot_biome_significance, ms.plot_within_region_sd):
        fig, ax = fn()
        assert fig is not None and ax is not None
        plt.close(fig)
    fig, ax = ms.plot_multi_indicator_scatter("Wheat|rf|roff|ct", "ecoregion")
    assert fig is not None
    plt.close(fig)
