"""Tests for the multi-indicator practice-change co-benefit analysis."""

import numpy as np

from sbtn_leaf.claude_analysis import practice_change as pc


# --------------------------------------------------------------------------- #
# Practice-flow resolution (light IO)
# --------------------------------------------------------------------------- #


def test_practice_flows_cereal_vs_tillage_only():
    wheat = pc.practice_flows("Wheat")
    assert wheat["kind"] == "cereal"
    assert wheat["baseline"] == "Wheat|rf|roff|ct"
    assert wheat["regen"] == "Wheat|rf|ron|rt"
    assert wheat["residue_only"] == "Wheat|rf|ron|ct"
    assert wheat["tillage_only"] == "Wheat|rf|roff|rt"

    soy = pc.practice_flows("Soybeans")
    assert soy["kind"] == "tillage_only"
    assert soy["baseline"] == "Soybeans|rf|na|ct"
    assert soy["regen"] == "Soybeans|rf|na|rt"
    assert soy["residue_only"] is None


# --------------------------------------------------------------------------- #
# Co-benefit table / summary
# --------------------------------------------------------------------------- #


def test_cobenefit_table_schema_and_content():
    t = pc.cobenefit_table("Wheat", "ecoregion")
    expected = {
        "commodity", "region_id", "soc_base", "soc_regen", "se_base", "se_regen",
        "d_soc", "d_soc_pct", "d_se_red", "d_se_red_pct", "win_win",
        "soc_rank", "se_rank", "priority", "biome", "realm",
    }
    assert expected.issubset(t.columns)
    assert len(t) > 100
    assert t["win_win"].dtype == bool
    assert np.isfinite(t["d_soc"]).all() and np.isfinite(t["d_se_red"]).all()
    # priority is a 0–1 composite of percentile ranks
    assert ((t["priority"] >= 0) & (t["priority"] <= 1)).all()
    # the switch reduces erosion everywhere (lower C-factor) -> non-negative
    assert (t["d_se_red"] >= -1e-9).all()


def test_cobenefit_summary_benefits_positive():
    s = pc.cobenefit_summary(["Wheat", "Soybeans"], "ecoregion")
    assert set(s["commodity"]) == {"Wheat", "Soybeans"}
    for _, row in s.iterrows():
        assert row["med_d_soc"] > 0          # carbon gained
        assert row["med_d_se_red"] > 0       # erosion avoided
        assert 0 <= row["pct_win_win"] <= 100
        assert -1 <= row["spearman_soc_se"] <= 1
    # the wheat switch is win-win for the large majority of regions
    assert s.loc[s["commodity"] == "Wheat", "pct_win_win"].iloc[0] > 80


def test_attribution_matches_manuscript_drivers():
    a = pc.practice_attribution_table(["Wheat"], "ecoregion").iloc[0]
    # residue retention drives the SOC gain; reduced tillage drives erosion control
    assert a["soc_residue"] > a["soc_tillage"]
    assert a["se_tillage"] > a["se_residue"]
    assert a["soc_driver"] == "residue" and a["se_driver"] == "tillage"


def test_priority_regions_sorted_and_bounded():
    p = pc.priority_regions("Wheat", "ecoregion", 10)
    assert len(p) == 10
    assert p["priority"].is_monotonic_decreasing
    assert ((p["priority"] >= 0) & (p["priority"] <= 1)).all()


# --------------------------------------------------------------------------- #
# Plot smoke tests (Agg backend, no display)
# --------------------------------------------------------------------------- #


def test_nonmap_plots_return_fig():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for fn, args in (
        (pc.plot_cobenefit_scatter, ("Wheat", "ecoregion")),
        (pc.plot_benefit_distributions, (["Wheat", "Maize"], "ecoregion")),
        (pc.plot_attribution, ("Wheat", "ecoregion")),
        (pc.plot_cobenefit_by_realm, ("Wheat", "ecoregion")),
    ):
        fig, _ = fn(*args)
        assert fig is not None
        plt.close(fig)


def test_map_degrades_or_renders_offline():
    # Offline-safe: with downloads disabled the map returns None when no geometry is
    # cached/local, or a (fig, ax) pair when it is. Either is acceptable.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res = pc.plot_cobenefit_map("Wheat", "ecoregion", "priority", allow_download=False)
    if res is not None:
        fig, ax = res
        assert fig is not None
        plt.close(fig)
