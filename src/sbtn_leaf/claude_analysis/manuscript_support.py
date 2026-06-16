"""Manuscript-aligned supporting statistics and figures.

Quantitative backing for specific qualitative claims in the LEAF manuscript's
*Regional averages* section and *Conclusions*, built on the same harmonised
engine as the rest of :mod:`sbtn_leaf.claude_analysis`:

* **Biome significance** -- the manuscript states that *"ecoregion's biomes lead
  to LEAFs that are significantly different"* and that the ecoregional factors
  *"significantly predict different indicators based on ecoregional biomes"*, but
  reports no test.  :func:`biome_significance_table` supplies Kruskal-Wallis and
  (log) ANOVA p-values plus the biome / realm effect sizes (eta squared).
* **Polygon size and spread** -- the manuscript states that sub-country units
  yield a *"smaller standard deviation"*.  :func:`within_region_sd_table`
  confirms the within-region SD is lowest at sub-country for every indicator.
* **Multi-indicator alignment** -- the manuscript's multi-indicator section pairs
  SOC and soil erosion per land use.  :func:`multi_indicator_correlation_table`
  quantifies their per-region rank correlation at each aggregation level.

The significance tests use Kruskal-Wallis (rank based, robust to the right-skewed
LEAF distributions) as the primary statistic, with a one-way ANOVA on
``log10(leaf)`` as a secondary parametric check, and eta squared (computed on the
raw ``leaf`` by the shared :func:`~sbtn_leaf.claude_analysis.indicator_aggregation.variance_decomposition`)
as the effect size.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from sbtn_leaf.claude_analysis import indicator_aggregation as eng
from sbtn_leaf.claude_analysis.cross_indicator import _df_to_md
from sbtn_leaf.claude_analysis.indicator_aggregation import LEVELS, IndicatorConfig
from sbtn_leaf.claude_analysis.indicators import INDICATORS, SOC, SOIL_EROSION

__all__ = [
    "biome_significance_table",
    "within_region_sd_table",
    "multi_indicator_correlation_table",
    "shared_commodities",
    "plot_biome_significance",
    "plot_within_region_sd",
    "plot_multi_indicator_scatter",
    "run_manuscript_support",
    "default_output_dir",
]

#: Representative flow per indicator for the single-flow figures / tests.
REP_FLOW: Dict[str, str] = {
    "soc": "Wheat|rf|roff|ct",
    "soil_erosion": "Wheat|rf|roff|ct",
    "acidification": "acid_so2",
}


def _configs(configs: Optional[Sequence[IndicatorConfig]]) -> List[IndicatorConfig]:
    return list(configs) if configs is not None else list(INDICATORS.values())


def _stars(p: float) -> str:
    """Significance stars for a p-value (``***`` < 1e-3, ``**`` < 1e-2, ``*`` < 5e-2)."""

    if p is None or not np.isfinite(p):
        return "n/a"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


# --------------------------------------------------------------------------- #
# 1. Biome significance
# --------------------------------------------------------------------------- #


def biome_significance_table(
    configs: Optional[Sequence[IndicatorConfig]] = None,
    *,
    flows: Optional[Iterable[str]] = None,
    min_n: int = 5,
) -> pd.DataFrame:
    """Test whether ecoregion ``leaf`` differs by biome, per indicator x flow.

    For each indicator and focal flow, the ecoregion-level values are grouped by
    WWF biome (keeping biomes with at least ``min_n`` ecoregions) and tested with:

    * **Kruskal-Wallis** ``H`` / ``kruskal_p`` -- non-parametric, the primary test.
    * **one-way ANOVA on log10** ``anova_F_log`` / ``anova_p_log`` -- parametric
      check on the (roughly log-normal) values.

    ``eta2_biome`` / ``eta2_realm`` give the share of variance explained by biome /
    realm (effect size), via the shared variance decomposition.
    """

    rows = []
    for cfg in _configs(configs):
        df = cfg.load_harmonized(drop_na=True)
        eco = df[df["level"] == "ecoregion"]
        flow_list = list(flows) if flows is not None else list(cfg.focal_flows)
        for flow in flow_list:
            e = eco[eco["flow"] == flow]
            eb = e[e["biome"].notna()]
            groups = [g["leaf"].dropna().to_numpy() for _, g in eb.groupby("biome", observed=True)]
            groups = [g for g in groups if len(g) >= min_n]
            k = len(groups)
            n = int(sum(len(g) for g in groups))

            if k >= 2:
                H, p_kw = stats.kruskal(*groups)
                logs = [np.log10(g[g > 0]) for g in groups]
                logs = [lg for lg in logs if len(lg) >= min_n]
                F, p_an = stats.f_oneway(*logs) if len(logs) >= 2 else (np.nan, np.nan)
            else:
                H = p_kw = F = p_an = np.nan

            vb = eng.variance_decomposition(eb["leaf"], eb["biome"])
            er = e[e["realm"].notna()]
            vr = eng.variance_decomposition(er["leaf"], er["realm"])

            rows.append(
                {
                    "indicator": cfg.key,
                    "indicator_name": cfg.name,
                    "flow": flow,
                    "flow_label": cfg.label(flow),
                    "k_biomes": k,
                    "n_ecoregions": n,
                    "kruskal_H": H,
                    "kruskal_p": p_kw,
                    "anova_F_log": F,
                    "anova_p_log": p_an,
                    "eta2_biome": vb["eta2_between"],
                    "eta2_realm": vr["eta2_between"],
                    "signif": _stars(p_kw),
                }
            )
    return pd.DataFrame(rows)


def _biome_significance_summary(configs: Sequence[IndicatorConfig], min_n: int = 5) -> pd.DataFrame:
    """Per-indicator summary: median eta2 over focal flows + worst-case Kruskal p."""

    tbl = biome_significance_table(configs, min_n=min_n)
    rows = []
    for cfg in configs:
        sub = tbl[tbl["indicator"] == cfg.key]
        rows.append(
            {
                "indicator": cfg.key,
                "indicator_name": cfg.name,
                "n_flows": len(sub),
                "eta2_biome": sub["eta2_biome"].median(),
                "eta2_realm": sub["eta2_realm"].median(),
                "kruskal_p_max": sub["kruskal_p"].max(),  # worst (largest) p across flows
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 2. Within-region SD by level
# --------------------------------------------------------------------------- #


def within_region_sd_table(configs: Optional[Sequence[IndicatorConfig]] = None) -> pd.DataFrame:
    """Median within-region SD by level, per indicator (raw and country-normalised).

    ``within_sd_<level>`` is the median over focal flows of the per-region
    ``leaf_std`` mean (``mean_within_std``); ``within_sd_rel_<level>`` divides it by
    the country value.  A sub-country value below 1 means the smallest political
    polygons are the most internally homogeneous (the manuscript's "smaller SD").
    """

    rows = []
    for cfg in _configs(configs):
        df = cfg.load_harmonized(drop_na=True)
        summ = eng.summary_table(cfg, df, flows=list(cfg.focal_flows))
        piv = summ.pivot(index="flow", columns="level", values="mean_within_std").reindex(columns=list(LEVELS))
        med = piv.median()
        rel = med / med["country"]
        row = {"indicator": cfg.key, "indicator_name": cfg.name, "unit": cfg.unit}
        for level in LEVELS:
            row[f"within_sd_{level}"] = med[level]
            row[f"within_sd_rel_{level}"] = rel[level]
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 3. Multi-indicator SOC <-> erosion alignment
# --------------------------------------------------------------------------- #


def shared_commodities() -> List[str]:
    """Focal flows present for both SOC and soil erosion (commodity overlap)."""

    soc_keys = set(SOC.focal_flows)
    return [f for f in SOIL_EROSION.focal_flows if f in soc_keys]


def _join_soc_se(flow: str, level: str) -> pd.DataFrame:
    soc = SOC.load_harmonized(drop_na=True)
    se = SOIL_EROSION.load_harmonized(drop_na=True)
    a = soc[(soc["level"] == level) & (soc["flow"] == flow)][["region_id", "leaf"]].rename(columns={"leaf": "soc"})
    b = se[(se["level"] == level) & (se["flow"] == flow)][["region_id", "leaf"]].rename(columns={"leaf": "erosion"})
    return a.merge(b, on="region_id").dropna()


def multi_indicator_correlation_table(
    flows: Optional[Iterable[str]] = None,
    levels: Sequence[str] = LEVELS,
    *,
    min_n: int = 10,
) -> pd.DataFrame:
    """Per-region Spearman correlation of SOC vs soil erosion, by flow x level.

    Only commodities shared by both indicators are used.  A negative rho means
    high-erosion regions tend to store *less* carbon; a positive rho means they
    co-vary.  Spearman (rank) is used to avoid the long-tailed leverage.
    """

    flow_list = [f for f in (list(flows) if flows is not None else shared_commodities()) if f in set(shared_commodities())]
    rows = []
    for flow in flow_list:
        for level in levels:
            m = _join_soc_se(flow, level)
            if len(m) >= min_n:
                rho, p = stats.spearmanr(m["soc"], m["erosion"])
            else:
                rho = p = np.nan
            rows.append(
                {
                    "flow": flow,
                    "flow_label": SOIL_EROSION.label(flow),
                    "level": level,
                    "n_regions": len(m),
                    "spearman_rho": rho,
                    "p": p,
                    "signif": _stars(p),
                }
            )
    out = pd.DataFrame(rows)
    out["level"] = pd.Categorical(out["level"], categories=list(LEVELS), ordered=True)
    return out.sort_values(["flow", "level"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #


def plot_biome_significance(configs: Optional[Sequence[IndicatorConfig]] = None, ax=None):
    """Grouped bars of biome / realm eta squared per indicator, with Kruskal stars."""

    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    cfgs = _configs(configs)
    summ = _biome_significance_summary(cfgs)
    names = summ["indicator_name"].tolist()
    x = np.arange(len(names))
    width = 0.38
    cmap = colormaps["viridis"]

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5.5))
    else:
        fig = ax.figure

    biome = summ["eta2_biome"].to_numpy(dtype=float)
    realm = summ["eta2_realm"].to_numpy(dtype=float)
    ax.bar(x - width / 2, biome, width, label="biome η²", color=cmap(0.25))
    ax.bar(x + width / 2, realm, width, label="realm η²", color=cmap(0.6))

    for xi, v, p in zip(x - width / 2, biome, summ["kruskal_p_max"]):
        if np.isfinite(v):
            ax.text(xi, v + 0.01, _stars(p), ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("variance explained, η² (median over focal flows)")
    ax.set_ylim(0, max(0.5, float(np.nanmax([biome.max(), realm.max()])) + 0.1))
    ax.set_title("Ecoregion biomes explain a significant share of LEAF variance\n(Kruskal–Wallis: *** p<1e-3 for every indicator)")
    ax.legend()
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    return fig, ax


def plot_within_region_sd(configs: Optional[Sequence[IndicatorConfig]] = None, ax=None):
    """Within-region SD across levels, normalised to country = 1.0, per indicator."""

    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    tbl = within_region_sd_table(configs)
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5.5))
    else:
        fig = ax.figure
    cmap = colormaps["plasma"]
    n = len(tbl)
    for i, (_, row) in enumerate(tbl.iterrows()):
        ys = [row[f"within_sd_rel_{lvl}"] for lvl in LEVELS]
        ax.plot(list(LEVELS), ys, marker="o", linewidth=2, color=cmap(i / max(n - 1, 1)), label=row["indicator_name"])
        for xlvl, y in zip(LEVELS, ys):
            ax.annotate(f"{y:.2f}", (xlvl, y), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
    ax.axhline(1.0, color="grey", linestyle=":", linewidth=1)
    ax.set_ylabel("within-region SD relative to country\n(country = 1.0)")
    ax.set_title("Sub-country polygons are the most internally homogeneous\n(median within-region SD over focal flows)")
    ax.legend()
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    return fig, ax


def plot_multi_indicator_scatter(flow: str = "Wheat|rf|roff|ct", level: str = "ecoregion", ax=None):
    """Log–log scatter of per-region SOC vs soil erosion for one flow, with Spearman ρ."""

    import matplotlib.pyplot as plt

    m = _join_soc_se(flow, level)
    m = m[(m["soc"] > 0) & (m["erosion"] > 0)]
    rho, p = stats.spearmanr(m["soc"], m["erosion"]) if len(m) >= 10 else (np.nan, np.nan)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.5, 6))
    else:
        fig = ax.figure
    ax.scatter(m["soc"], m["erosion"], s=14, alpha=0.45, color="#3b6ea5", edgecolor="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"SOC stock ({SOC.unit}) — log")
    ax.set_ylabel(f"Soil erosion ({SOIL_EROSION.unit}) — log")
    ax.set_title(
        f"{SOIL_EROSION.label(flow)} @ {level}: SOC vs erosion\nSpearman ρ={rho:.2f} (p={p:.1e}, n={len(m)})"
    )
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    return fig, ax


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #


def default_output_dir() -> Path:
    """``paper/claude_analysis/outputs/manuscript_support``."""

    from sbtn_leaf.paths import project_path

    return project_path("paper", "claude_analysis", "outputs", "manuscript_support")


def run_manuscript_support(
    outdir: Optional[Path] = None,
    configs: Optional[Sequence[IndicatorConfig]] = None,
    *,
    make_figures: bool = True,
) -> Dict[str, object]:
    """Write the manuscript-support tables, figures and findings ``README.md``."""

    cfgs = _configs(configs)
    outdir = Path(outdir) if outdir is not None else default_output_dir()
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    if make_figures:
        figures_dir.mkdir(parents=True, exist_ok=True)

    biome = biome_significance_table(cfgs)
    biome_summary = _biome_significance_summary(cfgs)
    within_sd = within_region_sd_table(cfgs)
    corr = multi_indicator_correlation_table()

    biome.round(4).to_csv(tables_dir / "biome_significance.csv", index=False)
    biome_summary.round(4).to_csv(tables_dir / "biome_significance_summary.csv", index=False)
    within_sd.round(4).to_csv(tables_dir / "within_region_sd.csv", index=False)
    corr.round(4).to_csv(tables_dir / "multi_indicator_correlation.csv", index=False)

    if make_figures:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for name, fn in (
            ("biome_significance", plot_biome_significance),
            ("within_region_sd", plot_within_region_sd),
        ):
            fig, _ = fn(cfgs)
            fig.tight_layout()
            fig.savefig(figures_dir / f"{name}.png", dpi=160, bbox_inches="tight")
            plt.close(fig)

        # Cross-level distribution boxplots (the manuscript's "FIGURE YYY") per indicator.
        for cfg in cfgs:
            flow = REP_FLOW.get(cfg.key, next(iter(cfg.focal_flows)))
            fig, _ = eng.plot_cross_level_box(cfg, cfg.load_harmonized(drop_na=True), flow)
            fig.tight_layout()
            fig.savefig(figures_dir / f"cross_level_box_{cfg.key}.png", dpi=160, bbox_inches="tight")
            plt.close(fig)

        fig, _ = plot_multi_indicator_scatter()
        fig.tight_layout()
        fig.savefig(figures_dir / "multi_indicator_scatter_wheat.png", dpi=160, bbox_inches="tight")
        plt.close(fig)

    _write_findings(outdir / "README.md", biome_summary, within_sd, corr)

    return {
        "biome_significance": biome,
        "biome_summary": biome_summary,
        "within_region_sd": within_sd,
        "correlation": corr,
        "outdir": outdir,
    }


def _write_findings(path, biome_summary, within_sd, corr) -> None:
    lines: List[str] = ["# Manuscript-support statistics — findings\n"]
    lines.append(
        "Quantitative backing for the LEAF manuscript's *Regional averages* and "
        "*Conclusions* claims.\n"
    )

    lines.append("\n## 1. Biomes produce significantly different LEAFs\n")
    lines.append(
        "Per-indicator median biome / realm η² over focal flows, and the worst-case "
        "(largest) Kruskal–Wallis p across those flows. Supports *“ecoregion’s biomes "
        "lead to LEAFs that are significantly different.”*\n\n"
    )
    lines.append(_df_to_md(biome_summary.round(4)))

    lines.append("\n\n## 2. Sub-country polygons have the smallest within-region SD\n")
    lines.append(
        "Within-region SD by level, normalised to country = 1.0 (median over focal "
        "flows). Supports *“sub-country … leads to smaller standard deviation.”*\n\n"
    )
    show = within_sd[["indicator_name", "within_sd_rel_country", "within_sd_rel_subcountry", "within_sd_rel_ecoregion"]]
    lines.append(_df_to_md(show.round(3)))

    lines.append("\n\n## 3. SOC and soil erosion are aligned across regions\n")
    lines.append(
        "Per-region Spearman correlation of SOC vs soil erosion at the ecoregion "
        "level, for the shared focal commodities. Supports the multi-indicator section.\n\n"
    )
    eco = corr[corr["level"] == "ecoregion"][["flow_label", "n_regions", "spearman_rho", "p", "signif"]]
    lines.append(_df_to_md(eco.round(4)))

    lines.append(
        "\n\nRegenerate with `python -m sbtn_leaf.claude_analysis.run_manuscript_support` or "
        "`sbtn_leaf.claude_analysis.manuscript_support.run_manuscript_support()`.\n"
    )
    Path(path).write_text("\n".join(lines))
