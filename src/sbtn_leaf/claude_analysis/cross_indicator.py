"""Cross-indicator comparison: how ecoregions reframe LEAF averages.

This is the analytical heart of the question *"do ecoregions represent SOC, soil
erosion and terrestrial acidification averages differently than sub-country or
country units?"*.  It runs the generic engine over all three indicators and then
contrasts them on two axes, mirroring sections 6 and 7 of the soil-erosion
aggregation notebook but *across* indicators:

* **Ecoregion significance** (section-6 analogue) -- for each indicator, how much
  of the ecoregion-level variance an ecological grouping (biome / realm) explains
  (``eta2``), and how much within-country variance a single national LEAF hides
  (``frac_within``).  A high biome ``eta2`` means ecoregions carry signal that
  political polygons cannot.
* **Level reframing** (section-7 analogue) -- how the central estimate, the tail
  and the between-region spread move as polygons shrink country -> subcountry ->
  ecoregion.  Because the three indicators use different units, the central
  estimate is normalised to the country mean (country = 1.0) so they share one
  axis.

All metrics are summarised across each indicator's focal flows (median, which is
robust to the long-tailed commodity flows) so a single number per indicator is
comparable side by side.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from sbtn_leaf.claude_analysis.indicator_aggregation import (
    LEVELS,
    IndicatorConfig,
    coverage_table,
    hierarchy_variance_table,
    plot_biome_box,
    plot_sensitivity_heatmap,
    summary_table,
)
from sbtn_leaf.claude_analysis.indicators import INDICATORS

__all__ = [
    "combined_summary_table",
    "significance_table",
    "level_reframing_table",
    "dispersion_table",
    "plot_significance_comparison",
    "plot_level_reframing",
    "plot_dispersion_by_level",
    "run_cross_indicator_analysis",
    "default_output_dir",
]


# --------------------------------------------------------------------------- #
# Cross-indicator tables
# --------------------------------------------------------------------------- #


def _configs(configs: Optional[Sequence[IndicatorConfig]]) -> List[IndicatorConfig]:
    return list(configs) if configs is not None else list(INDICATORS.values())


def combined_summary_table(configs: Optional[Sequence[IndicatorConfig]] = None) -> pd.DataFrame:
    """Stack every indicator's focal ``flow x level`` summary into one long table."""

    frames = []
    for cfg in _configs(configs):
        df = cfg.load_harmonized(drop_na=True)
        summ = summary_table(cfg, df, flows=list(cfg.focal_flows))
        summ.insert(0, "indicator_name", cfg.name)
        summ.insert(0, "unit", cfg.unit)
        frames.append(summ)
    return pd.concat(frames, ignore_index=True)


def significance_table(
    configs: Optional[Sequence[IndicatorConfig]] = None, *, agg: str = "median"
) -> pd.DataFrame:
    """One row per indicator: ecoregion significance summarised over focal flows.

    ``eco_biome_eta2`` / ``eco_realm_eta2`` -- share of ecoregion-level variance
    explained by biome / realm (ecoregion significance, 0-1).
    ``subcty_within_country_frac`` -- share of admin-1 variance a single national
    LEAF hides inside its borders (political coarsening cost, 0-1).
    """

    rows = []
    for cfg in _configs(configs):
        df = cfg.load_harmonized(drop_na=True)
        var = hierarchy_variance_table(cfg, df, flows=list(cfg.focal_flows))
        agg_fn = getattr(var[["eco_between_biome_eta2", "eco_between_realm_eta2", "subcty_within_country_frac"]], agg)
        summary = agg_fn()
        rows.append(
            {
                "indicator": cfg.key,
                "indicator_name": cfg.name,
                "unit": cfg.unit,
                "n_focal_flows": var["flow"].nunique(),
                "eco_biome_eta2": summary["eco_between_biome_eta2"],
                "eco_realm_eta2": summary["eco_between_realm_eta2"],
                "subcty_within_country_frac": summary["subcty_within_country_frac"],
            }
        )
    return pd.DataFrame(rows)


def level_reframing_table(configs: Optional[Sequence[IndicatorConfig]] = None) -> pd.DataFrame:
    """How the central estimate moves across levels, per indicator.

    For each focal flow the between-region ``mean_leaf`` at each level is divided
    by its *country* mean (country = 1.0), then the ratios are aggregated (median)
    across focal flows.  Values > 1 mean finer polygons surface a *higher* average
    than the national figure (hotspots that countries average away); < 1 means
    coarsening inflates the estimate.
    """

    rows = []
    for cfg in _configs(configs):
        df = cfg.load_harmonized(drop_na=True)
        summ = summary_table(cfg, df, flows=list(cfg.focal_flows))
        wide = summ.pivot(index="flow", columns="level", values="mean_leaf").reindex(columns=list(LEVELS))
        med = summ.pivot(index="flow", columns="level", values="median_leaf").reindex(columns=list(LEVELS))
        # Normalise each flow to its own country mean / median, then median across flows.
        rel_mean = wide.div(wide["country"], axis=0)
        rel_median = med.div(med["country"], axis=0)
        row = {"indicator": cfg.key, "indicator_name": cfg.name, "unit": cfg.unit}
        for level in LEVELS:
            row[f"rel_mean_{level}"] = rel_mean[level].median()
            row[f"rel_median_{level}"] = rel_median[level].median()
        rows.append(row)
    return pd.DataFrame(rows)


def dispersion_table(configs: Optional[Sequence[IndicatorConfig]] = None) -> pd.DataFrame:
    """Between-region coefficient of variation by level, per indicator.

    Median (across focal flows) of the between-region CV at each level.  Rising CV
    from country to ecoregion means finer polygons *expose* more cross-region
    contrast; falling CV means they smooth it.
    """

    rows = []
    for cfg in _configs(configs):
        df = cfg.load_harmonized(drop_na=True)
        summ = summary_table(cfg, df, flows=list(cfg.focal_flows))
        cv = summ.pivot(index="flow", columns="level", values="cv").reindex(columns=list(LEVELS))
        row = {"indicator": cfg.key, "indicator_name": cfg.name}
        for level in LEVELS:
            row[f"cv_{level}"] = cv[level].median()
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Cross-indicator figures
# --------------------------------------------------------------------------- #


def plot_significance_comparison(configs: Optional[Sequence[IndicatorConfig]] = None, ax=None):
    """Grouped bars of ecoregion significance metrics across indicators."""

    import matplotlib.pyplot as plt

    tbl = significance_table(configs)
    metrics = [
        ("eco_biome_eta2", "biome η² (ecoregion)"),
        ("eco_realm_eta2", "realm η² (ecoregion)"),
        ("subcty_within_country_frac", "within-country var. hidden"),
    ]
    names = tbl["indicator_name"].tolist()
    x = np.arange(len(names))
    width = 0.25

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5.5))
    else:
        fig = ax.figure
    from matplotlib import colormaps

    cmap = colormaps["viridis"]
    for i, (col, label) in enumerate(metrics):
        ax.bar(x + (i - 1) * width, tbl[col].to_numpy(dtype=float), width, label=label, color=cmap(0.15 + 0.35 * i))
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("variance fraction (0–1)")
    ax.set_ylim(0, 1)
    ax.set_title("Ecoregion significance by indicator\n(median across focal flows)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    for i, (col, _) in enumerate(metrics):
        for xi, v in zip(x + (i - 1) * width, tbl[col].to_numpy(dtype=float)):
            if np.isfinite(v):
                ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    return fig, ax


def plot_level_reframing(configs: Optional[Sequence[IndicatorConfig]] = None, ax=None):
    """Line chart of the normalised central estimate (country = 1) across levels."""

    import matplotlib.pyplot as plt

    tbl = level_reframing_table(configs)
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5.5))
    else:
        fig = ax.figure
    from matplotlib import colormaps

    cmap = colormaps["plasma"]
    n = len(tbl)
    for i, (_, row) in enumerate(tbl.iterrows()):
        ys = [row[f"rel_mean_{lvl}"] for lvl in LEVELS]
        ax.plot(list(LEVELS), ys, marker="o", linewidth=2, color=cmap(i / max(n - 1, 1)), label=row["indicator_name"])
        for x, y in zip(LEVELS, ys):
            ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
    ax.axhline(1.0, color="grey", linestyle=":", linewidth=1)
    ax.set_ylabel("mean relative to country mean\n(country = 1.0)")
    ax.set_title("How the central estimate is reframed by polygon size\n(median across focal flows)")
    ax.legend()
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    return fig, ax


def plot_dispersion_by_level(configs: Optional[Sequence[IndicatorConfig]] = None, ax=None):
    """Grouped bars of between-region CV at each level, per indicator."""

    import matplotlib.pyplot as plt

    tbl = dispersion_table(configs)
    names = tbl["indicator_name"].tolist()
    x = np.arange(len(names))
    width = 0.25
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5.5))
    else:
        fig = ax.figure
    from matplotlib import colormaps

    cmap = colormaps["viridis"]
    for i, level in enumerate(LEVELS):
        ax.bar(x + (i - 1) * width, tbl[f"cv_{level}"].to_numpy(dtype=float), width, label=level, color=cmap(0.15 + 0.35 * i))
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("between-region CV (median over focal flows)")
    ax.set_title("Between-region dispersion by aggregation level")
    ax.legend(title="level")
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    return fig, ax


# --------------------------------------------------------------------------- #
# One-shot pipeline
# --------------------------------------------------------------------------- #


def default_output_dir() -> Path:
    """``paper/claude_analysis/outputs`` -- where the deliverable artefacts live."""

    from sbtn_leaf.paths import project_path

    return project_path("paper", "claude_analysis", "outputs")


def run_cross_indicator_analysis(
    outdir: Optional[Path] = None,
    configs: Optional[Sequence[IndicatorConfig]] = None,
    *,
    make_figures: bool = True,
) -> Dict[str, object]:
    """Run the whole cross-indicator comparison, writing tables (+ figures + README).

    Writes, under ``outdir`` (default :func:`default_output_dir`):

    * ``tables/`` -- per-indicator coverage / focal summary / hierarchy variance,
      plus the cross-indicator significance / reframing / dispersion tables.
    * ``figures/`` -- per-indicator biome boxplot + sensitivity heatmap, and the
      three cross-indicator comparison figures.
    * ``README.md`` -- the findings narrative.

    Returns a dict of the key DataFrames for interactive use.
    """

    cfgs = _configs(configs)
    outdir = Path(outdir) if outdir is not None else default_output_dir()
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    if make_figures:
        figures_dir.mkdir(parents=True, exist_ok=True)

    # -- per-indicator tables -------------------------------------------------
    per_indicator: Dict[str, Dict[str, pd.DataFrame]] = {}
    for cfg in cfgs:
        raw = cfg.load_harmonized(drop_na=False)
        df = raw[raw["leaf"].notna()].reset_index(drop=True)
        cov = coverage_table(raw)
        summ = summary_table(cfg, df, flows=list(cfg.focal_flows))
        var = hierarchy_variance_table(cfg, df, flows=list(cfg.focal_flows))
        cov.to_csv(tables_dir / f"{cfg.key}_coverage.csv", index=False)
        summ.round(4).to_csv(tables_dir / f"{cfg.key}_focal_summary.csv", index=False)
        var.round(4).to_csv(tables_dir / f"{cfg.key}_hierarchy_variance.csv", index=False)
        per_indicator[cfg.key] = {"raw": raw, "df": df, "coverage": cov, "summary": summ, "variance": var}

    # -- cross-indicator tables ----------------------------------------------
    combined = combined_summary_table(cfgs)
    significance = significance_table(cfgs)
    reframing = level_reframing_table(cfgs)
    dispersion = dispersion_table(cfgs)
    combined.round(4).to_csv(tables_dir / "combined_focal_summary.csv", index=False)
    significance.round(4).to_csv(tables_dir / "ecoregion_significance.csv", index=False)
    reframing.round(4).to_csv(tables_dir / "level_reframing.csv", index=False)
    dispersion.round(4).to_csv(tables_dir / "dispersion_by_level.csv", index=False)

    # -- figures --------------------------------------------------------------
    if make_figures:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for cfg in cfgs:
            df = per_indicator[cfg.key]["df"]
            focal0 = next(iter(cfg.focal_flows))
            fig, _ = plot_biome_box(cfg, df, focal0)
            fig.tight_layout()
            fig.savefig(figures_dir / f"{cfg.key}_biome_box.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

            fig, _ = plot_sensitivity_heatmap(cfg, df)
            fig.tight_layout()
            fig.savefig(figures_dir / f"{cfg.key}_sensitivity_heatmap.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

        for name, fn in (
            ("significance_comparison", plot_significance_comparison),
            ("level_reframing", plot_level_reframing),
            ("dispersion_by_level", plot_dispersion_by_level),
        ):
            fig, _ = fn(cfgs)
            fig.tight_layout()
            fig.savefig(figures_dir / f"{name}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

    _write_findings(outdir / "README.md", cfgs, per_indicator, significance, reframing, dispersion)

    return {
        "per_indicator": per_indicator,
        "combined": combined,
        "significance": significance,
        "reframing": reframing,
        "dispersion": dispersion,
        "outdir": outdir,
    }


# --------------------------------------------------------------------------- #
# Markdown findings writer (no extra deps)
# --------------------------------------------------------------------------- #


def _fmt(v) -> str:
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def _df_to_md(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    head = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = ["| " + " | ".join(_fmt(v) for v in row) + " |" for row in df.itertuples(index=False, name=None)]
    return "\n".join([head, sep, *body])


def _write_findings(path, cfgs, per_indicator, significance, reframing, dispersion) -> None:
    lines: List[str] = ["# Cross-indicator ecoregion aggregation — findings\n"]
    lines.append(
        "How three SBTN-Land LEAFs — **SOC**, **soil erosion** and **terrestrial "
        "acidification** — behave as polygons shrink country → subcountry → "
        "ecoregion, and how much ecological grouping (biome / realm) explains.\n"
    )
    lines.append("Indicators and units:\n")
    lines.append(
        _df_to_md(
            pd.DataFrame(
                [
                    {"indicator": c.name, "unit": c.unit, "flows": c.flow_kind, "direction": c.direction_word}
                    for c in cfgs
                ]
            )
        )
    )

    lines.append("\n\n## Coverage by level (share of region × flow cells with a value)\n")
    cov = pd.concat(
        [per_indicator[c.key]["coverage"].assign(indicator=c.name) for c in cfgs], ignore_index=True
    )
    lines.append(_df_to_md(cov[["indicator", "level", "n_regions", "n_flows", "coverage_pct"]].round(1)))

    lines.append("\n\n## Ecoregion significance (median across focal flows)\n")
    lines.append(
        "`eco_biome_eta2` / `eco_realm_eta2` = share of ecoregion-level variance "
        "explained by biome / realm (higher ⇒ ecoregions carry signal political "
        "units miss). `subcty_within_country_frac` = share of admin-1 variance a "
        "single national LEAF hides.\n\n"
    )
    lines.append(_df_to_md(significance[["indicator_name", "n_focal_flows", "eco_biome_eta2", "eco_realm_eta2", "subcty_within_country_frac"]].round(3)))

    lines.append("\n\n## Level reframing of the central estimate (country mean = 1.0)\n")
    lines.append(_df_to_md(reframing.round(3)))

    lines.append("\n\n## Between-region dispersion (median CV by level)\n")
    lines.append(_df_to_md(dispersion.round(3)))

    lines.append(
        "\n\nTables live in `tables/`; figures in `figures/`. Regenerate everything "
        "with `python -m sbtn_leaf.claude_analysis.run_cross_indicator` or "
        "`sbtn_leaf.claude_analysis.cross_indicator.run_cross_indicator_analysis()`.\n"
    )
    Path(path).write_text("\n".join(lines))
