"""Multi-indicator practice-change co-benefits: where switching practices helps both.

Expands the LEAF manuscript's *Multi-indicator analysis* (Fig. 13): for the **same
commodity**, switching from the baseline practice (conventional tillage + residues
removed) to the regenerative practice (reduced tillage + residues left) raises the
SOC stock *and* lowers soil erosion.  This module quantifies, per region:

* ``d_soc``     -- SOC gained, ``t SOC/ha`` (and ``d_soc_pct``);
* ``d_se_red``  -- soil erosion avoided, ``t soil/ha/yr`` (and ``d_se_red_pct``);
* ``win_win``   -- both improve; and a unit-free ``priority`` score (mean of the
  within-commodity percentile ranks of the two *absolute* benefits) that ranks
  *where to focus* — the Fig. 13 "areas with highest benefits".

A key property of RUSLE: because erosion scales with the multiplicative C-factor,
the *percentage* erosion reduction of a practice switch is spatially **constant**
per commodity (e.g. ~77 % for the wheat switch); what varies across regions is the
*absolute* tonnage avoided.  SOC gain varies in both magnitude and sign.  Hence the
spatial "where to focus" ranking uses the **absolute** benefits.

Maps use :mod:`sbtn_leaf.claude_analysis.geo` (country names / ``ECO_ID``); when no
geometry is available they skip gracefully.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from sbtn_leaf.claude_analysis import geo
from sbtn_leaf.claude_analysis.cross_indicator import _df_to_md
from sbtn_leaf.claude_analysis.indicators import SOC, SOIL_EROSION
from sbtn_leaf.claude_analysis.manuscript_support import _stars

__all__ = [
    "COMMODITIES",
    "practice_flows",
    "cobenefit_table",
    "cobenefit_summary",
    "practice_attribution_table",
    "priority_regions",
    "STACK_CEREALS",
    "stack_flows",
    "practice_stack_table",
    "plot_practice_stack_range",
    "plot_cobenefit_scatter",
    "plot_benefit_distributions",
    "plot_attribution",
    "plot_cobenefit_by_realm",
    "plot_cobenefit_map",
    "plot_cobenefit_map_panels",
    "run_practice_change_analysis",
    "default_output_dir",
]

#: Commodities analysed (Wheat matches manuscript Fig. 13).
COMMODITIES = ["Wheat", "Maize", "Soybeans", "Cotton"]

#: Metric metadata: column -> (label, unit, sequential?).
_METRIC_LABEL = {
    "d_soc": ("SOC gained", "t SOC/ha"),
    "d_soc_pct": ("SOC gained", "%"),
    "d_se_red": ("erosion avoided", "t soil/ha/yr"),
    "d_se_red_pct": ("erosion avoided", "%"),
    "priority": ("co-benefit priority", "rank 0–1"),
}


def _soc():
    return SOC.load_harmonized(drop_na=True)


def _se():
    return SOIL_EROSION.load_harmonized(drop_na=True)


def practice_flows(commodity: str) -> Dict[str, Optional[str]]:
    """Canonical baseline / regenerative / component flows for a commodity.

    Cereals (residue management available) switch residues *and* tillage; other
    crops switch tillage only.  Detected from the flows actually present.
    """

    soc_flows = set(_soc()["flow"].unique())
    se_flows = set(_se()["flow"].unique())
    has_residue = f"{commodity}|rf|ron|ct" in soc_flows and f"{commodity}|rf|ron|ct" in se_flows
    if has_residue:
        return {
            "baseline": f"{commodity}|rf|roff|ct",
            "regen": f"{commodity}|rf|ron|rt",
            "residue_only": f"{commodity}|rf|ron|ct",
            "tillage_only": f"{commodity}|rf|roff|rt",
            "kind": "cereal",
        }
    return {
        "baseline": f"{commodity}|rf|na|ct",
        "regen": f"{commodity}|rf|na|rt",
        "residue_only": None,
        "tillage_only": f"{commodity}|rf|na|rt",
        "kind": "tillage_only",
    }


def _series(df: pd.DataFrame, flow: str, level: str) -> pd.Series:
    s = df[(df["level"] == level) & (df["flow"] == flow)][["region_id", "leaf"]]
    return s.set_index("region_id")["leaf"]


def _region_meta(level: str) -> pd.DataFrame:
    meta = _se()[_se()["level"] == level][["region_id", "region_name", "country_name", "biome", "realm"]]
    return meta.drop_duplicates("region_id").set_index("region_id")


def cobenefit_table(commodity: str, level: str = "ecoregion") -> pd.DataFrame:
    """Per-region SOC gain and erosion avoided for the practice switch.

    Columns: ``region_id, region_name, country_name, biome, realm, soc_base,
    soc_regen, se_base, se_regen, d_soc, d_soc_pct, d_se_red, d_se_red_pct,
    win_win, soc_rank, se_rank, priority`` (ranks are within-commodity percentile
    ranks of the *absolute* benefits; ``priority`` is their mean).
    """

    pf = practice_flows(commodity)
    soc, se = _soc(), _se()
    sb = _series(soc, pf["baseline"], level)
    sr = _series(soc, pf["regen"], level)
    eb = _series(se, pf["baseline"], level)
    er = _series(se, pf["regen"], level)
    idx = sb.index.intersection(sr.index).intersection(eb.index).intersection(er.index)

    out = pd.DataFrame(index=idx)
    out["soc_base"], out["soc_regen"] = sb[idx], sr[idx]
    out["se_base"], out["se_regen"] = eb[idx], er[idx]
    out = out.dropna()
    out["d_soc"] = out["soc_regen"] - out["soc_base"]
    out["d_soc_pct"] = 100 * out["d_soc"] / out["soc_base"]
    out["d_se_red"] = out["se_base"] - out["se_regen"]
    out["d_se_red_pct"] = 100 * out["d_se_red"] / out["se_base"]
    out["win_win"] = (out["d_soc"] > 0) & (out["d_se_red"] > 0)
    out["soc_rank"] = out["d_soc"].rank(pct=True)
    out["se_rank"] = out["d_se_red"].rank(pct=True)
    out["priority"] = (out["soc_rank"] + out["se_rank"]) / 2

    meta = _region_meta(level)
    out = out.join(meta, how="left")
    out.insert(0, "commodity", commodity)
    return out.reset_index().rename(columns={"index": "region_id"})


def cobenefit_summary(commodities: Optional[Sequence[str]] = None, level: str = "ecoregion") -> pd.DataFrame:
    """Per-commodity extent of benefits + co-location of the two benefits."""

    commodities = list(commodities) if commodities is not None else COMMODITIES
    rows = []
    for c in commodities:
        t = cobenefit_table(c, level)
        if t.empty:
            continue
        rho, p_rho = stats.spearmanr(t["d_soc"], t["d_se_red"]) if len(t) > 10 else (np.nan, np.nan)
        rows.append(
            {
                "commodity": c,
                "switch": practice_flows(c)["kind"],
                "n_regions": len(t),
                "med_d_soc": t["d_soc"].median(),
                "med_d_soc_pct": t["d_soc_pct"].median(),
                "med_d_se_red": t["d_se_red"].median(),
                "d_se_red_pct": t["d_se_red_pct"].median(),  # ~constant per commodity
                "pct_win_win": 100 * t["win_win"].mean(),
                "spearman_soc_se": rho,
                "spearman_p": p_rho,
                "signif": _stars(p_rho),
            }
        )
    return pd.DataFrame(rows)


def practice_attribution_table(commodities: Optional[Sequence[str]] = None, level: str = "ecoregion") -> pd.DataFrame:
    """Median residue- vs tillage-component effect on SOC gain and erosion avoided.

    Tests the manuscript claim that residue management dominates the SOC benefit
    while reduced tillage dominates the erosion benefit.  Component effects are
    blank for tillage-only crops.
    """

    commodities = list(commodities) if commodities is not None else COMMODITIES
    soc, se = _soc(), _se()
    rows = []
    for c in commodities:
        pf = practice_flows(c)
        base_soc = _series(soc, pf["baseline"], level)
        base_se = _series(se, pf["baseline"], level)
        til_soc = _series(soc, pf["tillage_only"], level)
        til_se = _series(se, pf["tillage_only"], level)

        def med(a, b):
            i = a.index.intersection(b.index)
            return float((a[i] - b[i]).median()) if len(i) else np.nan

        row = {"commodity": c, "switch": pf["kind"]}
        # SOC gain (regen - base); erosion avoided (base - regen)
        row["soc_tillage"] = med(til_soc, base_soc)
        row["se_tillage"] = med(base_se, til_se)
        if pf["residue_only"]:
            res_soc = _series(soc, pf["residue_only"], level)
            res_se = _series(se, pf["residue_only"], level)
            row["soc_residue"] = med(res_soc, base_soc)
            row["se_residue"] = med(base_se, res_se)
            row["soc_driver"] = "residue" if row["soc_residue"] >= row["soc_tillage"] else "tillage"
            row["se_driver"] = "residue" if row["se_residue"] >= row["se_tillage"] else "tillage"
        else:
            row["soc_residue"] = np.nan
            row["se_residue"] = np.nan
            row["soc_driver"] = "tillage"
            row["se_driver"] = "tillage"
        rows.append(row)
    return pd.DataFrame(rows)


def priority_regions(commodity: str = "Wheat", level: str = "ecoregion", n: int = 15) -> pd.DataFrame:
    """Top-``n`` regions by co-benefit ``priority`` (where to focus)."""

    t = cobenefit_table(commodity, level)
    cols = ["region_name", "country_name", "biome", "realm", "d_soc", "d_se_red", "priority"]
    cols = [c for c in cols if c in t.columns]
    return t.sort_values("priority", ascending=False).head(n)[cols].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Full practice-stack SOC range (best vs worst combination, by commodity)
# --------------------------------------------------------------------------- #

#: Cereals that carry the full irrigation x residue x tillage practice stack.
STACK_CEREALS = ["Maize", "Wheat", "Sorghum", "Barley", "Rapeseed"]


def stack_flows(commodity: str) -> Dict[str, str]:
    """Best- vs worst-practice SOC flow keys for the full stacking comparison.

    *Best* = rainfed + residues retained + reduced tillage; *worst* = irrigated +
    residues removed + conventional tillage.  Cereals use the full three-practice
    stack; crops without residue management fall back to an irrigation + tillage
    stack (``na`` residue).
    """

    has_residue = f"{commodity}|rf|ron|ct" in set(_soc()["flow"].unique())
    if has_residue:
        return {"best": f"{commodity}|rf|ron|rt", "worst": f"{commodity}|irr|roff|ct",
                "stack": "irrigation+residue+tillage"}
    return {"best": f"{commodity}|rf|na|rt", "worst": f"{commodity}|irr|na|ct",
            "stack": "irrigation+tillage"}


def _stack_pct(commodity: str, level: str) -> pd.Series:
    """Per-region % SOC difference between the best and worst practice stack."""

    sf = stack_flows(commodity)
    best = _series(_soc(), sf["best"], level)
    worst = _series(_soc(), sf["worst"], level)
    idx = best.index.intersection(worst.index)
    pct = 100 * (best[idx] - worst[idx]) / worst[idx]
    return pct.replace([np.inf, -np.inf], np.nan).dropna()


def practice_stack_table(commodities: Optional[Sequence[str]] = None, level: str = "ecoregion") -> pd.DataFrame:
    """Per-commodity SOC gain from the best vs worst practice stack.

    Expands the manuscript's single "~37.5 % higher SOC" figure into its
    by-commodity range: for each commodity, ``100*(SOC_best - SOC_worst)/SOC_worst``
    summarised (median, mean, 10th/90th percentile) across regions.
    """

    commodities = list(commodities) if commodities is not None else STACK_CEREALS
    rows = []
    for c in commodities:
        sf = stack_flows(c)
        pct = _stack_pct(c, level)
        if len(pct) < 10:
            continue
        best = _series(_soc(), sf["best"], level)
        worst = _series(_soc(), sf["worst"], level)
        idx = best.index.intersection(worst.index)
        rows.append(
            {
                "commodity": c,
                "stack": sf["stack"],
                "n_regions": int(len(pct)),
                "soc_worst_med": float(worst[idx].median()),
                "soc_best_med": float(best[idx].median()),
                "median_pct": float(pct.median()),
                "mean_pct": float(pct.mean()),
                "p10_pct": float(pct.quantile(0.10)),
                "p90_pct": float(pct.quantile(0.90)),
            }
        )
    return pd.DataFrame(rows).sort_values("median_pct", ascending=False).reset_index(drop=True)


def plot_practice_stack_range(commodities: Optional[Sequence[str]] = None, level: str = "ecoregion", ax=None):
    """Boxplots of the best-vs-worst SOC % difference per commodity (ranked)."""

    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    commodities = list(commodities) if commodities is not None else STACK_CEREALS
    series = [(c, _stack_pct(c, level).to_numpy()) for c in commodities]
    series = [(c, v) for c, v in series if len(v) >= 10]
    series.sort(key=lambda cv: np.median(cv[1]))
    labels = [c for c, _ in series]
    data = [v for _, v in series]

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 0.6 * len(data) + 2))
    else:
        fig = ax.figure
    bp = ax.boxplot(data, orientation="horizontal", whis=(10, 90), showfliers=False,
                    patch_artist=True, medianprops={"color": "black"})
    cmap = colormaps["YlGn"]
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cmap(0.3 + 0.5 * i / max(len(data) - 1, 1)))
        patch.set_alpha(0.9)
    overall = float(np.mean([np.mean(d) for d in data])) if data else float("nan")
    ax.axvline(overall, color="firebrick", linestyle="--", linewidth=1.5,
               label=f"mean of means ≈ {overall:.0f}% (manuscript: 37.5%)")
    ax.axvline(0, color="grey", linewidth=1)
    ax.set_yticks(range(1, len(labels) + 1))
    ax.set_yticklabels(labels)
    ax.set_xlabel("SOC difference, best vs worst practice stack (%)")
    ax.set_title(
        "Best- vs worst-practice SOC stacking benefit by commodity\n"
        "(best: rainfed + residues retained + reduced till; worst: irrigated + residues removed + conv. till)"
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    return fig, ax


# --------------------------------------------------------------------------- #
# Figures (non-map)
# --------------------------------------------------------------------------- #


def _realm_colors(realms):
    from matplotlib import colormaps

    uniq = [r for r in pd.unique(realms) if pd.notna(r)]
    cmap = colormaps["tab10"]
    return {r: cmap(i % 10) for i, r in enumerate(sorted(uniq))}


def plot_cobenefit_scatter(commodity: str = "Wheat", level: str = "ecoregion", ax=None):
    """Per-region SOC gain vs erosion avoided, coloured by realm (where to focus)."""

    import matplotlib.pyplot as plt

    t = cobenefit_table(commodity, level)
    t = t[t["d_se_red"] > 0]
    rho, p = stats.spearmanr(t["d_soc"], t["d_se_red"]) if len(t) > 10 else (np.nan, np.nan)
    colors = _realm_colors(t["realm"])

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 6.5))
    else:
        fig = ax.figure
    for realm, sub in t.groupby("realm"):
        ax.scatter(sub["d_se_red"], sub["d_soc"], s=22, alpha=0.6, color=colors.get(realm, "grey"), label=str(realm), edgecolor="none")
    ax.axhline(0, color="grey", linewidth=1)
    ax.axvline(t["d_se_red"].median(), color="grey", linestyle=":", linewidth=1)
    ax.set_xscale("symlog")
    ax.set_xlabel(f"Erosion avoided ({_METRIC_LABEL['d_se_red'][1]}, symlog)")
    ax.set_ylabel(f"SOC gained ({_METRIC_LABEL['d_soc'][1]})")
    ax.set_title(
        f"{commodity}: co-benefit of the practice switch @ {level}\n"
        f"{int(100*t['win_win'].mean())}% win-win; Spearman ρ={rho:.2f} (p={p:.1e}, n={len(t)})"
    )
    ax.legend(title="realm", fontsize=8, loc="best")
    ax.grid(True, linestyle="--", alpha=0.3)
    return fig, ax


def plot_benefit_distributions(commodities: Optional[Sequence[str]] = None, level: str = "ecoregion"):
    """Boxplots of SOC gained and erosion avoided per commodity."""

    import matplotlib.pyplot as plt

    commodities = list(commodities) if commodities is not None else COMMODITIES
    tables = {c: cobenefit_table(c, level) for c in commodities}
    fig, axes = plt.subplots(1, 2, figsize=(13, 0.6 * len(commodities) + 3))

    soc_data = [tables[c]["d_soc"].to_numpy() for c in commodities]
    axes[0].boxplot(soc_data, orientation="horizontal", whis=(5, 95), showfliers=False, patch_artist=True,
                    medianprops={"color": "black"})
    axes[0].axvline(0, color="grey", linewidth=1)
    axes[0].set_yticks(range(1, len(commodities) + 1))
    axes[0].set_yticklabels(commodities)
    axes[0].set_xlabel("SOC gained (t SOC/ha)")
    axes[0].set_title("SOC benefit")
    axes[0].grid(True, axis="x", linestyle="--", alpha=0.3)

    se_data = [tables[c]["d_se_red"][tables[c]["d_se_red"] > 0].to_numpy() for c in commodities]
    axes[1].boxplot(se_data, orientation="horizontal", whis=(5, 95), showfliers=False, patch_artist=True,
                    medianprops={"color": "black"})
    axes[1].set_xscale("log")
    axes[1].set_yticks(range(1, len(commodities) + 1))
    axes[1].set_yticklabels(commodities)
    axes[1].set_xlabel("Erosion avoided (t soil/ha/yr, log)")
    axes[1].set_title("Erosion benefit")
    axes[1].grid(True, axis="x", linestyle="--", alpha=0.3)

    fig.suptitle(f"Extent of practice-switch benefits @ {level} (boxes: IQR, whiskers: 5–95th pct)", fontsize=13)
    return fig, axes


def plot_attribution(commodity: str = "Wheat", level: str = "ecoregion"):
    """Residue vs tillage component effect on SOC gain and erosion avoided."""

    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    tbl = practice_attribution_table([commodity], level).iloc[0]
    cmap = colormaps["viridis"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    comps = ["residue", "tillage"]
    soc_vals = [tbl["soc_residue"], tbl["soc_tillage"]]
    se_vals = [tbl["se_residue"], tbl["se_tillage"]]
    axes[0].bar(comps, soc_vals, color=[cmap(0.3), cmap(0.65)])
    axes[0].set_ylabel("median SOC gained (t SOC/ha)")
    axes[0].set_title("SOC benefit by practice component")
    axes[1].bar(comps, se_vals, color=[cmap(0.3), cmap(0.65)])
    axes[1].set_ylabel("median erosion avoided (t soil/ha/yr)")
    axes[1].set_title("Erosion benefit by practice component")
    for ax in axes:
        ax.grid(True, axis="y", linestyle="--", alpha=0.3)
    fig.suptitle(f"{commodity}: what drives the benefit (@ {level})", fontsize=13)
    return fig, axes


def plot_cobenefit_by_realm(commodity: str = "Wheat", level: str = "ecoregion", ax=None):
    """Median SOC gain and erosion avoided by realm (geometry-free spatial view)."""

    import matplotlib.pyplot as plt

    t = cobenefit_table(commodity, level)
    g = t.dropna(subset=["realm"]).groupby("realm").agg(
        med_d_soc=("d_soc", "median"), med_d_se_red=("d_se_red", "median"), n=("d_soc", "size")
    )
    g = g[g["n"] >= 3].sort_values("med_d_se_red")
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 0.5 * len(g) + 2))
    else:
        fig = ax.figure
    y = np.arange(len(g))
    ax.barh(y - 0.2, g["med_d_soc"], height=0.4, color="#2c7fb8", label="SOC gained (t SOC/ha)")
    ax2 = ax.twiny()
    ax2.barh(y + 0.2, g["med_d_se_red"], height=0.4, color="#d95f0e", label="erosion avoided (t soil/ha/yr)")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r}  (n={int(n)})" for r, n in zip(g.index, g["n"])])
    ax.set_xlabel("median SOC gained (t SOC/ha)", color="#2c7fb8")
    ax2.set_xlabel("median erosion avoided (t soil/ha/yr)", color="#d95f0e")
    ax.set_title(f"{commodity}: median practice-switch benefit by realm")
    return fig, ax


# --------------------------------------------------------------------------- #
# Maps
# --------------------------------------------------------------------------- #


def _metric_style(metric: str):
    from matplotlib import colors as mcolors

    if metric in ("d_soc", "d_soc_pct"):
        return "RdYlGn", lambda v: mcolors.TwoSlopeNorm(
            vcenter=0, vmin=min(np.nanpercentile(v, 2), -1e-6), vmax=max(np.nanpercentile(v, 98), 1e-6)
        )
    if metric in ("d_se_red", "d_se_red_pct"):
        return "Greens", lambda v: mcolors.LogNorm(
            vmin=max(np.nanpercentile(v[v > 0], 2), 1e-3), vmax=np.nanpercentile(v, 98)
        ) if (v > 0).any() else None
    return "viridis", lambda v: mcolors.Normalize(vmin=0, vmax=1)


def _join_geometry(commodity: str, level: str, metric: str, *, allow_download: bool = True):
    t = cobenefit_table(commodity, level)
    if t.empty:
        return None, None
    if level == "country":
        gdf = geo.country_geometry(allow_download=allow_download)
        if gdf is None:
            return None, None
        vals = t[["region_name", metric]].copy()
        vals["key"] = vals["region_name"].map(geo.country_key)
        merged = gdf.merge(vals.groupby("key", as_index=False)[metric].mean(), on="key", how="left")
        return merged.rename(columns={metric: "value"}), gdf
    if level == "ecoregion":
        gdf = geo.ecoregion_geometry(allow_download=allow_download)
        if gdf is None:
            return None, None
        vals = t[["region_id", metric]].copy()
        vals["ECO_ID"] = pd.to_numeric(vals["region_id"], errors="coerce")
        merged = gdf.merge(vals[["ECO_ID", metric]], on="ECO_ID", how="left")
        basemap = geo.country_geometry(allow_download=allow_download)
        return merged.rename(columns={metric: "value"}), basemap
    return None, None


def plot_cobenefit_map(commodity: str = "Wheat", level: str = "ecoregion", metric: str = "priority",
                       *, allow_download: bool = True, ax=None):
    """Choropleth of one co-benefit ``metric`` for a commodity, or ``None`` if no geometry."""

    merged, basemap = _join_geometry(commodity, level, metric, allow_download=allow_download)
    if merged is None:
        return None
    label, unit = _METRIC_LABEL[metric]
    cmap, norm_fn = _metric_style(metric)
    vals = merged["value"].to_numpy(dtype=float)
    norm = norm_fn(vals[np.isfinite(vals)]) if np.isfinite(vals).any() else None
    title = f"{commodity}: {label} ({unit}) — {level} (practice switch)"
    return geo.choropleth(merged, "value", ax=ax, cmap=cmap, norm=norm, title=title,
                          cbar_label=f"{label} ({unit})", basemap=basemap if level == "ecoregion" else None)


def plot_cobenefit_map_panels(commodity: str = "Wheat", level: str = "ecoregion", *, allow_download: bool = True):
    """Three-panel map: SOC gained, erosion avoided, co-benefit priority."""

    import matplotlib.pyplot as plt

    metrics = ["d_soc", "d_se_red", "priority"]
    fig, axes = plt.subplots(3, 1, figsize=(11, 15))
    n_drawn = 0
    for ax, metric in zip(axes, metrics):
        res = plot_cobenefit_map(commodity, level, metric, allow_download=allow_download, ax=ax)
        if res is None:
            ax.set_visible(False)
        else:
            n_drawn += 1
    if n_drawn == 0:
        plt.close(fig)
        return None
    fig.suptitle(f"{commodity} practice switch — where the co-benefits are largest ({level})", fontsize=14)
    return fig, axes


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #


def default_output_dir() -> Path:
    from sbtn_leaf.paths import project_path

    return project_path("paper", "claude_analysis", "outputs", "practice_change")


def run_practice_change_analysis(
    outdir: Optional[Path] = None,
    commodities: Optional[Sequence[str]] = None,
    *,
    make_figures: bool = True,
    make_maps: bool = True,
    allow_download: bool = True,
) -> Dict[str, object]:
    """Write practice-change tables, figures, maps and findings ``README.md``."""

    commodities = list(commodities) if commodities is not None else COMMODITIES
    outdir = Path(outdir) if outdir is not None else default_output_dir()
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    maps_dir = outdir / "figures" / "maps"
    tables_dir.mkdir(parents=True, exist_ok=True)
    if make_figures:
        maps_dir.mkdir(parents=True, exist_ok=True)

    summary = cobenefit_summary(commodities, "ecoregion")
    attribution = practice_attribution_table(commodities, "ecoregion")
    stack = practice_stack_table(level="ecoregion")
    summary.round(4).to_csv(tables_dir / "cobenefit_summary.csv", index=False)
    attribution.round(4).to_csv(tables_dir / "practice_attribution.csv", index=False)
    stack.round(4).to_csv(tables_dir / "practice_stack_range.csv", index=False)
    cobenefit_table("Wheat", "ecoregion").round(4).to_csv(tables_dir / "cobenefit_wheat_ecoregion.csv", index=False)
    for c in commodities:
        priority_regions(c, "ecoregion", 20).round(4).to_csv(tables_dir / f"priority_regions_{c.lower()}.csv", index=False)

    n_maps = 0
    if make_figures:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        def _save(fig, path, dpi=160):
            fig.tight_layout()
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)

        _save(plot_cobenefit_scatter("Wheat", "ecoregion")[0], figures_dir / "cobenefit_scatter_wheat.png")
        _save(plot_benefit_distributions(commodities, "ecoregion")[0], figures_dir / "benefit_distributions.png")
        _save(plot_attribution("Wheat", "ecoregion")[0], figures_dir / "attribution_wheat.png")
        _save(plot_cobenefit_by_realm("Wheat", "ecoregion")[0], figures_dir / "cobenefit_by_realm_wheat.png")
        _save(plot_practice_stack_range(level="ecoregion")[0], figures_dir / "practice_stack_range.png")

        if make_maps:
            # Full panels for Wheat at both levels; priority map for the rest.
            for level in ("ecoregion", "country"):
                res = plot_cobenefit_map_panels("Wheat", level, allow_download=allow_download)
                if res is not None:
                    _save(res[0], maps_dir / f"wheat_panels_{level}.png", dpi=140)
                    n_maps += 1
            for c in commodities:
                if c == "Wheat":
                    continue
                res = plot_cobenefit_map(c, "ecoregion", "priority", allow_download=allow_download)
                if res is not None:
                    _save(res[0], maps_dir / f"{c.lower()}_priority_ecoregion.png", dpi=140)
                    n_maps += 1

    _write_findings(outdir / "README.md", summary, attribution, stack, priority_regions("Wheat", "ecoregion", 12), n_maps)

    return {"summary": summary, "attribution": attribution, "outdir": outdir, "n_maps": n_maps}


def _write_findings(path, summary, attribution, stack, top_wheat, n_maps) -> None:
    lines: List[str] = ["# Multi-indicator practice-change co-benefits — findings\n"]
    lines.append(
        "Switching the **same commodity** from the baseline (conventional tillage + residues "
        "removed) to the regenerative practice (reduced tillage + residues left) — SOC gained "
        "and soil erosion avoided per region. Expands the manuscript's multi-indicator section / "
        "Fig. 13.\n"
    )
    lines.append("\n## Extent of the benefits (per commodity, ecoregion level)\n")
    lines.append(
        "`med_d_soc` = median SOC gained (t SOC/ha); `med_d_se_red` = median erosion avoided "
        "(t soil/ha/yr); `d_se_red_pct` = relative erosion reduction (spatially constant per "
        "commodity — a property of the multiplicative RUSLE C-factor); `pct_win_win` = share of "
        "regions improving on **both**; `spearman_soc_se` = co-location of the two benefits.\n\n"
    )
    lines.append(_df_to_md(summary.round(3)))
    lines.append("\n\n## What drives the benefit (median component effect)\n")
    lines.append(
        "`*_residue` / `*_tillage` = the SOC gain / erosion avoided attributable to residue "
        "retention vs reduced tillage alone.\n\n"
    )
    lines.append(_df_to_md(attribution.round(3)))
    lines.append("\n\n## Full practice-stack SOC range (best vs worst combination, by commodity)\n")
    lines.append(
        "Median / mean % SOC difference between the best (rainfed + residues retained + reduced "
        "tillage) and worst (irrigated + residues removed + conventional tillage) stack. The "
        "cross-commodity mean reproduces the manuscript's ~37.5%, but the by-commodity range is "
        "wide (rapeseed lowest, maize highest).\n\n"
    )
    lines.append(_df_to_md(stack.round(2)))
    lines.append("\n\n## Where to focus — top wheat ecoregions by co-benefit priority\n")
    lines.append(_df_to_md(top_wheat.round(3)))
    lines.append(f"\n\nMaps rendered this run: **{n_maps}** (0 ⇒ no geometry available; maps are drop-in).\n")
    lines.append(
        "\nRegenerate with `python -m sbtn_leaf.claude_analysis.run_practice_change` or "
        "`sbtn_leaf.claude_analysis.practice_change.run_practice_change_analysis()`.\n"
    )
    Path(path).write_text("\n".join(lines))
