"""Generic spatial-aggregation comparison engine for SBTN-Land LEAFs.

This module generalises the soil-erosion-specific analysis in
:mod:`sbtn_leaf.claude_analysis.se_aggregation_analysis` so the *same* machinery
can compare **SOC**, **soil erosion** and **terrestrial acidification** LEAFs
across the three published spatial aggregation levels:

``country``     -- FAO ADM0 boundaries (coarsest)
``subcountry``  -- FAO ADM1 boundaries
``ecoregion``   -- WWF 2017 ecoregions (carry biome / realm grouping)

Every indicator is harmonised onto one tidy long schema (:data:`STD_COLS`) with a
canonical ``flow`` key and per-region ``leaf`` (mean), ``leaf_median`` and
``leaf_std`` columns.  Once harmonised, all summary statistics, variance
decompositions and figures are indicator-agnostic and parameterised only by the
:class:`IndicatorConfig` (which carries the unit, the human-readable flow labels
and the representative *focal* flows).

The design answers two questions, for each indicator and then *across* them:

* **Polygon size** -- how do the mean / median / standard deviation of the
  estimate change as polygons shrink from countries to ecoregions?
* **Ecoregion significance** -- do ecological (ecoregion / biome / realm)
  boundaries capture variation that political (country / admin-1) boundaries
  average away?

Statistics depend only on ``pandas`` / ``numpy``; the figures additionally use
``matplotlib``.  The core variance maths is reused from the soil-erosion module
so the two analyses cannot drift apart.
"""

from __future__ import annotations

import warnings
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

# Reuse the validated one-way ANOVA partition so SOC / acidification variance
# decompositions are computed by the *exact* same code as soil erosion.
from sbtn_leaf.claude_analysis.se_aggregation_analysis import variance_decomposition

__all__ = [
    "LEVELS",
    "STD_COLS",
    "IndicatorConfig",
    "load_long_level",
    "coverage_table",
    "validate_flow_coverage",
    "summary_table",
    "hierarchy_variance_table",
    "plot_cross_level_box",
    "plot_within_std_box",
    "plot_distribution_overlay",
    "plot_biome_box",
    "plot_sensitivity_heatmap",
    "plot_commodity_ranking",
    "variance_decomposition",
]

#: Aggregation levels, ordered coarse -> fine (an ordered proxy for polygon size).
LEVELS = ("country", "subcountry", "ecoregion")

#: The harmonised long schema shared by every indicator and level.
STD_COLS = [
    "level",
    "region_id",
    "region_name",
    "country_name",
    "biome",
    "realm",
    "flow",
    "flow_raw",
    "leaf",
    "leaf_median",
    "leaf_std",
]

#: Per-level region key + descriptive columns in the source CSVs.  ``ecoregion``
#: additionally carries the biome / realm grouping fields (joined when absent).
_LEVEL_KEYS: Dict[str, Dict[str, Optional[str]]] = {
    "country": {"id": "ADM0_NAME", "name": "ADM0_NAME", "country": "ADM0_NAME"},
    "subcountry": {"id": "ADM1_CODE", "name": "ADM1_NAME", "country": "ADM0_NAME"},
    "ecoregion": {"id": "ECO_ID", "name": "ECO_NAME", "country": None},
}


# --------------------------------------------------------------------------- #
# Indicator configuration
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class IndicatorConfig:
    """Everything indicator-specific the generic engine needs.

    Parameters
    ----------
    key, name, unit:
        Short identifier, display name and physical unit of the ``leaf`` value.
    flow_kind:
        ``"commodity"`` (a land use x management, e.g. SOC / soil erosion) or
        ``"pollutant"`` (a region-wide characterisation factor, e.g. the three
        acidifying gases).  Controls only wording in titles / findings.
    higher_is_better:
        ``True`` for SOC (more soil carbon is desirable), ``False`` for soil
        erosion and acidification (higher = worse).  Used only for narrative.
    focal_flows:
        Ordered ``canonical flow -> label`` map of the representative flows the
        figures and focal tables highlight.
    loader:
        Zero-argument callable returning the *raw* harmonised long frame for this
        indicator (all three levels concatenated, no-data rows kept).
    labeller:
        Maps a canonical ``flow`` key to a human-readable label.
    """

    key: str
    name: str
    unit: str
    flow_kind: str
    higher_is_better: bool
    focal_flows: "OrderedDict[str, str]"
    loader: Callable[[], pd.DataFrame]
    labeller: Callable[[str], str]

    # -- convenience ----------------------------------------------------------
    def label(self, flow: str) -> str:
        return self.labeller(flow)

    def load_harmonized(self, *, drop_na: bool = True) -> pd.DataFrame:
        """Return the harmonised long frame; drop no-data ``leaf`` rows by default."""

        raw = self.loader()
        if drop_na:
            raw = raw[raw["leaf"].notna()].reset_index(drop=True)
        return raw

    @property
    def direction_word(self) -> str:
        return "higher = more carbon stored" if self.higher_is_better else "higher = worse"


# --------------------------------------------------------------------------- #
# Loading & harmonisation (long "metric"/"variable" tables)
# --------------------------------------------------------------------------- #


def load_long_level(
    path: Path,
    level: str,
    *,
    pivot_col: str,
    value_map: Dict[str, str],
    canonical: Callable[[str], str],
    eco_meta: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Harmonise one long-format aggregation CSV onto :data:`STD_COLS`.

    Both the SOC tables (``variable`` in ``{leaf, leaf_median, leaf_std}``) and
    the acidification tables (``metric`` in ``{cf_mean, cf_median, cf_std}``) are
    region x flow x statistic long tables; this single reader pivots either onto
    one ``leaf`` / ``leaf_median`` / ``leaf_std`` row per region x flow.

    Parameters
    ----------
    path:
        CSV to read.
    level:
        One of :data:`LEVELS`.
    pivot_col:
        Column holding the statistic name (``"variable"`` or ``"metric"``).
    value_map:
        Maps the ``pivot_col`` values onto ``leaf`` / ``leaf_median`` /
        ``leaf_std``.
    canonical:
        Maps the raw ``flow_name`` onto the canonical ``flow`` key.
    eco_meta:
        Optional ``ECO_ID -> (ECO_NAME, BIOME_NAME, REALM)`` lookup, joined for
        ecoregion tables that lack those columns (the SOC ``v1.0`` files do).
    """

    if level not in _LEVEL_KEYS:
        raise ValueError(f"Unknown level {level!r}; expected one of {LEVELS}.")

    df = pd.read_csv(path)
    keys = _LEVEL_KEYS[level]
    id_col = keys["id"]

    # groupby().unstack() keeps every existing region x flow cell -- including the
    # all-NaN "no data" cells -- without exploding into a full cartesian product,
    # so coverage stays directly comparable across levels and indicators.
    wide = (
        df.groupby([id_col, "flow_name", pivot_col])["value"]
        .first()
        .unstack(pivot_col)
        .rename(columns=value_map)
        .reset_index()
    )

    out = pd.DataFrame(index=wide.index)
    out["level"] = level
    out["region_id"] = wide[id_col].astype("string")
    out["flow_raw"] = wide["flow_name"].astype("string")
    out["flow"] = wide["flow_name"].map(canonical).astype("string")
    for col in ("leaf", "leaf_median", "leaf_std"):
        out[col] = wide[col] if col in wide.columns else np.nan

    if level == "country":
        out["region_name"] = wide[id_col]
        out["country_name"] = wide[id_col]
        out["biome"] = pd.NA
        out["realm"] = pd.NA
    elif level == "subcountry":
        meta = df[[id_col, keys["name"], keys["country"]]].drop_duplicates(id_col)
        merged = wide.merge(meta, on=id_col, how="left")
        out["region_name"] = merged[keys["name"]]
        out["country_name"] = merged[keys["country"]]
        out["biome"] = pd.NA
        out["realm"] = pd.NA
    else:  # ecoregion
        if eco_meta is not None:
            meta = eco_meta.drop_duplicates("ECO_ID")
        else:
            meta = df[["ECO_ID", "ECO_NAME", "BIOME_NAME", "REALM"]].drop_duplicates("ECO_ID")
        merged = wide.merge(meta, on="ECO_ID", how="left")
        out["region_name"] = merged["ECO_NAME"]
        out["country_name"] = pd.NA
        out["biome"] = merged["BIOME_NAME"]
        out["realm"] = merged["REALM"]

    return out.reindex(columns=STD_COLS)


def coverage_table(df_raw: pd.DataFrame) -> pd.DataFrame:
    """Per-level region / flow counts and ``leaf`` non-null coverage."""

    rows = []
    for level in LEVELS:
        sub = df_raw[df_raw["level"] == level]
        rows.append(
            {
                "level": level,
                "n_regions": sub["region_id"].nunique(),
                "n_flows": sub["flow"].nunique(),
                "n_rows": len(sub),
                "coverage_pct": 100.0 * sub["leaf"].notna().mean() if len(sub) else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def validate_flow_coverage(df_raw: pd.DataFrame, focal_flows: Iterable[str] = ()) -> pd.DataFrame:
    """Report how canonical flows overlap across the three levels.

    Returns one row per canonical flow with a boolean per level plus an
    ``in_all`` flag; warns if any ``focal_flows`` key is absent from a level.
    """

    present = {level: set(df_raw.loc[df_raw["level"] == level, "flow"].unique()) for level in LEVELS}
    all_flows = sorted(set().union(*present.values())) if present else []
    table = pd.DataFrame(
        {level: [f in present[level] for f in all_flows] for level in LEVELS},
        index=all_flows,
    )
    table["in_all"] = table[list(LEVELS)].all(axis=1) if len(table) else []
    table.index.name = "flow"

    missing = [f for f in focal_flows if f in table.index and not table.loc[f, "in_all"]]
    if missing:
        warnings.warn(f"Focal flows missing from at least one level: {missing}")
    return table


# --------------------------------------------------------------------------- #
# Summary statistics
# --------------------------------------------------------------------------- #


def summary_table(
    cfg: IndicatorConfig,
    df: pd.DataFrame,
    flows: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Per ``flow`` x ``level`` summary of the regional ``leaf`` distribution.

    Columns: ``n_regions, mean_leaf, median_leaf, std_leaf`` (the *between-region*
    spread), ``cv``, ``p10, p90, max_leaf`` and ``mean_within_std`` (mean of the
    per-region ``leaf_std``, i.e. typical *within-region* heterogeneity).
    """

    if flows is not None:
        df = df[df["flow"].isin(list(flows))]

    grouped = df.groupby(["flow", "level"], observed=True)
    leaf = grouped["leaf"].agg(
        n_regions="count",
        mean_leaf="mean",
        median_leaf="median",
        std_leaf="std",
        p10=lambda s: s.quantile(0.10),
        p90=lambda s: s.quantile(0.90),
        max_leaf="max",
    )
    leaf["mean_within_std"] = grouped["leaf_std"].mean()
    leaf["cv"] = leaf["std_leaf"] / leaf["mean_leaf"]
    out = leaf.reset_index()
    out["level"] = pd.Categorical(out["level"], categories=list(LEVELS), ordered=True)
    out["flow_label"] = out["flow"].map(cfg.label)
    out["indicator"] = cfg.key
    return out.sort_values(["flow", "level"]).reset_index(drop=True)


def hierarchy_variance_table(
    cfg: IndicatorConfig,
    df: pd.DataFrame,
    flows: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Nested variance partitions that quantify aggregation-driven information loss.

    For each flow:

    * **within-country** -- partition the *subcountry* ``leaf`` values by their
      parent country.  ``frac_within`` is the share of variance a single national
      LEAF hides inside its admin-1 units (the cost of political coarsening).
    * **between-biome** / **between-realm** -- partition the *ecoregion* ``leaf``
      values by biome / realm.  ``eta2_between`` is how much an ecological grouping
      explains (ecoregion significance).
    """

    flows = list(flows) if flows is not None else list(cfg.focal_flows)
    sub = df[df["level"] == "subcountry"]
    eco = df[df["level"] == "ecoregion"]
    rows = []
    for flow in flows:
        s = sub[sub["flow"] == flow]
        e = eco[eco["flow"] == flow]
        wc = variance_decomposition(s["leaf"], s["country_name"])
        bb = variance_decomposition(e["leaf"], e["biome"])
        br = variance_decomposition(e["leaf"], e["realm"])
        rows.append(
            {
                "indicator": cfg.key,
                "flow": flow,
                "flow_label": cfg.label(flow),
                "subcty_within_country_frac": wc["frac_within"],
                "subcty_n_units": wc["n_obs"],
                "eco_between_biome_eta2": bb["eta2_between"],
                "eco_between_realm_eta2": br["eta2_between"],
                "eco_n_units": bb["n_obs"],
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Plot helpers (matplotlib only; every function returns (fig, ax))
# --------------------------------------------------------------------------- #


def _level_colors() -> List[tuple]:
    from matplotlib import colormaps

    cmap = colormaps["viridis"]
    return [cmap(x) for x in (0.15, 0.5, 0.85)]


def _positive(s: pd.Series) -> np.ndarray:
    s = pd.to_numeric(s, errors="coerce").dropna()
    return s[s > 0].to_numpy()


def _ylabel(cfg: IndicatorConfig, value: str) -> str:
    if value == "leaf_std":
        return f"Within-region SD ({cfg.unit})"
    return f"{cfg.name} ({cfg.unit})"


def plot_cross_level_box(cfg: IndicatorConfig, df: pd.DataFrame, flow: str, value: str = "leaf", ax=None):
    """Boxplots of a flow's regional ``value`` at each aggregation level."""

    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    data, labels = [], []
    for level in LEVELS:
        arr = _positive(df.loc[(df["flow"] == flow) & (df["level"] == level), value])
        data.append(arr)
        labels.append(f"{level}\n(n={len(arr)})")

    bp = ax.boxplot(
        data, whis=(5, 95), showfliers=False, patch_artist=True, widths=0.6, medianprops={"color": "black"}
    )
    for patch, color in zip(bp["boxes"], _level_colors()):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    if all((a > 0).all() for a in data if len(a)):
        ax.set_yscale("log")
    ax.set_xticks(range(1, len(LEVELS) + 1))
    ax.set_xticklabels(labels)
    ax.set_ylabel(_ylabel(cfg, value))
    ax.set_title(
        f"{cfg.label(flow)} — {cfg.name} by aggregation level\n(boxes: IQR, whiskers: 5–95th pct)"
    )
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    return fig, ax


def plot_within_std_box(cfg: IndicatorConfig, df: pd.DataFrame, flow: str, ax=None):
    """Boxplots of per-region ``leaf_std`` (within-region heterogeneity) by level."""

    return plot_cross_level_box(cfg, df, flow, value="leaf_std", ax=ax)


def plot_distribution_overlay(cfg: IndicatorConfig, df: pd.DataFrame, flow: str, bins: int = 40, ax=None):
    """Overlaid log-x histograms of the regional ``leaf`` distribution per level."""

    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 6))
    else:
        fig = ax.figure

    series = {level: _positive(df.loc[(df["flow"] == flow) & (df["level"] == level), "leaf"]) for level in LEVELS}
    allv = (
        np.concatenate([v for v in series.values() if len(v)])
        if any(len(v) for v in series.values())
        else np.array([1.0])
    )
    edges = np.logspace(np.log10(max(allv.min(), 1e-6)), np.log10(allv.max()), bins)
    for (level, arr), color in zip(series.items(), _level_colors()):
        if len(arr):
            ax.hist(arr, bins=edges, density=True, histtype="step", linewidth=2, color=color, label=f"{level} (n={len(arr)})")
            ax.axvline(np.median(arr), color=color, linestyle=":", linewidth=1.5)
    ax.set_xscale("log")
    ax.set_xlabel(f"{cfg.name} ({cfg.unit}) — log scale")
    ax.set_ylabel("Density")
    ax.set_title(f"{cfg.label(flow)} — regional distribution by level\n(dotted lines: medians)")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)
    return fig, ax


def plot_biome_box(cfg: IndicatorConfig, df: pd.DataFrame, flow: str, ax=None):
    """Horizontal boxplots of ecoregion ``leaf`` by WWF biome for one flow (ranked)."""

    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    sub = df[(df["level"] == "ecoregion") & (df["flow"] == flow)].copy()
    sub = sub[sub["biome"].notna()]
    order = sub.groupby("biome", observed=True)["leaf"].median().sort_values().index.tolist()
    data = [_positive(sub.loc[sub["biome"] == b, "leaf"]) for b in order]
    counts = [len(d) for d in data]

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 0.5 * len(order) + 2))
    else:
        fig = ax.figure
    bp = ax.boxplot(data, orientation="horizontal", whis=(5, 95), showfliers=False, patch_artist=True, medianprops={"color": "black"})
    cmap = colormaps["YlGn" if cfg.higher_is_better else "YlOrBr"]
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cmap(0.2 + 0.7 * i / max(len(order) - 1, 1)))
        patch.set_alpha(0.9)
    ax.set_yticks(range(1, len(order) + 1))
    ax.set_yticklabels([f"{b}  (n={c})" for b, c in zip(order, counts)])
    if all((d > 0).all() for d in data if len(d)):
        ax.set_xscale("log")
    ax.set_xlabel(f"{cfg.name} ({cfg.unit})")
    ax.set_title(f"{cfg.label(flow)} — ecoregion {cfg.name} by biome")
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    return fig, ax


def plot_sensitivity_heatmap(cfg: IndicatorConfig, df: pd.DataFrame, flows: Optional[Iterable[str]] = None, metric: str = "mean_leaf", ax=None):
    """Heatmap of a summary ``metric`` across focal flows (rows) x levels (cols)."""

    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors

    flows = list(flows) if flows is not None else list(cfg.focal_flows)
    summ = summary_table(cfg, df, flows=flows)
    mat = summ.pivot(index="flow", columns="level", values=metric).reindex(index=flows, columns=list(LEVELS))

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 0.6 * len(flows) + 2))
    else:
        fig = ax.figure
    values = mat.to_numpy(dtype=float)
    norm = mcolors.LogNorm(vmin=np.nanmin(values[values > 0]), vmax=np.nanmax(values)) if (values > 0).any() else None
    im = ax.imshow(values, aspect="auto", cmap="viridis", norm=norm)
    ax.set_xticks(range(len(LEVELS)))
    ax.set_xticklabels(LEVELS)
    ax.set_yticks(range(len(flows)))
    ax.set_yticklabels([cfg.label(f) for f in flows])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2g}", ha="center", va="center", color="white", fontsize=8)
    ax.set_title(f"{cfg.name} {metric} ({cfg.unit}) by flow x level")
    fig.colorbar(im, ax=ax, label=f"{metric} ({cfg.unit})")
    return fig, ax


def plot_commodity_ranking(cfg: IndicatorConfig, df: pd.DataFrame, level: str = "ecoregion", flows: Optional[Iterable[str]] = None, ax=None):
    """Horizontal boxplots of ``leaf`` per focal flow at one level, ranked by median."""

    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    flows = list(flows) if flows is not None else list(cfg.focal_flows)
    sub = df[(df["level"] == level) & (df["flow"].isin(flows))]
    order = sub.groupby("flow", observed=True)["leaf"].median().sort_values().index.tolist()
    data = [_positive(sub.loc[sub["flow"] == f, "leaf"]) for f in order]

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 0.6 * len(order) + 2))
    else:
        fig = ax.figure
    bp = ax.boxplot(data, orientation="horizontal", whis=(5, 95), showfliers=False, patch_artist=True, medianprops={"color": "black"})
    cmap = colormaps["viridis"]
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cmap(i / max(len(order) - 1, 1)))
        patch.set_alpha(0.85)
    ax.set_yticks(range(1, len(order) + 1))
    ax.set_yticklabels([cfg.label(f) for f in order])
    if all((d > 0).all() for d in data if len(d)):
        ax.set_xscale("log")
    ax.set_xlabel(f"{cfg.name} ({cfg.unit}) — log scale")
    ax.set_title(f"{cfg.name} ranking at {level} level\n(boxes: IQR, whiskers: 5–95th pct)")
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    return fig, ax
