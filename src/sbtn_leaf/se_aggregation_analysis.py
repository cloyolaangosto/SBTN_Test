"""Compare soil-erosion LEAF statistics across spatial aggregation levels.

The soil-erosion LEAFs are published as three aggregation tables that clip the
*same* 25 km erosion raster to three different polygon sets:

``se_country_clipped_2.csv``     -- FAO country boundaries (ADM0)
``se_subcountry_clipped_2.csv``  -- FAO admin-1 boundaries (ADM1)
``se_ecoregions_clipped_2.csv``  -- WWF 2017 ecoregions

This module harmonises those tables onto a common schema and a common
``flow_name`` key, then provides summary statistics, variance decompositions and
matplotlib figures used to study two questions:

* **Polygon size** -- how do the mean / median / standard deviation of the
  erosion estimate change as polygons shrink from countries to ecoregions?
* **Ecoregion significance** -- do ecological (ecoregion / biome) boundaries
  capture erosion variation that political boundaries average away?

The three files use *different* ``flow_name`` conventions (the ecoregion file
mostly uses descriptive names such as ``Rainfed_Wheat_residues_removed_from_the_field``
and ``Broadleaf_Deciduous_Tropical`` while the country / subcountry files use
short codes such as ``Wheat_rf_roff`` and ``BRDC_Tropical``).  :func:`canonical_flow`
maps every convention onto a single ``commodity|water|residue|tillage`` key (or a
``BRDC_/NEEV_`` land-cover key), which yields a complete 106-flow match across the
three levels.

Statistics and plots depend only on ``pandas`` / ``numpy`` / ``scipy`` /
``matplotlib``.  The geometry-dependent maps lazily import ``geopandas`` and the
boundary files; when those are unavailable (e.g. DVC data not pulled) the map
helpers return ``None`` with a warning so the rest of the analysis still runs.
"""

from __future__ import annotations

import math
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from sbtn_leaf.paths import leaf_path

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

#: Prefix the ecoregion source rasters (and therefore ecoregion ``flow_name``s)
#: carry; stripped during canonicalisation.
ECO_PREFIX = "se_rate_25km_clipped_"

#: Aggregation levels, ordered coarse -> fine (an ordered proxy for polygon size).
LEVELS: Tuple[str, str, str] = ("country", "subcountry", "ecoregion")

#: Soil-erosion aggregation tables, keyed by level.
FILES: Dict[str, str] = {
    "country": "se_country_clipped_2.csv",
    "subcountry": "se_subcountry_clipped_2.csv",
    "ecoregion": "se_ecoregions_clipped_2.csv",
}

#: Physical unit of the ``leaf`` value.
UNIT = "t soil/ha/yr"

#: Representative focal commodities (canonical key -> display label).  A spread of
#: major traded commodities plus two land covers, all on a rainfed / conventional
#: baseline so they are directly comparable.  Every key is present at all three
#: levels (validated by :func:`validate_flow_coverage`).
FOCAL_FLOWS: "OrderedDict[str, str]" = OrderedDict(
    [
        ("Wheat|rf|roff|ct", "Wheat (rainfed)"),
        ("Maize|rf|roff|ct", "Maize (rainfed)"),
        ("Soybeans|rf|na|ct", "Soybeans (rainfed)"),
        ("Oil_palm|rf|na|ct", "Oil palm (rainfed)"),
        ("Coffee|rf|na|ct", "Coffee (rainfed)"),
        ("Sugarcane|rf|na|ct", "Sugarcane (rainfed)"),
        ("Cotton|rf|na|ct", "Cotton (rainfed)"),
        ("Grassland", "Grassland"),
        ("BRDC_Tropical", "Broadleaf-decid. forest (tropical)"),
    ]
)

#: Geometry join keys and boundary-file region columns, by level.
_GEOMETRY = {
    "country": {"join": "ADM0_NAME", "region": "ADM0_NAME"},
    "subcountry": {"join": "ADM1_CODE", "region": "ADM1_NAME"},
    "ecoregion": {"join": "ECO_ID", "region": "ECO_NAME"},
}

_WATER = {"irr": "irrigated", "rf": "rainfed"}
_RESIDUE = {"roff": "residues removed", "ron": "residues retained"}

# --------------------------------------------------------------------------- #
# Flow-name canonicalisation
# --------------------------------------------------------------------------- #


def strip_prefix(raw: str) -> str:
    """Remove the ecoregion raster prefix / ``.tif`` suffix from ``raw``."""

    s = str(raw)
    if s.startswith(ECO_PREFIX):
        s = s[len(ECO_PREFIX) :]
    if s.endswith(".tif"):
        s = s[:-4]
    return s


def canonical_flow(raw: str) -> str:
    """Map any ``flow_name`` spelling onto a single canonical key.

    Crops collapse to ``"{Commodity}|{water}|{residue}|{tillage}"`` with
    ``water`` in ``{irr, rf}``, ``residue`` in ``{roff, ron, na}`` and
    ``tillage`` in ``{ct, rt}``.  Land covers collapse to ``"BRDC_{zone}"`` /
    ``"NEEV_{zone}"`` or the bare ``"Grassland"`` / ``"Urban"`` token.
    """

    s = strip_prefix(raw)

    if s in ("Grassland", "Urban"):
        return s
    if s.startswith("Broadleaf_Deciduous_"):
        return "BRDC_" + s[len("Broadleaf_Deciduous_") :]
    if s.startswith("Needleleaf_Evergreen_"):
        return "NEEV_" + s[len("Needleleaf_Evergreen_") :]
    if s.startswith("BRDC_") or s.startswith("NEEV_"):
        return s

    # Descriptive crop, e.g. ``Rainfed_Wheat_residues_removed_from_the_field``.
    if s.startswith("Irrigated_") or s.startswith("Rainfed_"):
        water = "irr" if s.startswith("Irrigated_") else "rf"
        rest = s.split("_", 1)[1]
        residue = "na"
        if rest.endswith("_residues_left_on_the_field"):
            residue, crop = "ron", rest[: -len("_residues_left_on_the_field")]
        elif rest.endswith("_residues_removed_from_the_field"):
            residue, crop = "roff", rest[: -len("_residues_removed_from_the_field")]
        else:
            crop = rest
        crop = crop[:1].upper() + crop[1:]
        return f"{crop}|{water}|{residue}|ct"

    # Short crop code, e.g. ``Wheat_rf_roff_rt`` / ``Soybeans_rf`` / ``cotton_rf``.
    marker = "_irr" if "_irr" in s else ("_rf" if "_rf" in s else None)
    if marker is not None:
        idx = s.find(marker)
        crop = s[:idx]
        water = "irr" if marker == "_irr" else "rf"
        toks = s[idx + len(marker) :].lstrip("_").split("_")
        toks = [t for t in toks if t]
        residue = "roff" if "roff" in toks else ("ron" if "ron" in toks else "na")
        tillage = "rt" if "rt" in toks else "ct"
        crop = crop[:1].upper() + crop[1:]
        return f"{crop}|{water}|{residue}|{tillage}"

    return f"UNPARSED:{s}"


def flow_label(canon: str) -> str:
    """Return a human-readable label for a canonical flow key."""

    if canon in FOCAL_FLOWS:
        return FOCAL_FLOWS[canon]
    if canon in ("Grassland", "Urban"):
        return canon
    if canon.startswith("BRDC_"):
        return "Broadleaf-decid. " + canon[5:].replace("_", " ")
    if canon.startswith("NEEV_"):
        return "Needleleaf-everg. " + canon[5:].replace("_", " ")
    if "|" in canon:
        crop, water, residue, tillage = canon.split("|")
        parts = [_WATER.get(water, water)]
        if residue != "na":
            parts.append(_RESIDUE[residue])
        if tillage == "rt":
            parts.append("reduced till.")
        return f"{crop.replace('_', ' ')} ({', '.join(parts)})"
    return canon


# --------------------------------------------------------------------------- #
# Loading & harmonisation
# --------------------------------------------------------------------------- #

_STD_COLS = [
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


def load_level(level: str, *, drop_na: bool = True) -> pd.DataFrame:
    """Load one aggregation table, harmonised to the common schema.

    Parameters
    ----------
    level:
        One of :data:`LEVELS`.
    drop_na:
        When ``True`` (default) rows whose ``leaf`` (mean) value is missing -- i.e.
        regions where the commodity has no modelled erosion -- are dropped.

    Returns
    -------
    pandas.DataFrame
        Columns: ``level, region_id, region_name, country_name, biome, realm,
        flow, flow_raw, leaf, leaf_median, leaf_std``.  ``flow`` is the canonical
        key; ``biome`` / ``realm`` are populated only for the ecoregion level.
    """

    if level not in FILES:
        raise ValueError(f"Unknown level {level!r}; expected one of {LEVELS}.")

    df = pd.read_csv(leaf_path("soil_erosion", FILES[level]))
    df["flow_raw"] = df["flow_name"].map(strip_prefix)
    df["flow"] = df["flow_name"].map(canonical_flow)

    if level == "ecoregion":
        # Reshape long (metric/value) -> wide on the true keys (ECO_ID, flow) only.
        # groupby().unstack() keeps every existing region x flow cell, including the
        # all-NaN "no data" cells, without exploding into a full cartesian product,
        # so coverage stays directly comparable to the country/subcountry grids.
        metrics = (
            df.groupby(["ECO_ID", "flow", "metric"])["value"]
            .first()
            .unstack("metric")
            .rename(columns={"cf_mean": "leaf", "cf_median": "leaf_median", "cf_std": "leaf_std"})
            .reset_index()
        )
        meta = df[["ECO_ID", "ECO_NAME", "BIOME_NAME", "REALM", "flow", "flow_raw"]].drop_duplicates(
            subset=["ECO_ID", "flow"]
        )
        wide = meta.merge(metrics, on=["ECO_ID", "flow"], how="left")
        wide["level"] = "ecoregion"
        wide["region_id"] = wide["ECO_ID"].astype("string")
        wide["region_name"] = wide["ECO_NAME"]
        wide["country_name"] = pd.NA
        wide["biome"] = wide["BIOME_NAME"]
        wide["realm"] = wide["REALM"]
        out = wide
    else:
        df["level"] = level
        if level == "country":
            df["region_id"] = df["ADM0_NAME"].astype("string")
            df["region_name"] = df["ADM0_NAME"]
            df["country_name"] = df["ADM0_NAME"]
        else:  # subcountry
            df["region_id"] = df["ADM1_CODE"].astype("string")
            df["region_name"] = df["ADM1_NAME"]
            df["country_name"] = df["ADM0_NAME"]
        df["biome"] = pd.NA
        df["realm"] = pd.NA
        out = df

    out = out.reindex(columns=_STD_COLS)
    if drop_na:
        out = out[out["leaf"].notna()]
    return out.reset_index(drop=True)


def load_harmonized(*, drop_na: bool = True) -> pd.DataFrame:
    """Concatenate all three levels into one tidy long DataFrame."""

    return pd.concat(
        (load_level(level, drop_na=drop_na) for level in LEVELS),
        ignore_index=True,
    )


def validate_flow_coverage(df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Report how canonical flows overlap across the three levels.

    Returns one row per canonical flow with a boolean for each level and an
    ``in_all`` flag.  Also warns if any :data:`FOCAL_FLOWS` key is missing
    somewhere.
    """

    if df is None:
        df = load_harmonized(drop_na=False)
    present = {
        level: set(df.loc[df["level"] == level, "flow"].unique()) for level in LEVELS
    }
    all_flows = sorted(set().union(*present.values()))
    table = pd.DataFrame(
        {level: [f in present[level] for f in all_flows] for level in LEVELS},
        index=all_flows,
    )
    table["in_all"] = table[list(LEVELS)].all(axis=1)
    table.index.name = "flow"

    missing = [f for f in FOCAL_FLOWS if not table.loc[f, "in_all"]] if len(table) else []
    if missing:
        warnings.warn(f"Focal flows missing from at least one level: {missing}")
    return table


def coverage_table(df_raw: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Per-level region counts and ``leaf`` non-null coverage."""

    if df_raw is None:
        df_raw = load_harmonized(drop_na=False)
    rows = []
    for level in LEVELS:
        sub = df_raw[df_raw["level"] == level]
        rows.append(
            {
                "level": level,
                "n_regions": sub["region_id"].nunique(),
                "n_flows": sub["flow"].nunique(),
                "n_rows": len(sub),
                "coverage_pct": 100.0 * sub["leaf"].notna().mean(),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Summary statistics
# --------------------------------------------------------------------------- #


def summary_table(df: pd.DataFrame, flows: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """Per ``flow`` x ``level`` summary of the regional ``leaf`` distribution.

    Columns: ``n_regions, mean_leaf, median_leaf, std_leaf`` (the *between-region*
    spread), ``cv`` (coefficient of variation), ``p10, p90, max_leaf`` and
    ``mean_within_std`` (mean of the per-region ``leaf_std``, i.e. typical
    *within-region* heterogeneity).
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
    out["flow_label"] = out["flow"].map(flow_label)
    return out.sort_values(["flow", "level"]).reset_index(drop=True)


def variance_decomposition(values: Sequence[float], groups: Sequence) -> Dict[str, float]:
    """One-way ANOVA variance partition of ``values`` by ``groups``.

    Returns the total / between-group / within-group sums of squares together
    with ``eta2_between`` (fraction of variance explained by the grouping) and
    ``frac_within`` (the residual fraction hidden *inside* the groups).
    """

    d = pd.DataFrame({"v": pd.to_numeric(pd.Series(values), errors="coerce"), "g": list(groups)})
    d = d.dropna(subset=["v", "g"])
    n = len(d)
    n_groups = d["g"].nunique()
    if n < 2 or n_groups < 2:
        return {
            "n_obs": n,
            "n_groups": int(n_groups),
            "ss_total": float("nan"),
            "ss_between": float("nan"),
            "ss_within": float("nan"),
            "eta2_between": float("nan"),
            "frac_within": float("nan"),
        }
    grand = d["v"].mean()
    gmean = d.groupby("g")["v"].transform("mean")
    ss_total = float(((d["v"] - grand) ** 2).sum())
    ss_between = float(((gmean - grand) ** 2).sum())
    ss_within = float(((d["v"] - gmean) ** 2).sum())
    eta2 = ss_between / ss_total if ss_total > 0 else float("nan")
    return {
        "n_obs": n,
        "n_groups": int(n_groups),
        "ss_total": ss_total,
        "ss_between": ss_between,
        "ss_within": ss_within,
        "eta2_between": eta2,
        "frac_within": (ss_within / ss_total if ss_total > 0 else float("nan")),
    }


def hierarchy_variance_table(
    df: pd.DataFrame, flows: Optional[Iterable[str]] = None
) -> pd.DataFrame:
    """Nested variance partitions that quantify aggregation-driven information loss.

    For each focal flow:

    * **within-country** -- partition the *subcountry* ``leaf`` values by their
      parent country.  ``frac_within`` is the share of erosion variance a single
      national LEAF hides inside its admin-1 units (political coarsening cost).
    * **between-biome** / **between-realm** -- partition the *ecoregion* ``leaf``
      values by biome / realm.  ``eta2_between`` is how much an ecological
      grouping explains (ecoregion significance).
    """

    flows = list(flows) if flows is not None else list(FOCAL_FLOWS)
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
                "flow": flow,
                "flow_label": flow_label(flow),
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


def plot_cross_level_box(df: pd.DataFrame, flow: str, value: str = "leaf", ax=None):
    """Boxplots of a focal flow's regional ``value`` at each aggregation level."""

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
    ax.set_yscale("log")
    ax.set_xticks(range(1, len(LEVELS) + 1))
    ax.set_xticklabels(labels)
    ylab = f"Soil erosion ({UNIT})" if value == "leaf" else f"Within-region SD ({UNIT})"
    ax.set_ylabel(ylab)
    ax.set_title(f"{flow_label(flow)} — {value} by aggregation level\n(boxes: IQR, whiskers: 5–95th pct)")
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    return fig, ax


def plot_within_std_box(df: pd.DataFrame, flow: str, ax=None):
    """Boxplots of per-region ``leaf_std`` (within-region heterogeneity) by level."""

    return plot_cross_level_box(df, flow, value="leaf_std", ax=ax)


def plot_distribution_overlay(df: pd.DataFrame, flow: str, bins: int = 40, ax=None):
    """Overlaid log-x histograms of the regional ``leaf`` distribution per level."""

    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 6))
    else:
        fig = ax.figure

    series = {level: _positive(df.loc[(df["flow"] == flow) & (df["level"] == level), "leaf"]) for level in LEVELS}
    allv = np.concatenate([v for v in series.values() if len(v)]) if any(len(v) for v in series.values()) else np.array([1.0])
    edges = np.logspace(np.log10(max(allv.min(), 1e-3)), np.log10(allv.max()), bins)
    for (level, arr), color in zip(series.items(), _level_colors()):
        if len(arr):
            ax.hist(arr, bins=edges, density=True, histtype="step", linewidth=2, color=color, label=f"{level} (n={len(arr)})")
            ax.axvline(np.median(arr), color=color, linestyle=":", linewidth=1.5)
    ax.set_xscale("log")
    ax.set_xlabel(f"Soil erosion ({UNIT}) — log scale")
    ax.set_ylabel("Density")
    ax.set_title(f"{flow_label(flow)} — regional distribution by level\n(dotted lines: medians)")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)
    return fig, ax


def plot_commodity_ranking(df: pd.DataFrame, level: str = "ecoregion", flows: Optional[Iterable[str]] = None, ax=None):
    """Horizontal boxplots of ``leaf`` per focal commodity at one level, ranked by median."""

    import matplotlib.pyplot as plt

    flows = list(flows) if flows is not None else list(FOCAL_FLOWS)
    sub = df[(df["level"] == level) & (df["flow"].isin(flows))]
    order = sub.groupby("flow", observed=True)["leaf"].median().sort_values().index.tolist()
    data = [_positive(sub.loc[sub["flow"] == f, "leaf"]) for f in order]

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 0.6 * len(order) + 2))
    else:
        fig = ax.figure
    bp = ax.boxplot(data, orientation="horizontal", whis=(5, 95), showfliers=False, patch_artist=True, medianprops={"color": "black"})
    from matplotlib import colormaps

    cmap = colormaps["viridis"]
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cmap(i / max(len(order) - 1, 1)))
        patch.set_alpha(0.85)
    ax.set_yticks(range(1, len(order) + 1))
    ax.set_yticklabels([flow_label(f) for f in order])
    ax.set_xscale("log")
    ax.set_xlabel(f"Soil erosion ({UNIT}) — log scale")
    ax.set_title(f"Commodity erosion ranking at {level} level\n(boxes: IQR, whiskers: 5–95th pct)")
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    return fig, ax


def plot_biome_box(df: pd.DataFrame, flow: str, ax=None):
    """Horizontal boxplots of ecoregion ``leaf`` by biome for one flow (ranked)."""

    import matplotlib.pyplot as plt

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
    from matplotlib import colormaps

    cmap = colormaps["YlOrBr"]
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cmap(0.2 + 0.7 * i / max(len(order) - 1, 1)))
        patch.set_alpha(0.9)
    ax.set_yticks(range(1, len(order) + 1))
    ax.set_yticklabels([f"{b}  (n={c})" for b, c in zip(order, counts)])
    ax.set_xscale("log")
    ax.set_xlabel(f"Soil erosion ({UNIT}) — log scale")
    ax.set_title(f"{flow_label(flow)} — ecoregion erosion by biome")
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    return fig, ax


def plot_sensitivity_heatmap(df: pd.DataFrame, flows: Optional[Iterable[str]] = None, metric: str = "mean_leaf", ax=None):
    """Heatmap of a summary ``metric`` across focal flows (rows) x levels (cols)."""

    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors

    flows = list(flows) if flows is not None else list(FOCAL_FLOWS)
    summ = summary_table(df, flows=flows)
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
    ax.set_yticklabels([flow_label(f) for f in flows])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", color="white", fontsize=8)
    ax.set_title(f"{metric} ({UNIT}) by commodity x level")
    fig.colorbar(im, ax=ax, label=metric)
    return fig, ax


# --------------------------------------------------------------------------- #
# Maps (lazy geopandas; activate when boundary geometry is available)
# --------------------------------------------------------------------------- #


def load_boundaries(level: str):
    """Return the boundary GeoDataFrame for ``level`` or ``None`` if unavailable.

    Reuses :mod:`sbtn_leaf.data_loader` for country / ecoregion boundaries and the
    conventional ``CountryLayers/SubCountry_Level1/g2015_2014_1.shp`` path for
    subcountry.  Missing geometry or a missing ``geopandas`` install degrade to a
    warning + ``None`` rather than an exception, so the statistical analysis is
    unaffected.
    """

    try:
        import geopandas as gpd  # noqa: F401  (import guard)

        if level == "country":
            from sbtn_leaf.data_loader import get_country_boundaries

            return get_country_boundaries()
        if level == "ecoregion":
            from sbtn_leaf.data_loader import get_ecoregions_shapefile

            return get_ecoregions_shapefile()
        if level == "subcountry":
            from sbtn_leaf.paths import data_path

            path = data_path("CountryLayers", "SubCountry_Level1", "g2015_2014_1.shp")
            if not path.exists():
                raise FileNotFoundError(path)
            return gpd.read_file(path)
        raise ValueError(f"Unknown level {level!r}.")
    except ValueError:
        raise
    except Exception as exc:  # pragma: no cover - env dependent (missing/corrupt geometry, driver errors)
        warnings.warn(f"Boundary geometry for '{level}' unavailable ({type(exc).__name__}: {exc}).")
        return None


def choropleth_geodata(df: pd.DataFrame, flow: str, level: str, value: str = "leaf"):
    """Join one flow's per-region ``value`` onto the boundary geometry.

    Returns a GeoDataFrame with a ``value`` column, or ``None`` if geometry is
    unavailable.
    """

    boundaries = load_boundaries(level)
    if boundaries is None:
        return None
    join = _GEOMETRY[level]["join"]
    vals = df[(df["flow"] == flow) & (df["level"] == level)][["region_id", value]].copy()
    gdf = boundaries.copy()
    gdf["_key"] = gdf[join].astype("string")
    merged = gdf.merge(vals.rename(columns={"region_id": "_key", value: "value"}), on="_key", how="left")
    return merged


def polygon_areas_km2(gdf, equal_area_crs: str = "EPSG:6933") -> pd.Series:
    """Polygon areas in km^2 computed in an equal-area projection."""

    return gdf.to_crs(equal_area_crs).area / 1e6


def plot_choropleth_panels(df: pd.DataFrame, flow: str, value: str = "leaf"):
    """Three-panel country/subcountry/ecoregion choropleth on a shared log scale.

    Returns ``(fig, axes)`` or ``None`` when no level has geometry available.
    Renders only the panels whose geometry is present.
    """

    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors

    geodata = {level: choropleth_geodata(df, flow, level, value) for level in LEVELS}
    available = {k: v for k, v in geodata.items() if v is not None}
    if not available:
        warnings.warn("No boundary geometry available for any level; skipping maps.")
        return None

    try:
        from sbtn_leaf.map_plotting import _get_world_map

        basemap = _get_world_map()
    except (FileNotFoundError, OSError, Exception):  # pragma: no cover - env dependent
        basemap = None

    allv = pd.concat([g["value"] for g in available.values()]).dropna()
    allv = allv[allv > 0]
    norm = mcolors.LogNorm(vmin=allv.min(), vmax=allv.max()) if len(allv) else None

    fig, axes = plt.subplots(1, len(available), figsize=(7 * len(available), 6))
    if len(available) == 1:
        axes = [axes]
    for ax, (level, gdf) in zip(axes, available.items()):
        if basemap is not None:
            basemap.plot(ax=ax, color="lightgrey", edgecolor="white", linewidth=0.3)
        gdf.plot(ax=ax, column="value", cmap="viridis", norm=norm, legend=False, missing_kwds={"color": "whitesmoke"})
        ax.set_title(f"{level}")
        ax.set_axis_off()
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=axes, shrink=0.6, label=f"{value} ({UNIT})")
    fig.suptitle(f"{flow_label(flow)} — soil erosion by aggregation level", fontsize=14)
    return fig, axes


def plot_area_vs_value(df: pd.DataFrame, level: str, flow: str, value: str = "leaf"):
    """Scatter of polygon area vs ``value`` (with Spearman rho), or ``None``."""

    import matplotlib.pyplot as plt
    from scipy import stats

    gdf = choropleth_geodata(df, flow, level, value)
    if gdf is None:
        return None
    gdf = gdf[gdf["value"].notna()].copy()
    gdf["area_km2"] = polygon_areas_km2(gdf)
    x = gdf["area_km2"].to_numpy()
    y = gdf["value"].to_numpy()
    mask = (x > 0) & (y > 0)
    x, y = x[mask], y[mask]
    rho, p = stats.spearmanr(x, y) if len(x) > 2 else (float("nan"), float("nan"))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, s=10, alpha=0.4, color=_level_colors()[1])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Polygon area (km², log)")
    ax.set_ylabel(f"{value} ({UNIT}, log)")
    ax.set_title(f"{flow_label(flow)} @ {level}: area vs {value}\nSpearman ρ={rho:.2f} (p={p:.1e}, n={len(x)})")
    ax.grid(True, linestyle="--", alpha=0.4)
    return fig, ax


# --------------------------------------------------------------------------- #
# One-shot pipeline
# --------------------------------------------------------------------------- #


def default_output_dir() -> Path:
    return leaf_path("soil_erosion", "aggregation_analysis")


def run_full_analysis(
    outdir: Optional[Path] = None,
    flows: Optional[Iterable[str]] = None,
    *,
    make_figures: bool = True,
) -> Dict[str, object]:
    """Run the whole comparison, writing tables (+ optionally figures) to ``outdir``.

    Produces ``tables/`` (coverage, full and focal summaries, hierarchy variance),
    ``figures/`` (per focal flow + ranking + heatmap) and ``README.md``.  Map
    figures are attempted and silently skipped when geometry is unavailable.  Only
    the tables and ``README.md`` are version-controlled (see ``LEAFs/.gitignore``);
    the figures are regenerable here and embedded in the analysis notebook.
    Returns a dict with the key DataFrames for interactive use.
    """

    outdir = Path(outdir) if outdir is not None else default_output_dir()
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    if make_figures:
        figures_dir.mkdir(parents=True, exist_ok=True)

    flows = list(flows) if flows is not None else list(FOCAL_FLOWS)

    raw = load_harmonized(drop_na=False)
    df = raw[raw["leaf"].notna()].reset_index(drop=True)

    coverage = coverage_table(raw)
    flow_cov = validate_flow_coverage(raw)
    full_summary = summary_table(df)
    focal_summary = summary_table(df, flows=flows)
    variance = hierarchy_variance_table(df, flows=flows)

    coverage.to_csv(tables_dir / "coverage_by_level.csv", index=False)
    flow_cov.to_csv(tables_dir / "flow_coverage_matrix.csv")
    full_summary.to_csv(tables_dir / "summary_all_flows.csv", index=False)
    focal_summary.to_csv(tables_dir / "summary_focal_flows.csv", index=False)
    variance.to_csv(tables_dir / "hierarchy_variance.csv", index=False)

    n_maps = 0
    if make_figures:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for flow in flows:
            for name, fn in (
                ("box_leaf", plot_cross_level_box),
                ("box_within_std", plot_within_std_box),
                ("dist_overlay", plot_distribution_overlay),
            ):
                fig, _ = fn(df, flow)
                fig.tight_layout()
                fig.savefig(figures_dir / f"{name}__{_safe(flow)}.png", dpi=150, bbox_inches="tight")
                plt.close(fig)

        fig, _ = plot_commodity_ranking(df, level="ecoregion", flows=flows)
        fig.tight_layout()
        fig.savefig(figures_dir / "commodity_ranking_ecoregion.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        fig, _ = plot_sensitivity_heatmap(df, flows=flows)
        fig.tight_layout()
        fig.savefig(figures_dir / "sensitivity_heatmap_mean.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        for flow in flows[:3]:
            fig, _ = plot_biome_box(df, flow)
            fig.tight_layout()
            fig.savefig(figures_dir / f"biome_box__{_safe(flow)}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

        # Maps (skipped gracefully when geometry is absent).
        for flow in flows[:3]:
            result = plot_choropleth_panels(df, flow)
            if result is not None:
                fig, _ = result
                fig.savefig(figures_dir / f"map_panels__{_safe(flow)}.png", dpi=150, bbox_inches="tight")
                plt.close(fig)
                n_maps += 1

    # Written as README.md so it is tracked under LEAFs/ (see LEAFs/.gitignore,
    # which keeps data files + README.md but ignores loose images / other .md).
    _write_findings(outdir / "README.md", coverage, focal_summary, variance, n_maps)

    return {
        "df": df,
        "coverage": coverage,
        "flow_coverage": flow_cov,
        "summary_all": full_summary,
        "summary_focal": focal_summary,
        "variance": variance,
        "outdir": outdir,
        "n_maps_rendered": n_maps,
    }


def _safe(flow: str) -> str:
    return flow.replace("|", "_").replace(" ", "_")


def _fmt_cell(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def _df_to_md(df: pd.DataFrame) -> str:
    """Render a DataFrame as a GitHub-flavoured markdown table (no extra deps)."""

    cols = [str(c) for c in df.columns]
    head = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = [
        "| " + " | ".join(_fmt_cell(v) for v in row) + " |"
        for row in df.itertuples(index=False, name=None)
    ]
    return "\n".join([head, sep, *body])


def _write_findings(path: Path, coverage, focal_summary, variance, n_maps) -> None:
    lines = ["# Soil-erosion aggregation comparison — findings\n"]
    lines.append(f"Unit: **{UNIT}**. Levels (coarse→fine): {', '.join(LEVELS)}.\n")
    lines.append("## Coverage by level\n")
    lines.append(_df_to_md(coverage.round(1)))
    lines.append("\n\n## Focal-commodity summary (mean / median / between-region SD)\n")
    show = focal_summary[
        ["flow_label", "level", "n_regions", "mean_leaf", "median_leaf", "std_leaf", "mean_within_std", "max_leaf"]
    ].round(2)
    lines.append(_df_to_md(show))
    lines.append("\n\n## Nested variance (aggregation information loss)\n")
    lines.append(
        "`subcty_within_country_frac` = share of admin-1 erosion variance hidden inside a single "
        "national LEAF. `eco_between_biome_eta2` / `eco_between_realm_eta2` = share of ecoregion "
        "erosion variance explained by ecological grouping.\n\n"
    )
    lines.append(_df_to_md(variance.round(3)))
    lines.append(f"\n\nMaps rendered this run: **{n_maps}** (0 ⇒ boundary geometry not present; map code is drop-in).\n")
    lines.append(
        "\nFigures live in `figures/` (regenerated by "
        "`sbtn_leaf.se_aggregation_analysis.run_full_analysis()` and embedded in "
        "`examples/SoilErosion_Aggregation_Comparison.ipynb`); they are not "
        "version-controlled per `LEAFs/.gitignore`.\n"
    )
    path.write_text("\n".join(lines))
