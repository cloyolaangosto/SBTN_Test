"""Generate the cross-indicator ecoregion-aggregation narrative notebook.

This is a build helper (not part of the analysis): it assembles the notebook
cells and writes the .ipynb, which is then executed with ``jupyter nbconvert``.
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def md(*lines):
    return {"cell_type": "markdown", "metadata": {}, "source": _src(lines)}


def code(*lines):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": _src(lines)}


def _src(lines):
    text = "\n".join(lines)
    parts = text.split("\n")
    return [p + "\n" for p in parts[:-1]] + [parts[-1]]


cells = []

cells.append(md(
    "# Ecoregions vs political units: how SOC, soil erosion & acidification LEAFs change with polygon size",
    "",
    "The SBTN-Land LEAFs are published at three spatial **aggregation levels** that clip the",
    "same underlying rasters to three different polygon sets:",
    "",
    "| level | polygons | count |",
    "|---|---|---|",
    "| **country** | FAO ADM0 | 276 |",
    "| **subcountry** | FAO ADM1 | 3 422 |",
    "| **ecoregion** | WWF 2017 ecoregions | ~830–847 |",
    "",
    "This notebook extends sections **6 (ecoregion significance — biomes & nested variance)**",
    "and **7 (sensitivity heatmap)** of [`SoilErosion_Aggregation_Comparison.ipynb`](SoilErosion_Aggregation_Comparison.ipynb)",
    "from one indicator to **three**, to answer one question directly:",
    "",
    "> **Do ecoregions represent SOC, soil-erosion and terrestrial-acidification averages",
    "> *differently* than sub-country or country units — and for which indicator do",
    "> ecological boundaries matter most?**",
    "",
    "| indicator | `leaf` value | unit | direction |",
    "|---|---|---|---|",
    "| **SOC** | 2030 soil-organic-carbon stock | t SOC/ha | higher = more carbon stored |",
    "| **Soil erosion** | RUSLE soil loss | t soil/ha/yr | higher = worse |",
    "| **Terrestrial acidification** | accumulated-exceedance CF (NOₓ, NH₃, SO₂) | kg SO₂-eq./kg | higher = worse |",
    "",
    "All logic lives in `sbtn_leaf.claude_analysis` (engine: `indicator_aggregation`,",
    "configs: `indicators`, comparison: `cross_indicator`); this notebook is the narrative",
    "front-end. Tables and figures are also written to `paper/claude_analysis/outputs/` by",
    "`xi.run_cross_indicator_analysis()` (last cell).",
    "",
    "> **A note on the three indicators.** SOC and soil erosion are *land-use* LEAFs — one",
    "> value per region **×** commodity/management — so they share the canonical flow key",
    "> (`Commodity|water|residue|tillage`). Acidification is a *characterisation factor*",
    "> defined region-wide for three acidifying gases, so its 'flows' are the gases. The",
    "> aggregation machinery is identical; only the flow taxonomy differs.",
))

cells.append(code(
    "import warnings; warnings.filterwarnings('ignore')",
    "%matplotlib inline",
    "import pandas as pd, matplotlib.pyplot as plt",
    "pd.set_option('display.width', 200, 'display.max_columns', 30)",
    "",
    "from sbtn_leaf.claude_analysis import indicator_aggregation as eng",
    "from sbtn_leaf.claude_analysis import cross_indicator as xi",
    "from sbtn_leaf.claude_analysis.indicators import SOC, SOIL_EROSION, ACIDIFICATION, INDICATORS",
    "",
    "# Harmonise every indicator onto the shared long schema. `raw` keeps the 'no data'",
    "# cells so we can measure coverage; `data` drops them for the statistics.",
    "raw  = {k: cfg.load_harmonized(drop_na=False) for k, cfg in INDICATORS.items()}",
    "data = {k: r[r['leaf'].notna()].reset_index(drop=True) for k, r in raw.items()}",
    "",
    "# A representative flow per indicator for the per-flow figures.",
    "reps = {'soc': 'Wheat|rf|roff|ct', 'soil_erosion': 'Wheat|rf|roff|ct', 'acidification': 'acid_so2'}",
    "data['soc'].head()",
))

cells.append(md(
    "## 1. Harmonisation & coverage",
    "",
    "Each indicator's source tables use different `flow_name` conventions and column layouts",
    "(SOC `v1.0` tables carry a `variable` column; acidification tables carry a `metric`",
    "column; the SOC ecoregion table even omits biome/realm, which we backfill). The engine",
    "maps every spelling onto one canonical key and one `leaf`/`leaf_median`/`leaf_std` row",
    "per region × flow, so the three indicators become directly comparable.",
))

cells.append(code(
    "coverage = pd.concat(",
    "    [eng.coverage_table(raw[k]).assign(indicator=cfg.name) for k, cfg in INDICATORS.items()],",
    "    ignore_index=True,",
    ")",
    "coverage[['indicator', 'level', 'n_regions', 'n_flows', 'coverage_pct']].round(1)",
))

cells.append(md(
    "Coverage (share of region × flow cells with a modelled value) **rises as polygons get",
    "finer** for the two land-use indicators — SOC 48 → 51 → 63 %, soil erosion 51 → 58 → 65 %",
    "— so ecologically defined polygons resolve the most land where commodities actually grow.",
    "Acidification is ~100 % at every level because its characterisation factor is defined",
    "for every region regardless of land use.",
))

cells.append(md(
    "## 2. Cross-level summary statistics (mean / median / SD by level)",
    "",
    "The core mean/median/SD comparison, for each indicator's focal flows. `mean_leaf`,",
    "`median_leaf` and `std_leaf` are the **between-region** statistics; `mean_within_std`",
    "is the mean per-region `leaf_std` (typical **within-region** heterogeneity); `cv` is the",
    "between-region coefficient of variation.",
))

cells.append(code(
    "combined = xi.combined_summary_table()",
    "combined[['indicator_name', 'flow_label', 'level', 'n_regions', 'mean_leaf',",
    "          'median_leaf', 'std_leaf', 'mean_within_std', 'cv']].round(2)",
))

cells.append(code(
    "# Cross-level boxplots for one representative flow per indicator",
    "for k, flow in reps.items():",
    "    fig, ax = eng.plot_cross_level_box(INDICATORS[k], data[k], flow)",
    "    display(fig); plt.close(fig)",
))

cells.append(md(
    "The boxes show how the **published value's spread depends on polygon size**. For the",
    "land-use indicators the country boxes are compressed toward the middle, while finer",
    "polygons expose a longer tail — the hotspots national averages hide.",
))

cells.append(md(
    "## 3. Ecoregion significance — biomes and nested variance",
    "",
    "*(the section-6 analogue, now across three indicators)*",
    "",
    "Two complementary views:",
    "",
    "* **Biome breakdown** — ecoregion `leaf` for the representative flow, split by WWF biome.",
    "  Clear separation between biomes means the ecological grouping carries real signal.",
    "* **Nested variance decomposition** —",
    "  * `subcty_within_country_frac`: the share of admin-1 variance a single **national** LEAF",
    "    hides inside its sub-units (the cost of political coarsening).",
    "  * `eco_between_biome_eta2` / `eco_between_realm_eta2`: the share of ecoregion variance",
    "    **explained** by biome / realm (ecoregion significance).",
))

cells.append(code(
    "for k, flow in reps.items():",
    "    fig, ax = eng.plot_biome_box(INDICATORS[k], data[k], flow)",
    "    display(fig); plt.close(fig)",
))

cells.append(code(
    "# Per-indicator nested variance across each indicator's focal flows",
    "variance = pd.concat(",
    "    [eng.hierarchy_variance_table(cfg, data[k], flows=list(cfg.focal_flows))",
    "     for k, cfg in INDICATORS.items()],",
    "    ignore_index=True,",
    ")",
    "variance.round(3)",
))

cells.append(code(
    "fig, ax = xi.plot_significance_comparison(); display(fig); plt.close(fig)",
    "xi.significance_table().round(3)",
))

cells.append(md(
    "**This is the headline result.** Summarised across each indicator's focal flows (median):",
    "",
    "| indicator | biome η² | realm η² | within-country hidden |",
    "|---|---|---|---|",
    "| SOC stock | 0.17 | 0.17 | 0.32 |",
    "| Soil erosion | 0.26 | 0.16 | 0.39 |",
    "| Acidification | **0.32** | **0.38** | 0.25 |",
    "",
    "* **Acidification is the most 'ecological' indicator.** Biome and realm explain the",
    "  *largest* share of ecoregion variance (realm η² for NOₓ reaches **0.52**), because the",
    "  acidification CF is driven by soil sensitivity and atmospheric fate, which track",
    "  climatic/ecological gradients closely. For acidification, an ecoregion average is",
    "  genuinely more representative than a national one.",
    "* **Soil erosion hides the most inside national borders** (within-country fraction 0.39):",
    "  a single country LEAF averages away the largest slice of real sub-national erosion",
    "  contrast — so country-level erosion LEAFs are the least representative of the three.",
    "* **SOC is the most locally driven.** Biome explains the least (0.17): SOC stock varies",
    "  more with local soil type than with broad ecological zones, so neither political nor",
    "  ecological polygons capture it cleanly — though oil-palm SOC still hides **0.51** of its",
    "  variance within countries (perennials in heterogeneous tropics).",
))

cells.append(md(
    "## 4. Sensitivity heatmaps",
    "",
    "*(the section-7 analogue)* — a compact view of mean `leaf` across flow × level (log",
    "colour scale) for each indicator. Reading left→right shows how the central estimate of",
    "every flow shifts as polygons shrink country → subcountry → ecoregion.",
))

cells.append(code(
    "for k, cfg in INDICATORS.items():",
    "    fig, ax = eng.plot_sensitivity_heatmap(cfg, data[k])",
    "    display(fig); plt.close(fig)",
))

cells.append(md(
    "Within an indicator the **ranking of commodities is stable across levels** (a column",
    "re-colours but the row order barely changes), confirming the level choice reframes the",
    "*magnitude* of the average more than the *ordering* of which land uses are most",
    "carbon-rich / erosive / acidifying.",
))

cells.append(md(
    "## 5. How ecoregions reframe the average",
    "",
    "The most direct answer to the question. Because the three indicators have different",
    "units, each focal flow's level means are normalised to its **own country mean",
    "(country = 1.0)** and then summarised (median) across focal flows.",
))

cells.append(code(
    "fig, ax = xi.plot_level_reframing(); display(fig); plt.close(fig)",
    "xi.level_reframing_table().round(3)",
))

cells.append(code(
    "fig, ax = xi.plot_dispersion_by_level(); display(fig); plt.close(fig)",
    "xi.dispersion_table().round(3)",
))

cells.append(md(
    "* **Ecoregion averages run higher than national ones** for all three indicators (mean",
    "  relative to country: SOC 1.07, erosion 1.05, acidification 1.16). National averaging",
    "  dilutes the high-value regions; ecoregions restore them. For acidification the gap",
    "  opens already at the subcountry level (1.19) — political aggregation understates the",
    "  CF by ~15–20 %.",
    "* **Sub-country is not simply 'between' country and ecoregion.** For SOC and soil erosion",
    "  the subcountry mean dips *below* the country mean (0.98 / 0.89) before the ecoregion",
    "  mean rises above it — admin-1 polygons fragment commodities into their lower-value",
    "  growing margins, while ecoregions re-aggregate along the ecological gradient.",
    "* **Dispersion behaves differently per indicator.** SOC's between-region CV is flat",
    "  (~0.39) — SOC stock is similarly variable at every scale. Soil erosion peaks at",
    "  subcountry (CV 1.50): admin-1 polygons slice through erosion hotspots. Acidification",
    "  dispersion dips at subcountry, consistent with its strong realm-level structure.",
))

cells.append(md(
    "## 6. Findings",
    "",
    "1. **Ecoregions are the most informative level for acidification, the most *necessary*",
    "   for soil erosion, and the least decisive for SOC.** Ecological grouping explains a",
    "   third of acidification variance but only a sixth of SOC variance; soil erosion is the",
    "   indicator a national LEAF most badly misrepresents (39 % of variance hidden",
    "   within-country).",
    "2. **Ecoregion averages are systematically higher than national averages** (5–16 %),",
    "   because country polygons average commodity hotspots together with the land around",
    "   them. Using a country LEAF therefore tends to *understate* SOC opportunity, erosion",
    "   risk and acidification impact alike.",
    "3. **Coverage favours ecoregions too**: for the land-use indicators ecoregions resolve",
    "   ~13–14 percentage points more of the region × commodity grid than countries do.",
    "4. **Practical guidance** (consistent with `LEAFs/README.md`'s 'country as last resort'):",
    "   prefer the ecoregion LEAF where an ecological signal dominates (acidification, and",
    "   biome-structured erosion); fall back to subcountry where commodity-specific local",
    "   detail matters more than the ecological gradient (much of SOC).",
))

cells.append(md(
    "## 7. Reproduce all artifacts",
    "",
    "Regenerates every table and figure under `paper/claude_analysis/outputs/` and rewrites",
    "the findings `README.md` (also available from the CLI:",
    "`python -m sbtn_leaf.claude_analysis.run_cross_indicator`).",
))

cells.append(code(
    "result = xi.run_cross_indicator_analysis()",
    "print('Wrote artifacts to:', result['outdir'])",
    "result['significance'].round(3)",
))

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out = HERE / "Ecoregion_Aggregation_CrossIndicator.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n")
print("wrote", out)
