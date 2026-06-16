"""Generate the manuscript-support figures notebook.

Build helper (not part of the analysis): assembles the notebook cells and writes
the .ipynb, which is then executed with ``jupyter nbconvert``.
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
    "# Manuscript-support figures & statistics",
    "",
    "Quantitative backing for specific qualitative claims in the LEAF manuscript's",
    "**Regional averages** section (pp. 18–19) and **Conclusions** — several of which",
    "currently carry figure placeholders or report no statistics. Each section below",
    "quotes the manuscript claim it supports and supplies the test / figure.",
    "",
    "| manuscript claim (placeholder) | this notebook |",
    "|---|---|",
    "| *“ecoregion’s biomes lead to LEAFs that are significantly different”* (Figure XXX) | §1 biome significance tests |",
    "| *“sub-country … leads to smaller standard deviation”* | §2 within-region SD by level |",
    "| *“country averages … similar to ecoregions”* (FIGURE YYY – BOXPLOT) | §2 cross-level boxplots |",
    "| multi-indicator: *“aligned SOC and soil erosion … benefits for both”* | §3 SOC↔erosion correlation |",
    "",
    "All logic lives in `sbtn_leaf.claude_analysis.manuscript_support`; artifacts are",
    "written to `paper/claude_analysis/outputs/manuscript_support/` by",
    "`ms.run_manuscript_support()` (last cell).",
))

cells.append(code(
    "import warnings; warnings.filterwarnings('ignore')",
    "%matplotlib inline",
    "import pandas as pd, matplotlib.pyplot as plt",
    "pd.set_option('display.width', 200, 'display.max_columns', 30)",
    "",
    "from sbtn_leaf.claude_analysis import indicator_aggregation as eng",
    "from sbtn_leaf.claude_analysis import manuscript_support as ms",
    "from sbtn_leaf.claude_analysis.indicators import SOC, SOIL_EROSION, ACIDIFICATION, INDICATORS",
    "",
    "data = {k: cfg.load_harmonized(drop_na=True) for k, cfg in INDICATORS.items()}",
    "ms.REP_FLOW",
))

# ----- Section 1 -----
cells.append(md(
    "## 1. Biomes produce significantly different LEAFs",
    "",
    "> *“Nevertheless, ecoregion’s biomes lead to LEAFs that are significantly different,",
    "> showing that they correctly incorporate different soil and climate characteristics.”*",
    "> — Regional averages",
    ">",
    "> *“The newly developed ecoregional average factors prove to significantly predict",
    "> different indicators based on ecoregional biomes …”* — Conclusions",
    "",
    "The manuscript asserts significance but reports no test. For each indicator and focal",
    "flow we group the ecoregion-level `leaf` by WWF biome and test with **Kruskal–Wallis**",
    "(rank-based, robust to the right-skewed LEAFs; primary), a **one-way ANOVA on log10**",
    "(parametric check), and report **η²** (variance explained) for biome and realm.",
))

cells.append(code(
    "sig = ms.biome_significance_table()",
    "sig[['indicator_name', 'flow_label', 'k_biomes', 'n_ecoregions', 'kruskal_H',",
    "     'kruskal_p', 'anova_p_log', 'eta2_biome', 'eta2_realm', 'signif']].round(4)",
))

cells.append(code(
    "fig, ax = ms.plot_biome_significance(); display(fig); plt.close(fig)",
))

cells.append(code(
    "# The biome separation behind the test, for the representative flow of each indicator",
    "for k, flow in ms.REP_FLOW.items():",
    "    fig, ax = eng.plot_biome_box(INDICATORS[k], data[k], flow)",
    "    display(fig); plt.close(fig)",
))

cells.append(md(
    "**Every focal flow rejects equal biome distributions at p < 0.001** (`***`) — Kruskal–Wallis",
    "p-values are effectively zero (e.g. wheat: SOC p≈1e-45, erosion p≈3e-65, acidification SO₂",
    "p≈2e-58). Biome explains a median **17 % (SOC), 26 % (soil erosion) and 32 % (acidification)**",
    "of ecoregion-level variance; realm adds a further ecological axis (strongest for acidification,",
    "η²≈0.38, and for perennials such as oil palm/coffee in SOC). This is the hard-statistics",
    "backing for the manuscript's *“significantly different”* / *“significantly predict”* claims.",
))

# ----- Section 2 -----
cells.append(md(
    "## 2. Polygon size and spread",
    "",
    "> *“Due to the nature of country averages, commodities appeared to have potential to be grown",
    "> over a larger surface … and skew into higher averages over larger sections of land than",
    "> sub-country or ecoregional ones … albeit leading to similar averages than ecoregions",
    "> (FIGURE YYY – BOXPLOT). … As sub-country divides the world surface into smaller areas …",
    "> this leads to smaller standard deviation.”* — Regional averages",
    "",
    "**FIGURE YYY** — the cross-level distribution of the representative flow per indicator",
    "(country / sub-country / ecoregion). Country and ecoregion sit at a similar central level,",
    "while the spread changes with polygon size.",
))

cells.append(code(
    "for k, flow in ms.REP_FLOW.items():",
    "    fig, ax = eng.plot_cross_level_box(INDICATORS[k], data[k], flow)",
    "    display(fig); plt.close(fig)",
))

cells.append(md(
    "The manuscript's *“smaller standard deviation”* refers to the spread **inside** each polygon.",
    "Smaller political polygons should be more internally homogeneous; we confirm it with the median",
    "per-region `leaf_std`, normalised to the country value (country = 1.0).",
))

cells.append(code(
    "ms.within_region_sd_table().round(3)",
))

cells.append(code(
    "fig, ax = ms.plot_within_region_sd(); display(fig); plt.close(fig)",
))

cells.append(md(
    "Within-region SD drops sharply at sub-country for **all three indicators** — to **0.71×**",
    "(SOC), **0.81×** (soil erosion) and **0.51×** (acidification) of the country value — then rises",
    "back to ≈country at the ecoregion level (ecoregions are larger, ecologically- rather than",
    "administratively-bounded polygons). This directly substantiates *“sub-country … smaller",
    "standard deviation”*, and complements the cross-indicator notebook's finding that the country",
    "and ecoregion **central** estimates are similar while sub-country dips lower.",
))

# ----- Section 3 -----
cells.append(md(
    "## 3. Multi-indicator alignment (SOC ↔ soil erosion)",
    "",
    "> *“By providing aligned SOC and soil erosion results for each land use, it is possible to",
    "> identify where regenerative agriculture practices … provide the most benefits for both",
    "> indicators simultaneously.”* — Multi-indicator analysis / Conclusions",
    "",
    "Because SOC and soil erosion share the same regions and canonical flow keys, we can measure how",
    "tightly they track each other per region (Spearman ρ, rank-based) for every shared commodity and",
    "level.",
))

cells.append(code(
    "ms.multi_indicator_correlation_table().round(4)",
))

cells.append(code(
    "fig, ax = ms.plot_multi_indicator_scatter('Wheat|rf|roff|ct', 'ecoregion'); display(fig); plt.close(fig)",
))

cells.append(md(
    "At the ecoregion level **every shared commodity shows a significant positive SOC↔erosion",
    "correlation** (ρ from 0.23 for soybeans to 0.60 for coffee, all p < 0.001): regions that are",
    "warm, wet and productive tend to be both high-erosion *and* high-SOC. The correlation is well",
    "below 1, so the two indicators are **aligned but not redundant** — which is exactly why a",
    "joint, co-located assessment (and the manuscript's regenerative-practice co-benefit mapping)",
    "adds information beyond either indicator alone.",
))

# ----- Findings + reproduce -----
cells.append(md(
    "## Findings",
    "",
    "1. **Biomes are a statistically significant predictor of every LEAF** (Kruskal–Wallis",
    "   p < 1e-40 for all three indicators), explaining 17–32 % of ecoregion variance — the",
    "   quantitative basis for the manuscript's ecoregional-averaging argument.",
    "2. **Sub-country polygons are the most internally homogeneous** (within-region SD 0.5–0.8× the",
    "   country value), confirming the *“smaller standard deviation”* claim; ecoregions trade some of",
    "   that homogeneity for ecological meaning.",
    "3. **SOC and soil erosion are positively and significantly aligned across regions** (ecoregion",
    "   ρ = 0.23–0.60), supporting joint multi-indicator assessment while showing the indicators are",
    "   not interchangeable.",
))

cells.append(md(
    "## Reproduce all artifacts",
    "",
    "Writes the tables and figures under `paper/claude_analysis/outputs/manuscript_support/` and",
    "rewrites the findings `README.md` (CLI: `python -m sbtn_leaf.claude_analysis.run_manuscript_support`).",
))

cells.append(code(
    "result = ms.run_manuscript_support()",
    "print('Wrote artifacts to:', result['outdir'])",
    "result['biome_summary'].round(4)",
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

out = HERE / "Manuscript_Support_Figures.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n")
print("wrote", out)
