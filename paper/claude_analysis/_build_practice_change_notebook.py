"""Generate the multi-indicator practice-change notebook.

Build helper (not part of the analysis): assembles the cells and writes the
.ipynb, which is then executed with ``jupyter nbconvert``.
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
    "# Multi-indicator analysis: practice-switch co-benefits for SOC **and** soil erosion",
    "",
    "Expands the LEAF manuscript's **Multi-indicator analysis** (p. 17, Fig. 13) and **Conclusions**:",
    "",
    "> *“By providing aligned SOC and soil erosion results for each land use, it is possible to identify",
    "> where regenerative agriculture practices, like reduced tillage and residue management, can",
    "> provide the most benefits for both indicators simultaneously.”*",
    "",
    "For the **same commodity** we switch from the baseline practice (**conventional tillage + residues",
    "removed**) to the regenerative practice (**reduced tillage + residues left**) — exactly the wheat",
    "switch the manuscript describes — and quantify, per region:",
    "",
    "| metric | meaning |",
    "|---|---|",
    "| `d_soc` | SOC gained, t SOC/ha (and `d_soc_pct`) |",
    "| `d_se_red` | soil erosion avoided, t soil/ha/yr (and `d_se_red_pct`) |",
    "| `win_win` | both improve |",
    "| `priority` | mean of the within-commodity percentile ranks of the two **absolute** benefits — *where to focus* |",
    "",
    "Cereals (Wheat, Maize) switch residues *and* tillage; Soybeans & Cotton switch tillage only",
    "(residue management is not an option). Logic lives in",
    "`sbtn_leaf.claude_analysis.practice_change`; maps use `…claude_analysis.geo` (country = Natural",
    "Earth, ecoregion = RESOLVE Ecoregions-2017 by `ECO_ID`, fetched once and cached).",
    "",
    "> **Why absolute benefits drive the maps.** Because RUSLE erosion scales with the multiplicative",
    "> C-factor, the *percentage* erosion reduction of a given switch is spatially **constant** per",
    "> commodity (≈77 % for wheat); only the *absolute* tonnage avoided varies by location. SOC gain",
    "> varies in both magnitude and sign. So the spatial \"where to focus\" ranking uses absolute benefits.",
))

cells.append(code(
    "import warnings; warnings.filterwarnings('ignore')",
    "%matplotlib inline",
    "import pandas as pd, matplotlib.pyplot as plt",
    "pd.set_option('display.width', 200, 'display.max_columns', 30)",
    "",
    "from sbtn_leaf.claude_analysis import practice_change as pc",
    "",
    "summary = pc.cobenefit_summary()",
    "summary.round(3)",
))

# ---- Section 1 ----
cells.append(md(
    "## 1. The extent of the benefits",
    "",
    "How large is the SOC gain and the erosion reduction from the switch, per commodity, summarised",
    "across the ecoregions where each crop grows.",
))

cells.append(code(
    "fig, axes = pc.plot_benefit_distributions(); display(fig); plt.close(fig)",
))

cells.append(md(
    "* **Cereals reap the largest SOC gain** because the switch also retains residues: the wheat",
    "  switch adds a median **+8.5 t SOC/ha (+25 %)** and avoids **~13 t soil/ha/yr (−77 %)**; maize",
    "  **+9.0 t SOC/ha** and **~31 t soil/ha/yr (−69 %)**.",
    "* **Tillage-only crops still cut erosion sharply** but gain little carbon: soybeans / cotton avoid",
    "  **~23–26 t soil/ha/yr (−65 %)** for only **+0.6–0.7 t SOC/ha**.",
    "* The switch is **win-win in 96–98 % of cereal regions** and ~72 % of soybean/cotton regions —",
    "  i.e. it essentially never makes either indicator worse.",
))

# ---- Section 2 ----
cells.append(md(
    "## 2. Where to focus — co-benefit priority",
    "",
    "> *“… highlight which location might lead to larger benefits for SOC and soil erosion at the same",
    "> time.”* — Conclusions",
    "",
    "First, do the two benefits co-locate? The scatter pairs each ecoregion's SOC gain with its erosion",
    "avoided (coloured by realm).",
))

cells.append(code(
    "fig, ax = pc.plot_cobenefit_scatter('Wheat', 'ecoregion'); display(fig); plt.close(fig)",
))

cells.append(md(
    "The Spearman correlation is **mildly negative** (ρ ≈ −0.18 for wheat, −0.23 maize): the regions",
    "with the *largest* erosion avoided (steep, wet, highly erosive tropics) are **not** generally the",
    "regions with the *largest* SOC gain. So there is a real spatial trade-off — which is exactly why a",
    "**combined priority score** (high only when *both* benefits rank high) is needed to locate the",
    "balanced sweet spots, rather than optimising either indicator alone.",
))

cells.append(md(
    "### Maps — the Fig. 13 analogue",
    "",
    "Three panels (SOC gained · erosion avoided · co-benefit priority) for the wheat switch, at the",
    "ecoregion level (the paper's focus) and at country level. *(Geometry is fetched once and cached;",
    "if unavailable the maps skip gracefully.)*",
))

cells.append(code(
    "res = pc.plot_cobenefit_map_panels('Wheat', 'ecoregion')",
    "if res: fig, _ = res; display(fig); plt.close(fig)",
    "else: print('Geometry unavailable — maps skipped.')",
))

cells.append(code(
    "res = pc.plot_cobenefit_map_panels('Wheat', 'country')",
    "if res: fig, _ = res; display(fig); plt.close(fig)",
    "else: print('Geometry unavailable — maps skipped.')",
))

cells.append(md(
    "The **priority** panel is the decision layer: dark = ecoregions/countries where the wheat switch",
    "delivers a high rank for SOC gain *and* erosion avoided together. Consistent with the manuscript's",
    "qualitative call-outs (western Mexico, parts of China/Japan, south-east Africa), the priority",
    "hotspots cluster in the warm, erosion-prone, residue-responsive croplands.",
))

cells.append(code(
    "# Top ecoregions to focus on for the wheat switch",
    "pc.priority_regions('Wheat', 'ecoregion', 12).round(2)",
))

cells.append(code(
    "# Geometry-free spatial view: median benefit by biogeographic realm",
    "fig, ax = pc.plot_cobenefit_by_realm('Wheat', 'ecoregion'); display(fig); plt.close(fig)",
))

# ---- Section 3 ----
cells.append(md(
    "## 3. What drives the benefit — residue vs tillage",
    "",
    "> *“… residue management is the largest factor contributing to increased SOC stock … reduced",
    "> tillage being the main driver behind [erosion reduction] across all crops.”* — Results",
    "",
    "Decomposing the cereal switch into its two components (residue retention alone vs reduced tillage",
    "alone) tests that claim directly.",
))

cells.append(code(
    "pc.practice_attribution_table().round(3)",
))

cells.append(code(
    "fig, axes = pc.plot_attribution('Wheat', 'ecoregion'); display(fig); plt.close(fig)",
))

cells.append(md(
    "Confirmed for both cereals: **residue retention dominates the SOC gain** (wheat: +6.6 vs +0.5",
    "t SOC/ha for tillage), while **reduced tillage dominates the erosion benefit** (wheat: −10.7 vs",
    "−2.7 t soil/ha/yr for residue). The two practices are therefore **complementary** — each is the",
    "primary lever for a *different* indicator, so stacking them is what produces the win-win.",
))

cells.append(md(
    "## Findings",
    "",
    "1. **The switch is almost always a win-win** (96–98 % of cereal regions): it raises SOC and lowers",
    "   erosion at once, with cereal SOC gains (~+25 %) far exceeding tillage-only crops (~+2 %).",
    "2. **But the two benefits peak in different places** (Spearman ρ < 0), so a combined **priority**",
    "   score is required to target locations that score high on *both* — the operational answer to the",
    "   manuscript's *“where to focus.”*",
    "3. **Residue retention drives SOC, reduced tillage drives erosion control** — complementary levers,",
    "   quantified here, matching the manuscript's qualitative claim.",
    "4. The maps turn the aligned LEAFs into a deployable decision layer at both ecoregion and country",
    "   resolution.",
))

cells.append(md(
    "## Reproduce all artifacts",
    "",
    "Writes tables, figures and maps under `paper/claude_analysis/outputs/practice_change/` and rewrites",
    "the findings `README.md` (CLI: `python -m sbtn_leaf.claude_analysis.run_practice_change`).",
))

cells.append(code(
    "result = pc.run_practice_change_analysis()",
    "print('Wrote artifacts to:', result['outdir'], '| maps rendered:', result['n_maps'])",
    "result['summary'].round(3)",
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

out = HERE / "MultiIndicator_PracticeChange.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n")
print("wrote", out)
