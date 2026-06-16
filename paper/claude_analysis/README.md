# Claude analysis — spatial aggregation of SBTN-Land LEAFs

Notebooks and scripts analysing how the **aggregation level** (country →
subcountry → ecoregion) shapes the published LEAF values, and how much
ecological grouping (biome / realm) explains. The reusable code lives in
[`src/sbtn_leaf/claude_analysis/`](../../src/sbtn_leaf/claude_analysis).

## Contents

| file | what it is |
|---|---|
| [`Ecoregion_Aggregation_CrossIndicator.ipynb`](Ecoregion_Aggregation_CrossIndicator.ipynb) | **Main deliverable.** In-depth comparison of how ecoregions represent **SOC**, **soil-erosion** and **terrestrial-acidification** averages differently than sub-country / country units. Extends sections 6–7 of the soil-erosion notebook across all three indicators. |
| [`Manuscript_Support_Figures.ipynb`](Manuscript_Support_Figures.ipynb) | Extra figures + statistics backing specific claims in the LEAF manuscript (biome significance tests, within-region SD by level, SOC↔erosion alignment). |
| [`MultiIndicator_PracticeChange.ipynb`](MultiIndicator_PracticeChange.ipynb) | Practice-switch co-benefits (Fig. 13): SOC gained **and** erosion avoided when switching to reduced tillage + residue retention — maps of where to focus, benefit distributions, and residue-vs-tillage attribution. |
| [`SoilErosion_Aggregation_Comparison.ipynb`](SoilErosion_Aggregation_Comparison.ipynb) | The original single-indicator (soil-erosion) aggregation comparison, relocated here. |
| `_build_*.py` | Build scripts that regenerate each notebook's cells. |
| `outputs/` | Generated tables (`tables/*.csv`), figures (`figures/*.png`) and findings `README.md`; manuscript-support in `outputs/manuscript_support/`, practice-change (incl. maps) in `outputs/practice_change/`. |

## What the cross-indicator analysis answers

> *Do ecoregions represent SOC, soil-erosion and acidification averages
> differently than sub-country or country units — and for which indicator do
> ecological boundaries matter most?*

Headline result (median across each indicator's focal flows):

| indicator | biome η² (ecoregion) | realm η² (ecoregion) | within-country variance hidden |
|---|---|---|---|
| SOC stock | 0.17 | 0.17 | 0.32 |
| Soil erosion | 0.26 | 0.16 | **0.39** |
| Acidification | **0.32** | **0.38** | 0.25 |

* **Acidification** is the most "ecological": biome/realm explain the most
  ecoregion-level variance (realm η² for NOₓ reaches 0.52).
* **Soil erosion** is the indicator a national LEAF most badly misrepresents
  (39 % of its variance is hidden inside countries).
* **SOC** is the most locally driven (biome explains the least).
* For all three, **ecoregion averages run higher than national ones** (5–16 %):
  country polygons average commodity hotspots together with surrounding land.

## Reproduce

From the repository root (with the project installed, e.g. `pip install -e .`):

```bash
# regenerate every table, figure and the findings README under outputs/
python -m sbtn_leaf.claude_analysis.run_cross_indicator
python -m sbtn_leaf.claude_analysis.run_manuscript_support
python -m sbtn_leaf.claude_analysis.run_practice_change   # maps fetch geometry once, then cache

# or re-run the narrative notebooks end-to-end
jupyter nbconvert --to notebook --execute --inplace \
    paper/claude_analysis/Ecoregion_Aggregation_CrossIndicator.ipynb
jupyter nbconvert --to notebook --execute --inplace \
    paper/claude_analysis/Manuscript_Support_Figures.ipynb
jupyter nbconvert --to notebook --execute --inplace \
    paper/claude_analysis/MultiIndicator_PracticeChange.ipynb
```

### Practice-switch co-benefits (`outputs/practice_change/`)

Switching the **same commodity** from conventional-till + residues-removed to reduced-till +
residues-left (manuscript Fig. 13):

| finding | statistic |
|---|---|
| extent of benefit (wheat) | median **+8.5 t SOC/ha (+25 %)** and **−13 t soil/ha/yr (−77 %)**; **96 % win-win** |
| relative erosion cut is location-independent | constant per commodity (RUSLE C-factor): wheat −77 %, maize −69 %, soy/cotton −65 % |
| the two benefits don't co-locate | SOC-gain ↔ erosion-avoided Spearman ρ ≈ **−0.18…−0.29** → use the combined **priority** score |
| what drives each | SOC ← **residue retention**; erosion ← **reduced tillage** (complementary levers) |

Maps (country via Natural Earth; ecoregion via RESOLVE Ecoregions-2017 by `ECO_ID`) are fetched once
and cached in `~/.cache/sbtn_leaf_geo`; all map code skips gracefully when geometry is unavailable.

### Manuscript-support statistics (`outputs/manuscript_support/`)

| manuscript claim | supporting statistic |
|---|---|
| *“ecoregion’s biomes lead to LEAFs that are significantly different”* | Kruskal–Wallis p < 1e-40 for SOC, erosion & acidification; biome η² = 0.17 / 0.26 / 0.32 |
| *“sub-country … smaller standard deviation”* | within-region SD = 0.71× / 0.81× / 0.51× the country value at sub-country |
| multi-indicator *“aligned SOC and soil erosion”* | ecoregion SOC↔erosion Spearman ρ = 0.23–0.60 (all p < 0.001) |

The aggregation and manuscript-support analyses are purely statistical and need no
geometry. Only the practice-change notebook draws maps, and it fetches/caches its
own geometry (above), so no DVC shapefiles are required for any of this.
