# Claude analysis — spatial aggregation of SBTN-Land LEAFs

Notebooks and scripts analysing how the **aggregation level** (country →
subcountry → ecoregion) shapes the published LEAF values, and how much
ecological grouping (biome / realm) explains. The reusable code lives in
[`src/sbtn_leaf/claude_analysis/`](../../src/sbtn_leaf/claude_analysis).

## Contents

| file | what it is |
|---|---|
| [`Ecoregion_Aggregation_CrossIndicator.ipynb`](Ecoregion_Aggregation_CrossIndicator.ipynb) | **Main deliverable.** In-depth comparison of how ecoregions represent **SOC**, **soil-erosion** and **terrestrial-acidification** averages differently than sub-country / country units. Extends sections 6–7 of the soil-erosion notebook across all three indicators. |
| [`SoilErosion_Aggregation_Comparison.ipynb`](SoilErosion_Aggregation_Comparison.ipynb) | The original single-indicator (soil-erosion) aggregation comparison, relocated here. |
| `_build_notebook.py` | Build script that regenerates the cross-indicator notebook's cells. |
| `outputs/` | Generated tables (`tables/*.csv`), figures (`figures/*.png`) and the findings `README.md`. |

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

# or re-run the narrative notebook end-to-end
jupyter nbconvert --to notebook --execute --inplace \
    paper/claude_analysis/Ecoregion_Aggregation_CrossIndicator.ipynb
```

Maps are intentionally omitted: the analysis is purely statistical and needs no
boundary geometry (the DVC-tracked shapefiles are not required).
