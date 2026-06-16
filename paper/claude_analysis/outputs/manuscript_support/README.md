# Manuscript-support statistics — findings

Quantitative backing for the LEAF manuscript's *Regional averages* and *Conclusions* claims.


## 1. Biomes produce significantly different LEAFs

Per-indicator median biome / realm η² over focal flows, and the worst-case (largest) Kruskal–Wallis p across those flows. Supports *“ecoregion’s biomes lead to LEAFs that are significantly different.”*


| indicator | indicator_name | n_flows | eta2_biome | eta2_realm | kruskal_p_max |
| --- | --- | --- | --- | --- | --- |
| soc | SOC stock | 9 | 0.167 | 0.1709 | 0 |
| soil_erosion | Soil erosion | 9 | 0.2595 | 0.1559 | 0 |
| acidification | Acidification CF | 3 | 0.3231 | 0.378 | 0 |


## 2. Sub-country polygons have the smallest within-region SD

Within-region SD by level, normalised to country = 1.0 (median over focal flows). Supports *“sub-country … leads to smaller standard deviation.”*


| indicator_name | within_sd_rel_country | within_sd_rel_subcountry | within_sd_rel_ecoregion |
| --- | --- | --- | --- |
| SOC stock | 1 | 0.706 | 1.063 |
| Soil erosion | 1 | 0.814 | 1.059 |
| Acidification CF | 1 | 0.507 | 1.046 |


## 3. SOC and soil erosion are aligned across regions

Per-region Spearman correlation of SOC vs soil erosion at the ecoregion level, for the shared focal commodities. Supports the multi-indicator section.


| flow_label | n_regions | spearman_rho | p | signif |
| --- | --- | --- | --- | --- |
| Broadleaf-decid. forest (tropical) | 418 | 0.4934 | 0 | *** |
| Coffee (rainfed) | 358 | 0.5978 | 0 | *** |
| Cotton (rainfed) | 497 | 0.4196 | 0 | *** |
| Grassland | 694 | 0.3519 | 0 | *** |
| Maize (rainfed) | 675 | 0.2684 | 0 | *** |
| Oil palm (rainfed) | 241 | 0.5936 | 0 | *** |
| Soybeans (rainfed) | 590 | 0.2292 | 0 | *** |
| Sugarcane (rainfed) | 411 | 0.2444 | 0 | *** |
| Wheat (rainfed) | 619 | 0.3335 | 0 | *** |


Regenerate with `python -m sbtn_leaf.claude_analysis.run_manuscript_support` or `sbtn_leaf.claude_analysis.manuscript_support.run_manuscript_support()`.
