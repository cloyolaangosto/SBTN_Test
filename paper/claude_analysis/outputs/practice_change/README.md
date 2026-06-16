# Multi-indicator practice-change co-benefits — findings

Switching the **same commodity** from the baseline (conventional tillage + residues removed) to the regenerative practice (reduced tillage + residues left) — SOC gained and soil erosion avoided per region. Expands the manuscript's multi-indicator section / Fig. 13.


## Extent of the benefits (per commodity, ecoregion level)

`med_d_soc` = median SOC gained (t SOC/ha); `med_d_se_red` = median erosion avoided (t soil/ha/yr); `d_se_red_pct` = relative erosion reduction (spatially constant per commodity — a property of the multiplicative RUSLE C-factor); `pct_win_win` = share of regions improving on **both**; `spearman_soc_se` = co-location of the two benefits.


| commodity | switch | n_regions | med_d_soc | med_d_soc_pct | med_d_se_red | d_se_red_pct | pct_win_win | spearman_soc_se | spearman_p | signif |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Wheat | cereal | 540 | 8.527 | 25.32 | 13.097 | 76.9 | 96.111 | -0.179 | 0 | *** |
| Maize | cereal | 650 | 8.974 | 26.852 | 30.794 | 69.2 | 97.538 | -0.232 | 0 | *** |
| Soybeans | tillage_only | 590 | 0.739 | 1.692 | 23.335 | 65 | 72.712 | -0.294 | 0 | *** |
| Cotton | tillage_only | 497 | 0.633 | 1.701 | 26.259 | 65 | 72.233 | -0.205 | 0 | *** |


## What drives the benefit (median component effect)

`*_residue` / `*_tillage` = the SOC gain / erosion avoided attributable to residue retention vs reduced tillage alone.


| commodity | switch | soc_tillage | se_tillage | soc_residue | se_residue | soc_driver | se_driver |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Wheat | cereal | 0.452 | 10.705 | 6.574 | 2.704 | residue | tillage |
| Maize | cereal | 0.292 | 28.947 | 7.969 | 5.344 | residue | tillage |
| Soybeans | tillage_only | 0.73 | 24.23 |  |  | tillage | tillage |
| Cotton | tillage_only | 0.633 | 27.202 |  |  | tillage | tillage |


## Where to focus — top wheat ecoregions by co-benefit priority

| region_name | country_name | biome | realm | d_soc | d_se_red | priority |
| --- | --- | --- | --- | --- | --- | --- |
| Llanos | <NA> | Tropical & Subtropical Grasslands, Savannas & Shrublands | Neotropic | 46.789 | 114.63 | 0.953 |
| Negro-Branco moist forests | <NA> | Tropical & Subtropical Moist Broadleaf Forests | Neotropic | 20.967 | 162.045 | 0.942 |
| Nansei Islands subtropical evergreen forests | <NA> | Tropical & Subtropical Moist Broadleaf Forests | Indomalayan | 27.697 | 98.486 | 0.914 |
| Honshu alpine conifer forests | <NA> | Temperate Conifer Forests | Palearctic | 18.74 | 121.29 | 0.913 |
| Nihonkai montane deciduous forests | <NA> | Temperate Broadleaf & Mixed Forests | Palearctic | 20.051 | 110.535 | 0.91 |
| Northeast Himalayan subalpine conifer forests | <NA> | Temperate Conifer Forests | Palearctic | 26.265 | 89.978 | 0.906 |
| Hainan Island monsoon rain forests | <NA> | Tropical & Subtropical Moist Broadleaf Forests | Indomalayan | 18.432 | 107.538 | 0.894 |
| Rock and Ice | <NA> |  |  | 16.729 | 116.736 | 0.885 |
| Oaxacan montane forests | <NA> | Tropical & Subtropical Moist Broadleaf Forests | Neotropic | 15.272 | 133.071 | 0.875 |
| Nihonkai evergreen forests | <NA> | Temperate Broadleaf & Mixed Forests | Palearctic | 15.847 | 113.795 | 0.871 |
| Bolivian Yungas | <NA> | Tropical & Subtropical Moist Broadleaf Forests | Neotropic | 14.205 | 165.351 | 0.869 |
| Western Himalayan alpine shrub and meadows | <NA> | Montane Grasslands & Shrublands | Palearctic | 15.563 | 115.183 | 0.867 |


Maps rendered this run: **5** (0 ⇒ no geometry available; maps are drop-in).


Regenerate with `python -m sbtn_leaf.claude_analysis.run_practice_change` or `sbtn_leaf.claude_analysis.practice_change.run_practice_change_analysis()`.
