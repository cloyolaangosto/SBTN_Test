# Cross-indicator ecoregion aggregation — findings

How three SBTN-Land LEAFs — **SOC**, **soil erosion** and **terrestrial acidification** — behave as polygons shrink country → subcountry → ecoregion, and how much ecological grouping (biome / realm) explains.

Indicators and units:

| indicator | unit | flows | direction |
| --- | --- | --- | --- |
| SOC stock | t SOC/ha | commodity | higher = more carbon stored |
| Soil erosion | t soil/ha/yr | commodity | higher = worse |
| Acidification CF | kg SO2-eq./kg | pollutant | higher = worse |


## Coverage by level (share of region × flow cells with a value)

| indicator | level | n_regions | n_flows | coverage_pct |
| --- | --- | --- | --- | --- |
| SOC stock | country | 276 | 110 | 48 |
| SOC stock | subcountry | 3422 | 110 | 51.1 |
| SOC stock | ecoregion | 796 | 110 | 63.2 |
| Soil erosion | country | 276 | 106 | 51.4 |
| Soil erosion | subcountry | 3422 | 106 | 58 |
| Soil erosion | ecoregion | 847 | 106 | 64.9 |
| Acidification CF | country | 276 | 3 | 100 |
| Acidification CF | subcountry | 3422 | 3 | 100 |
| Acidification CF | ecoregion | 829 | 3 | 100 |


## Ecoregion significance (median across focal flows)

`eco_biome_eta2` / `eco_realm_eta2` = share of ecoregion-level variance explained by biome / realm (higher ⇒ ecoregions carry signal political units miss). `subcty_within_country_frac` = share of admin-1 variance a single national LEAF hides.


| indicator_name | n_focal_flows | eco_biome_eta2 | eco_realm_eta2 | subcty_within_country_frac |
| --- | --- | --- | --- | --- |
| SOC stock | 9 | 0.167 | 0.171 | 0.323 |
| Soil erosion | 9 | 0.26 | 0.156 | 0.385 |
| Acidification CF | 3 | 0.323 | 0.378 | 0.252 |


## Level reframing of the central estimate (country mean = 1.0)

| indicator | indicator_name | unit | rel_mean_country | rel_median_country | rel_mean_subcountry | rel_median_subcountry | rel_mean_ecoregion | rel_median_ecoregion |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soc | SOC stock | t SOC/ha | 1 | 1 | 0.982 | 0.975 | 1.067 | 1.033 |
| soil_erosion | Soil erosion | t soil/ha/yr | 1 | 1 | 0.893 | 0.887 | 1.048 | 1.045 |
| acidification | Acidification CF | kg SO2-eq./kg | 1 | 1 | 1.193 | 1.263 | 1.164 | 1.104 |


## Between-region dispersion (median CV by level)

| indicator | indicator_name | cv_country | cv_subcountry | cv_ecoregion |
| --- | --- | --- | --- | --- |
| soc | SOC stock | 0.393 | 0.386 | 0.387 |
| soil_erosion | Soil erosion | 1.338 | 1.505 | 1.294 |
| acidification | Acidification CF | 0.979 | 0.846 | 0.952 |


Tables live in `tables/`; figures in `figures/`. Regenerate everything with `python -m sbtn_leaf.claude_analysis.run_cross_indicator` or `sbtn_leaf.claude_analysis.cross_indicator.run_cross_indicator_analysis()`.
