# Acidification: inter-biome vs inter-gas variability — findings

Tests the manuscript's claim (p.8, Fig. 2): *“inter-biome variability is larger than inter-gas variability.”* A one-way η² makes the two look equal because each factor's effect falls into the other's residual; the honest comparison is a **two-way (biome × gas)** decomposition of the ecoregion-level CFs.


## Two-way variance decomposition

| scale | eta2_biome | eta2_gas | eta2_interaction | eta2_residual | biome_gas_ratio | biome_median_spread | gas_median_spread |
| --- | --- | --- | --- | --- | --- | --- | --- |
| log10 | 0.2177 | 0.2202 | 0.0138 | 0.5482 | 0.9885 | 0.9421 | 0.6879 |
| linear | 0.1871 | 0.148 | 0.0896 | 0.5754 | 1.2645 | 8.7522 | 4.8742 |


## One-way biome η² within each gas (log10)

| gas | gas_label | n | biome_eta2_log | kruskal_p | signif |
| --- | --- | --- | --- | --- | --- |
| acid_nh3 | NH₃ | 828 | 0.2735 | 0 | *** |
| acid_nox | NOₓ | 828 | 0.3061 | 0 | *** |
| acid_so2 | SO₂ | 828 | 0.3231 | 0 | *** |


## Verdict

The claim is **defensible**. On a log scale the two main effects are comparable (biome η²=0.22 vs gas η²=0.22), but the spread of **biome medians** (0.94 dex) exceeds that of the gases (0.69 dex), and on the **absolute (linear)** scale biome dominates (η²=0.19 vs 0.15; median range 8.8× vs 4.9×). Within any single gas biome explains 0.27–0.32 of the variance, and the biome×gas interaction is small (η²=0.014), so the biome pattern is consistent across gases — a robust axis of variation, consistent with Roy et al. (2014).
