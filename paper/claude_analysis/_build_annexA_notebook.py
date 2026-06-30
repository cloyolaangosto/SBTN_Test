"""Generate 'Annex A - RothC further methodological details' (Supplementary Info)."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def md(*lines):
    text = "\n".join(lines)
    parts = text.split("\n")
    src = [p + "\n" for p in parts[:-1]] + [parts[-1]]
    return {"cell_type": "markdown", "metadata": {}, "source": src}


def code(*lines):
    text = "\n".join(lines)
    parts = text.split("\n")
    src = [p + "\n" for p in parts[:-1]] + [parts[-1]]
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": src}


cells = []

cells.append(md(
    "# Annex A — RothC further methodological details",
    "",
    "*Supplementary Information referenced in Methods → SOC: “Plant residue and other needed",
    "calculations for each commodity type are described in detail in Supplementary Information.”*",
    "",
    "This annex documents (A1) the RothC carbon model as implemented here, (A2) pool initialisation,",
    "(A3) decomposition and carbon partitioning, (A4) the rate-modifying factors, (A5) **plant-residue",
    "and carbon-input calculations for each commodity type** (annual & permanent crops, forests,",
    "grasslands), (A6) the reduced-tillage modification, and (A7–A8) inputs, constants and a code map.",
    "Every formula and constant is taken directly from the repository implementation",
    "(`src/sbtn_leaf/RothC_Core.py` and `src/sbtn_leaf/cropcalcs.py`); the parameter tables below are",
    "loaded live from the bundled data so this annex is fully reproducible.",
))

# ---- A1 ----
cells.append(md(
    "## A1. The RothC model",
    "",
    "SOC is simulated with the Rothamsted Carbon model (RothC; Coleman, Prout & Milne, 2024), which",
    "partitions soil organic carbon into four decomposing pools and one inert pool:",
    "",
    "- **DPM** — decomposable plant material",
    "- **RPM** — resistant plant material",
    "- **BIO** — microbial biomass",
    "- **HUM** — humified organic matter",
    "- **IOM** — inert organic matter (does not decompose)",
    "",
    "with `SOC = DPM + RPM + BIO + HUM + IOM`. The model is run in **monthly steps (Δt = 1/12 yr)** for",
    "14 years (2016–2030). Each month, the four active pools decay by first-order kinetics, the",
    "respired fraction leaves as CO₂, the remainder is redistributed to BIO and HUM, and fresh carbon",
    "from plant residues (and any manure) is added to DPM/RPM. Decay is modulated by temperature,",
    "moisture, soil cover and — under conservation tillage — a tillage modifier.",
))

# ---- A2 ----
cells.append(md(
    "## A2. Pool initialisation from SOC and clay",
    "",
    "Initial pools are derived from total SoilGrids SOC (`soc`, t C/ha) and clay content (`clay`, %)",
    "using the standard RothC/IPCC allometric partitioning (`RothC_Core.initialize_pools`):",
    "",
    "```",
    "IOM = 0.049 · soc**1.139",
    "RPM = (0.1847·soc + 0.1555) · (clay + 1.275)**(-0.1158)",
    "HUM = (0.7148·soc)         · (clay + 0.3421)**( 0.0184)",
    "BIO = (0.014·soc + 0.0075) · (clay + 8.8473)**( 0.0567)",
    "DPM = soc - IOM - RPM - HUM - BIO",
    "```",
))

# ---- A3 ----
cells.append(md(
    "## A3. Decomposition and carbon partitioning",
    "",
    "Each pool decays exponentially with a baseline annual rate constant `k` modulated by the combined",
    "rate modifier `rate_m = a·b·c` (temperature × moisture × plant cover; §A4) and, for reduced",
    "tillage, a per-pool tillage modifier `TRM` (§A6):",
    "",
    "```",
    "pool(t+1) = pool(t) · exp( -rate_m · TRM_pool · k_pool · Δt )",
    "k_DPM = 10.0   k_RPM = 0.3   k_BIO = 0.66   k_HUM = 0.02   (yr^-1)   # TRM_pool = 1 under conventional tillage",
    "```",
    "",
    "The carbon lost from each pool is split between CO₂ (respired) and the (BIO+HUM) pools using a",
    "clay-dependent ratio:",
    "",
    "```",
    "x          = 1.67 · (1.85 + 1.60 · exp(-0.0786·clay))   # CO2 : (BIO+HUM) ratio",
    "resp_frac  = x / (x + 1)        # fraction of each pool's loss respired as CO2",
    "to_BIOHUM  = 1 / (x + 1)        # fraction retained, of which 46% -> BIO and 54% -> HUM",
    "```",
    "",
    "Fresh carbon inputs are then added. Plant-residue carbon `C_in` (§A5) is split between DPM and RPM",
    "by the commodity-specific **DPM/RPM ratio** `r`; farmyard manure (`FYM`, where applicable) follows",
    "fixed RothC fractions:",
    "",
    "```",
    "DPM += r/(r+1) · C_in        RPM += 1/(r+1) · C_in",
    "DPM += 0.49·FYM   RPM += 0.49·FYM   HUM += 0.02·FYM",
    "```",
    "",
    "DPM/RPM ratio by commodity type (controls the decomposability of incoming residue):",
    "",
    "| commodity type | DPM/RPM | % DPM / % RPM |",
    "|---|---|---|",
    "| annual crops | 1.44 | 59 / 41 |",
    "| grassland | 1.44 | 59 / 41 |",
    "| permanent crops | 1.00 | 50 / 50 |",
    "| forest (deciduous / tropical woodland) | 0.25 | 20 / 80 |",
))

# ---- A4 ----
cells.append(md(
    "## A4. Rate-modifying factors",
    "",
    "**Temperature** `a` (°C; `RMF_Tmp`), zero below −5 °C:",
    "```",
    "a = 47.91 / ( exp(106.06 / (T + 18.27)) + 1 )      for T >= -5 ; a = 0 otherwise",
    "```",
    "",
    "**Moisture** `b` (`RMF_Moist`) tracks a soil-moisture deficit (SMD). The maximum deficit for a",
    "23 cm layer is scaled to the modelled depth, open-pan evaporation is multiplied by 0.75, and the",
    "soil-water content (SWC) is updated and bounded each month:",
    "```",
    "SMDmax  = -(20 + 1.3·clay - 0.01·clay**2)          # for 23 cm; scaled by depth/23",
    "b = 1.0                                            if SWC > 0.444·SMD",
    "b = 0.2 + 0.8 · (SMD - SWC) / (SMD - 0.444·SMD)    otherwise",
    "```",
    "",
    "**Plant cover** `c` (`RMF_PC`) slows decomposition under a canopy:",
    "```",
    "c = 1.0  (bare soil)        c = 0.6  (vegetated)",
    "```",
    "The monthly bare/vegetated state per commodity comes from the crop-cover window",
    "(`create_plant_cover_monthly_curve`, derived from the crop-coefficient phenology, §A5.2).",
))

# ---- A5 ----
cells.append(md(
    "## A5. Plant-residue and carbon-input calculations by commodity type",
    "",
    "This is the calculation the Methods refers to. Annual carbon input to RothC is the residue carbon",
    "returned to the soil; it is computed differently for crops, forests and grasslands, then",
    "distributed across the 12 months.",
))

cells.append(md(
    "### A5.1 Annual & permanent crops — yield → residue carbon",
    "",
    "Per pixel, the commodity yield (t/ha; from FAO-GAEZ v5 suitability-scaled yields, with FAO/SPAM",
    "fallbacks) is converted to residue carbon via `cropcalcs.calculate_crop_residues`, which selects",
    "one of three branches depending on the data available for the crop:",
    "",
    "```",
    "# Branch 1 - allometric regression (crops with an above-ground regression):",
    "ABG = slope · yield + intercept ;  BG = RS · ABG",
    "# Branch 2 - above-ground ratio (crops with R_AG and a root:shoot ratio RS > 0):",
    "ABG = R_AG · yield ;               BG = RS · ABG",
    "# Branch 3 - total-residue ratio (otherwise):",
    "Res = R_T · yield · DRY · C_content",
    "",
    "# Branches 1-2 then convert dry matter to carbon:",
    "Res = (ABG + BG) · DRY · C_content",
    "```",
    "",
    "where `DRY` is dry-matter fraction, `RS` root:shoot ratio, `R_AG`/`R_T` residue-to-yield ratios,",
    "and `C_content` the carbon fraction of dry matter (default **0.5** for crops). The parameter",
    "tables are loaded below.",
))

cells.append(code(
    "import warnings; warnings.filterwarnings('ignore')",
    "from sbtn_leaf import data_loader as dl",
    "",
    "# Residue-to-yield ratios per crop: R_AG (above-ground:yield), RS (root:shoot),",
    "# DRY (dry-matter fraction), R_T (total residue:yield, used only when RS == 0).",
    "dl.get_crop_residue_ratio_table()",
))

cells.append(code(
    "# Above-ground residue regressions (ABG = Slope*yield + Intercept), in t/ha dry matter.",
    "dl.get_crop_ag_residue_table()",
))

cells.append(md(
    "### A5.2 Monthly distribution of crop residue carbon",
    "",
    "Annual residue carbon is spread across months from the crop's phenology (planting date and cycle",
    "length by thermal-climate zone, from the crop-coefficient table; harvest month derived via the",
    "absolute-day table). The split concentrates inputs around harvest",
    "(`cropcalcs._distribute_residue_monthly`):",
    "",
    "```",
    "Annual crops:    50% at the harvest month ; remaining 50% split equally over the 3 months before",
    "                 harvest (1/6 each); 0 elsewhere.",
    "Permanent crops: 70% at the harvest month ; remaining 30% split equally over the 4 months before",
    "                 harvest (7.5% each); 0 elsewhere.",
    "```",
))

cells.append(md(
    "### A5.3 Forests — age-dependent litter",
    "",
    "Forest carbon input is litterfall, which rises with stand age towards a mature rate over a",
    "transition period `TP` (`cropcalcs.get_forest_litter_monthlyrate_fromda`):",
    "",
    "```",
    "litter_annual = min( res_rate ,  (res_rate / TP) · (forest_age + offset) )",
    "litter_month  = litter_annual / 12",
    "```",
    "",
    "`res_rate` (mature litter, t C/ha/yr) and `TP` (years) are looked up by IPCC climate and forest",
    "type — **BD** = broadleaf-deciduous, **NE** = needle-leaf-evergreen — with `_mean/_min/_max`",
    "supporting an optional triangular-distribution sampling. The IPCC default transition period",
    "(`IPCC_TP`) is 20 years. The lookup table:",
))

cells.append(code(
    "import polars as pl",
    "from sbtn_leaf.paths import data_path",
    "",
    "# Forest litter rates (t C/ha/yr) and transition periods (yr) by IPCC climate and forest type.",
    "pl.read_excel(data_path('forest', 'forest_residues_IPCC.xlsx'))",
))

cells.append(md(
    "### A5.4 Grasslands",
    "",
    "Grassland carbon input is the above- plus below-ground residue, looked up by climate zone and",
    "scaled to carbon (`cropcalcs.generate_grassland_residue_map`):",
    "",
    "```",
    "Res = (Residue_Above + Residue_Below) · 0.5 · C_content      # C_content default 0.47",
    "```",
    "",
    "with optional per-pixel stochastic sampling around the tabulated standard errors. (Grazing-animal",
    "variants — goat, sheep, cattle and combinations — give the six grassland flows reported in the",
    "main text.) The grassland residue lookup is climate-zone-indexed in the SOC support workflow.",
))

# ---- A6 ----
cells.append(md(
    "## A6. Reduced-tillage modification (Hyun & Yoo, 2024)",
    "",
    "Conservation tillage is represented by per-pool **tillage-rate modifiers (TRM)** that multiply the",
    "baseline decay constants `k` (§A3). Following Hyun & Yoo (2024), each pixel is classified into one",
    "of four nodes by a decision tree on **sand content (%)** and **SOC (t C/ha)**, and the four pools",
    "receive node-specific multipliers (`RothC_Core.RMF_TRM`):",
    "",
    "```",
    "if sand > 37.6:  node = 1 if SOC > 75.7 else 2",
    "else:            node = 3 if sand > 35.0 else 4",
    "```",
    "",
    "A TRM > 1 accelerates and a TRM < 1 slows that pool's decomposition under reduced tillage,",
    "depending on soil texture and carbon status. The coefficient matrix (rows = pools, columns =",
    "nodes 1–4):",
))

cells.append(code(
    "import pandas as pd",
    "from sbtn_leaf.RothC_Core import _TRM_COEFFICIENTS",
    "",
    "pd.DataFrame(",
    "    _TRM_COEFFICIENTS,",
    "    index=['DPM', 'RPM', 'BIO', 'HUM'],",
    "    columns=['node1 (sand>37.6, SOC>75.7)', 'node2 (sand>37.6, SOC<=75.7)',",
    "             'node3 (35<sand<=37.6)', 'node4 (sand<=35)'],",
    ")",
))

# ---- A7 ----
cells.append(md(
    "## A7. Inputs and key constants",
    "",
    "| input | source | resolution |",
    "|---|---|---|",
    "| monthly precipitation | NASA GPM IMERG Final L3 (Huffman et al., 2019), 2014–2023 mean | 0.1° → UHTU |",
    "| monthly temperature | GLDAS Catchment LSM L4 (Li et al., 2020), 2014–2023 mean | 1° → UHTU |",
    "| SOC (0–30 cm), clay & sand (15–30 cm) | SoilGrids 2.0 (Poggio et al., 2021) | → UHTU |",
    "| crop growing area / yield | FAO-GAEZ v5 crop suitability (FAO & IIASA, 2025) | → UHTU |",
    "| forest / grassland area | SOC raster of Morais et al. (2019) | → UHTU |",
    "",
    "| constant | value | role |",
    "|---|---|---|",
    "| k (DPM, RPM, BIO, HUM) | 10.0, 0.3, 0.66, 0.02 yr⁻¹ | baseline decay |",
    "| BIO : HUM split of retained C | 46 % : 54 % | partition of non-respired loss |",
    "| FYM split (DPM/RPM/HUM) | 0.49 / 0.49 / 0.02 | manure partition |",
    "| C content (crops / grassland) | 0.5 / 0.47 | dry matter → carbon |",
    "| plant-cover modifier (bare / veg) | 1.0 / 0.6 | decomposition rate |",
    "| open-pan evaporation scalar | 0.75 | moisture factor |",
    "| run length / step | 14 yr (2016–2030) / monthly | simulation |",
    "| averaging | winsorized at 1st/99th percentile per region | country / sub-country / ecoregion |",
))

# ---- A8 ----
cells.append(md(
    "## A8. Implementation map",
    "",
    "| component | function | file |",
    "|---|---|---|",
    "| pool initialisation | `initialize_pools` | `src/sbtn_leaf/RothC_Core.py` |",
    "| monthly decomposition + C partition | `decomp` | `src/sbtn_leaf/RothC_Core.py` |",
    "| temperature / moisture / cover modifiers | `RMF_Tmp`, `RMF_Moist`, `RMF_PC` | `src/sbtn_leaf/RothC_Core.py` |",
    "| tillage modifier | `RMF_TRM`, `_TRM_COEFFICIENTS` | `src/sbtn_leaf/RothC_Core.py` |",
    "| equilibrium / simulation | `run_equilibrium`, `run_simulation` | `src/sbtn_leaf/RothC_Core.py` |",
    "| crop residue from yield | `calculate_crop_residues` | `src/sbtn_leaf/cropcalcs.py` |",
    "| monthly residue distribution | `_distribute_residue_monthly` | `src/sbtn_leaf/cropcalcs.py` |",
    "| forest litter | `get_forest_litter_monthlyrate_fromda` | `src/sbtn_leaf/cropcalcs.py` |",
    "| grassland residue | `generate_grassland_residue_map` | `src/sbtn_leaf/cropcalcs.py` |",
    "| raster RothC driver | `run_RothC_crops` / `_forest` / `_grassland` | `src/sbtn_leaf/RothC_Raster.py` |",
    "",
    "*References: Coleman, Prout & Milne (2024); Hyun & Yoo (2024); Morais, Teixeira & Domingos (2019);",
    "Teixeira, Morais & Domingos (2021); Huffman et al. (2019); Li et al. (2020); Poggio et al. (2021);",
    "FAO & IIASA (2025) — full entries in the main reference list.*",
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

out = HERE / "AnnexA_RothC_Methodological_Details.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n")
print("wrote", out)
