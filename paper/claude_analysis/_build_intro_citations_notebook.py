"""Generate the Introduction-citations proposal notebook.

Build helper: assembles the cells and writes the .ipynb. The notebook is pure
markdown (a citation proposal for the manuscript Introduction), so no execution
is required, but it is nbconvert-validated like the others.
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def md(*lines):
    return {"cell_type": "markdown", "metadata": {}, "source": _src(lines)}


def _src(lines):
    text = "\n".join(lines)
    parts = text.split("\n")
    return [p + "\n" for p in parts[:-1]] + [parts[-1]]


cells = []

cells.append(md(
    "# Introduction — proposed citations for the `[citation]` placeholders",
    "",
    "The manuscript Introduction (pp. 1–3 of *LEAF_Manuscript_Jun29*) carries **17 unresolved",
    "citation markers** — `[citation]`, `[citations]`, `[paper_1, paper_2, paper_3]` and",
    "`[citation of the flows]`. This notebook proposes a reference for each, quoting the sentence it",
    "supports, and flags how confident the suggestion is.",
    "",
    "Confidence legend:",
    "",
    "| flag | meaning |",
    "|---|---|",
    "| **✓ verified** | the reference is real and topically apt — checked against the literature |",
    "| **◐ canonical** | a standard, widely-cited source for this claim (high confidence, not separately web-checked) |",
    "| **○ suggested** | plausible candidate the authors should confirm or swap for their preferred source |",
    "",
    "> **How to read this.** Several Introduction sentences are *framing* statements (e.g. the “five",
    "> options” taxonomy) for which there is no single canonical paper; for those I propose a candidate",
    "> and say so. The science-method claims (GWP, regionalised LCA, land-use LCIA, terrestrial",
    "> acidification) map onto well-established references and are given with DOIs where I could confirm",
    "> them. Nothing here was invented — uncertain DOIs are marked “confirm”.",
))

# --------------------------------------------------------------------------- #
cells.append(md(
    "## 1. The five options & the practice-change counterfactual (p. 1)",
    "",
    "> *“They can explore source reduction or efficiency improvements (Raimbekov et al., 2023),",
    "> change inputs (e.g., feedstocks) **[citation]**, change technologies **[citation]**, change",
    "> sourcing geographies **[citation]**, or change production practices (Bocken et al. 2014). … The",
    "> fifth option … requires a quantified counterfactual upon which companies can make informed",
    "> sustainability decisions **[citation]**.”*",
    "",
    "These four slots sit inside one taxonomy of corporate impact-reduction levers. There is **no",
    "single canonical paper** for “change inputs / technologies / geographies” individually — the",
    "cleanest fix is to anchor the whole sentence to one framework and drop the per-item markers.",
    "",
    "| slot | proposed reference | flag |",
    "|---|---|---|",
    "| change inputs (feedstocks) | Bocken, N.M.P., Short, S.W., Rana, P., Evans, S. (2014). *J. Cleaner Production* 65, 42–56 — the “sustainable business model archetypes” already cited two clauses later cover input/technology substitution; reuse it as the umbrella, **or** cite a cleaner-production strategy review. | ○ suggested |",
    "| change technologies | Bocken et al. (2014) (technology archetypes) **or** a cleaner-production / eco-efficiency reference. | ○ suggested |",
    "| change sourcing geographies | a sustainable-sourcing / supply-chain reconfiguration reference; Raimbekov et al. (2023, already cited) touches supply channels. | ○ suggested |",
    "| informed decisions / quantified counterfactual | Greenhouse Gas Protocol (2011), *Corporate Value Chain (Scope 3) Accounting and Reporting Standard*, WRI & WBCSD — the standard reference for baselines/counterfactuals in corporate accounting. | ◐ canonical |",
    "",
    "**Recommendation:** rephrase to a single framework citation for the five options (e.g. Bocken et",
    "al. 2014 as the umbrella) rather than four separate markers; keep a counterfactual/accounting",
    "citation (GHG Protocol Scope 3) for the last clause.",
))

# --------------------------------------------------------------------------- #
cells.append(md(
    "## 2. Climate accounting: 1.5 °C, the CF/GWP, and location-independence (p. 1)",
    "",
    "> *“… align corporate climate targets with the objective to restrict global warming to under",
    "> 1.5 °C **[citation]**. … translated into kg CO2-eq by the use of characterization factor (CF)",
    "> **[citation]**. … the GWP of 1 kg of CH₄ is of 28.4 kg of CO2-eq … this CF is independent of",
    "> where it’s emitted: GHGs emitted in China will have the same effect than those emitted in",
    "> Hawaii **[citation]**.”*",
    "",
    "| slot | proposed reference | flag |",
    "|---|---|---|",
    "| restrict warming to < 1.5 °C | UNFCCC (2015), *Paris Agreement* — and/or IPCC (2018), *Global Warming of 1.5 °C* (Special Report SR1.5), V. Masson-Delmotte et al. (eds.), Cambridge University Press. | ◐ canonical |",
    "| characterization factor / GWP (CH₄ = 28.4) | IPCC (2021), *Climate Change 2021: The Physical Science Basis* (AR6 WG1), **Forster et al., Chapter 7** (GWP-100 table, Table 7.15 — biogenic CH₄ ≈ 27–30). Optionally the LCIA-CF concept: Hauschild, Rosenbaum & Olsen (2018), *Life Cycle Assessment: Theory and Practice*, Springer. | ✓ verified |",
    "| GWP is location-independent (well-mixed GHGs) | the regionalisation literature contrasting global vs local impacts — Hauschild (2006), *Int. J. LCA* 11, 11–13 (“Spatial differentiation in LCA”) or Patouillard et al. (2018, below). | ◐ canonical |",
    "",
    "*Note:* the value **28.4** matches IPCC AR6’s methane GWP-100 (~28 biogenic / ~30 fossil, including",
    "the oxidation term); citing AR6 WG1 Ch. 7 makes the exact figure traceable.",
))

# --------------------------------------------------------------------------- #
cells.append(md(
    "## 3. Nature dependency & the maturity gap in land accounting (p. 2)",
    "",
    "> *“… companies … rely … on the ecological condition of land, freshwater, oceans, and",
    "> biodiversity **[citation]**. … barely developed for some land impacts … and slightly more",
    "> developed for terrestrial acidification **[citation]** … 1) lack of … impact assessment methods",
    "> **[citation]** and 2) lack of regionalized data … to support them **[citation]**.”*",
    "",
    "| slot | proposed reference | flag |",
    "|---|---|---|",
    "| economy depends on nature & biodiversity | Dasgupta, P. (2021), *The Economics of Biodiversity: The Dasgupta Review*, HM Treasury, London (“our economies are embedded within nature”). Optionally IPBES (2019), *Global Assessment Report on Biodiversity and Ecosystem Services*, and/or Rockström et al. (2009)/Steffen et al. (2015) planetary boundaries. | ✓ verified |",
    "| acidification more developed than SOC/erosion | the LCIA best-practice maturity assessment — Hauschild et al. (2013), *Int. J. LCA* 18, 683–697 (“Identifying best existing practice for characterization modeling in LCIA”). | ◐ canonical |",
    "| lack of implemented (land/soil) impact-assessment methods | Koellner et al. (2013), *Int. J. LCA* 18(6), 1188–1202 — UNEP-SETAC land-use LCIA guideline; and/or Curran et al. (2016) on land-use/biodiversity LCIA gaps. | ✓ verified |",
    "| lack of regionalized data | Patouillard et al. (2018), *J. Cleaner Production* 177, 398–412 (critical review of spatial differentiation in LCA). | ✓ verified |",
))

# --------------------------------------------------------------------------- #
cells.append(md(
    "## 4. LCA, regionalised methods, acidification CFs, and land-use flows (pp. 2–3)",
    "",
    "> *“Life Cycle Assessment … the expansion of carbon accounting into … multi-indicator accounting",
    "> of products and services **[citation]**, has made … progress in … regionalized impact",
    "> assessment methods **[citations]**. Regionalized terrestrial acidification CFs … developed by",
    "> **[paper_1, paper_2, paper_3]** … Teixeira, Morais & Domingos (2021) developed … generalized LCA",
    "> land-use flows **[citation of the flows]** …”*",
    "",
    "| slot | proposed reference | flag |",
    "|---|---|---|",
    "| LCA = multi-indicator accounting (definition) | ISO 14040:2006 & ISO 14044:2006 (*Environmental management — Life cycle assessment*). Optionally Hauschild, Rosenbaum & Olsen (2018) or Guinée et al. (2011), *Environ. Sci. Technol.* 45, 90–96. | ◐ canonical |",
    "| regionalized impact-assessment methods | Mutel & Hellweg (2009), *Environ. Sci. Technol.* 43, 5797–5803 (DOI 10.1021/es803002j); Patouillard et al. (2018); Frischknecht et al. (2019) UNEP global guidance. | ✓ verified |",
    "| **[paper_1, paper_2, paper_3]** regionalized acidification CFs | Roy et al. (2012), *Environ. Sci. Technol.* 46, 8270–8278 (DOI 10.1021/es3013563); Roy et al. (2014), *Sci. Total Environ.* 500, 270–276; **Azevedo et al. (2013)**, *Environ. Pollution* 174, 10–15. (Roy 2012/2014 are already in the reference list.) Alternatives: van Zelm et al. (2007); Seppälä et al. (2006); Posch et al. (2008). | ✓ verified |",
    "| **[citation of the flows]** generalised land-use flows | Koellner et al. (2013), *Int. J. LCA* 18(6), 1203–1215, “Principles for life cycle inventories of land use on a global scale” (DOI 10.1007/s11367-013-0580-6) — the land-use flow typology; plus Teixeira et al. (2021, already cited) who applied it. | ✓ verified |",
))

# --------------------------------------------------------------------------- #
cells.append(md(
    "## 5. Interpretation difficulty & implementation barriers (p. 3)",
    "",
    "> *“… different impact assessment methods delivering different measures despite sharing similar",
    "> naming **[citations]**. … challenging implementation at a more granular level due to lack of",
    "> more regionalized flows in databases and rising computational needs **[citations]**.”*",
    "",
    "| slot | proposed reference | flag |",
    "|---|---|---|",
    "| same name, different measures across methods | Owsianiak et al. (2014), *Int. J. LCA* 19, 1007–1021 (comparison of ReCiPe / IMPACT 2002+ / ILCD); Hauschild et al. (2013). | ◐ canonical |",
    "| rising computational needs of regionalisation | Mutel & Hellweg (2009), *Environ. Sci. Technol.* 43, 5797–5803 (DOI 10.1021/es803002j); Patouillard et al. (2018). | ✓ verified |",
))

# --------------------------------------------------------------------------- #
cells.append(md(
    "## Paste-ready reference entries (additions only)",
    "",
    "Formatted to match the manuscript’s existing reference list. Entries already present (Bocken",
    "2014; Roy 2012/2014; Teixeira 2021; Raimbekov 2023) are not repeated.",
    "",
    "- Azevedo, L. B., van Zelm, R., Hendriks, A. J., Bobbink, R., & Huijbregts, M. A. J. (2013).",
    "  Global assessment of the effects of terrestrial acidification on plant species richness.",
    "  *Environmental Pollution*, 174, 10–15.",
    "- Curran, M., de Baan, L., De Schryver, A. M., et al. (2016). Toward meaningful end points of",
    "  biodiversity in life cycle assessment. *Environmental Science & Technology*, 50(6), 2782–2795.",
    "- Dasgupta, P. (2021). *The Economics of Biodiversity: The Dasgupta Review.* HM Treasury, London.",
    "- Frischknecht, R., Jolliet, O., et al. (Eds.) (2019). *Global Guidance for Life Cycle Impact",
    "  Assessment Indicators, Volume 2.* UNEP/SETAC Life Cycle Initiative.",
    "- Greenhouse Gas Protocol (2011). *Corporate Value Chain (Scope 3) Accounting and Reporting",
    "  Standard.* World Resources Institute & WBCSD.",
    "- Guinée, J. B., Heijungs, R., Huppes, G., et al. (2011). Life cycle assessment: past, present,",
    "  and future. *Environmental Science & Technology*, 45(1), 90–96.",
    "- Hauschild, M. Z. (2006). Spatial differentiation in life cycle impact assessment: a decade of",
    "  method development to increase the environmental realism of LCIA. *Int. J. LCA*, 11(S1), 11–13.",
    "- Hauschild, M. Z., Goedkoop, M., Guinée, J., et al. (2013). Identifying best existing practice",
    "  for characterization modeling in life cycle impact assessment. *Int. J. LCA*, 18(3), 683–697.",
    "- Hauschild, M. Z., Rosenbaum, R. K., & Olsen, S. I. (Eds.) (2018). *Life Cycle Assessment:",
    "  Theory and Practice.* Springer.",
    "- IPBES (2019). *Global Assessment Report on Biodiversity and Ecosystem Services.* IPBES",
    "  secretariat, Bonn, Germany.",
    "- IPCC (2018). *Global Warming of 1.5 °C* (Special Report). Masson-Delmotte, V., et al. (Eds.).",
    "  Cambridge University Press.",
    "- IPCC (2021). *Climate Change 2021: The Physical Science Basis* (AR6, WG1). Forster, P., et al.",
    "  Chapter 7. Cambridge University Press.",
    "- ISO (2006). *ISO 14040: Environmental management — Life cycle assessment — Principles and",
    "  framework.* International Organization for Standardization.",
    "- Koellner, T., de Baan, L., Beck, T., et al. (2013). UNEP-SETAC guideline on global land use",
    "  impact assessment on biodiversity and ecosystem services in LCA. *Int. J. LCA*, 18(6),",
    "  1188–1202.",
    "- Koellner, T., de Baan, L., Beck, T., et al. (2013). Principles for life cycle inventories of",
    "  land use on a global scale. *Int. J. LCA*, 18(6), 1203–1215.",
    "- Mutel, C. L., & Hellweg, S. (2009). Regionalized life cycle assessment: computational",
    "  methodology and application to inventory databases. *Environ. Sci. Technol.*, 43(15),",
    "  5797–5803.",
    "- Owsianiak, M., Laurent, A., Bjørn, A., & Hauschild, M. Z. (2014). IMPACT 2002+, ReCiPe 2008 and",
    "  ILCD’s recommended practice for characterization modelling … a case study-based comparison.",
    "  *Int. J. LCA*, 19(5), 1007–1021.",
    "- Patouillard, L., Bulle, C., Querleu, C., et al. (2018). Critical review and practical",
    "  recommendations to integrate the spatial dimension into life cycle assessment. *J. Cleaner",
    "  Production*, 177, 398–412.",
    "- UNFCCC (2015). *Paris Agreement.* United Nations Framework Convention on Climate Change.",
    "- van Zelm, R., Huijbregts, M. A. J., van Jaarsveld, H. A., et al. (2007). Time horizon dependent",
    "  characterization factors for acidification in life-cycle assessment based on … Europe.",
    "  *Environ. Sci. Technol.*, 41(3), 922–927.",
    "",
    "### Caveats",
    "",
    "1. The “five options” markers (change inputs / technologies / geographies) have **no single",
    "   canonical source**; the recommendation is to anchor the sentence to one framework citation",
    "   rather than four separate ones.",
    "2. DOIs are given only where confirmed; for ISO standards, IPCC/IPBES reports and the Dasgupta",
    "   Review, the institutional citation is the standard form.",
    "3. For `[paper_1, paper_2, paper_3]`, Roy et al. (2012) and (2014) are already in the reference",
    "   list — pick a third (Azevedo 2013 recommended) or expand to the fuller acidification-CF set.",
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

out = HERE / "Introduction_Citations.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n")
print("wrote", out)
