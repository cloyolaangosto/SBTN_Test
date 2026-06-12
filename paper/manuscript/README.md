# LEAFs manuscript

A full research-article manuscript describing the SBTN Land LEAFs dataset and the
`sbtn_leaf` pipeline, covering all three soil-quality indicators (soil organic carbon,
soil erosion, terrestrial acidification).

## Files
- `main.tex` — the manuscript source (LaTeX, `article` class).
- `references.bib` — BibTeX bibliography.
- `figures/` — figures referenced by the manuscript.

## Compiling
Requires a standard TeX distribution (e.g. TeX Live) with `natbib`, `mhchem`,
`siunitx`, `booktabs`, `authblk`, and `hyperref`.

```bash
cd paper/manuscript
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

This produces `main.pdf`.

## Figures: real vs. placeholder
The following figures are real renders copied from the repository:
- `K_Curve_example.png` — crop-coefficient curve (from `documentation/support_files/`).
- `SOC_LEAF_Example_25_0.png`, `SOC_LEAF_Example_31_2.png` — SOC map and time series
  (from `examples/SOC_LEAF_Example_files/`). Other `SOC_LEAF_Example_*.png` renders are
  also copied here and can be swapped into the manuscript if preferred.

One figure is a **labelled placeholder** because it exists only as code in the analysis
notebook:
- Figure "acidification by realm" (boxplots) — regenerate by running
  `paper/Paper_Graphs_Figures.ipynb` (section *Terrestrial Acidification*), export the
  boxplot to `figures/acid_by_realm.png`, and replace the `\fbox{...}` placeholder in
  `main.tex` with `\includegraphics{acid_by_realm.png}`.

The 20 MB `Paper_Graphs_Figures.ipynb` was not executed when preparing this manuscript;
the reported numbers in the Results section were computed directly from the CSV tables in
`LEAFs/`.

## Numbers in the Results section
All quantitative statements (SOC mean/median, erosion mean/median, acidification factors
by gas and realm, reduced-tillage gain) were computed from the published factor tables in
`LEAFs/SOC/`, `LEAFs/soil_erosion/`, and `LEAFs/acidification/`. Re-running those
aggregations will reproduce the figures cited in the text.
