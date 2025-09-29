# SBTN Land's Land Environmental Assessment Factors (LEAFs)

## Project Purpose
The SBTN Land LEAFs project provides reference factors and reproducible workflows to help organizations evaluate how land
management practices affect soil health and ecosystem resilience. It centralizes the source data, theory documentation,
notebooks, and processing code required to derive Soil Organic Carbon (SOC), Soil Erosion, and other land-focused LEAF
outputs.

## Repository Structure
- [`data/`](data/): Lightweight input datasets, starter templates, and sample configurations used by the tutorials and
  processing pipelines.
- [`documentation/`](documentation/): Theory references and step-by-step guides, including the
  [SOC documentation](documentation/SOC_Documentation.md), [Soil Erosion guide](documentation/Soil_Erosion_Documentation.md),
  and data harmonization notes.
- [`examples/`](examples/): Executable notebooks and scripts that demonstrate how to prepare inputs and compute LEAF values.
- [`src/`](src/): Reusable Python modules that implement the harmonization, calculation, and export utilities used across the
  examples.
- [`LEAFs/`](LEAFs/): Output directory where generated LEAF tables, maps, and reports are stored after running the processing
  workflows. See [LEAFs/README.md](LEAFs/README.md) for guidance on organizing and versioning run artifacts.

## Data and Documentation Resources
### Lightweight datasets
Input samples and reduced-size datasets are stored in [`data/`](data/). Each subfolder includes a `README` or inline metadata
to explain provenance and preprocessing steps so you can substitute your own data when scaling up analyses.

### Theory and methodology documentation
Detailed theory, assumptions, and methodology walkthroughs for each LEAF category live in [`documentation/`](documentation/).
Start with the [SOC documentation](documentation/SOC_Documentation.md) and [Soil Erosion documentation](documentation/Soil_Erosion_Documentation.md)
for practice-specific instructions, and consult the [data harmonization notes](documentation/Data_Harmonization.md) for
information on preparing shared inputs.

### Example notebooks and scripts
The [`examples/`](examples/) folder contains Jupyter notebooks and Python scripts that mirror the workflows described in the
documentation. Each example references the lightweight datasets in [`data/`](data/) and saves results into
[`LEAFs/`](LEAFs/).

### LEAF outputs
Run artifacts—including intermediate harmonized datasets, final factor tables, and any exported visualizations—are written to
[`LEAFs/`](LEAFs/). Organize outputs by experiment name or date to keep multiple iterations side by side. When collaborating,
consult [LEAFs/README.md](LEAFs/README.md) for recommendations on which outputs to version (for example lightweight CSV
summaries such as [`sample_soc_summary.csv`](LEAFs/sample_soc_summary.csv)) and how to keep large rasters or proprietary
datasets out of Git.

## Quick Start
1. Install the project dependencies:
   ```bash
   pip install -e .
   ```
2. Explore the theory and data preparation guidance in [`documentation/`](documentation/).
3. Open an example notebook from [`examples/`](examples/) and update the configuration cells to point at the desired inputs
   inside [`data/`](data/).
4. Execute the notebook or script. Generated factors and diagnostic charts will appear under [`LEAFs/`](LEAFs/).

## Running the Examples
To work through the notebooks:
1. Launch Jupyter Lab or Notebook from the project root:
   ```bash
   jupyter lab
   ```
2. Open one of the example notebooks in [`examples/`](examples/) and follow the embedded instructions. Each example
   demonstrates how to load the lightweight datasets, call the processing utilities from [`src/`](src/), and export
   results.
3. After the run completes, inspect the generated outputs in [`LEAFs/`](LEAFs/) to validate results or share them with your
   team. Larger derived datasets should be added to `.gitignore` before committing.

## Contributing
We welcome improvements to the data pipelines, documentation, and examples. To contribute:
1. Fork the repository and create a feature branch.
2. Make your changes, ensuring new notebooks or datasets are placed in the appropriate top-level folders.
3. Add or update documentation in [`documentation/`](documentation/) and reference new assets from the relevant sections
   above.
4. Run the provided examples or tests to confirm your changes, then open a pull request describing the updates and any
   expected LEAF outputs that reviewers should verify in [`LEAFs/`](LEAFs/).
