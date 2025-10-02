# SBTN Land's Land Environmental Assessment Factors (LEAFs)

## Project Purpose
The SBTN Land LEAFs repository provides default ecoregion factors for a series of land uses and land managment options to help them evaluate against ecoregion's threshold and bring them into safe operating conditions, being part of SBTN Land v2 Targets. Full documentation on how to apply LEAFs into corporate accounting can be found in [AGILE](https://sciencebasedtargetsnetwork.org/wp-content/uploads/2025/04/SBTN-Land-Accounting-Guidelines-Draft-for-Public-Consultation.pdf) and the target setting process in [SBTN Land v2 Targets](https://sciencebasedtargetsnetwork.org/wp-content/uploads/2025/04/SBTN-Step-3-Land-Technical-Guidance-V2-DRAFT.pdf).

The repository also offers a reproducible workflows to help organizations evaluate new land uses, new land management options, or incorporate primary data to develop customized LEAFs. It centralizes the source data, theoretical documentation, example notebooks, and processing code required to derive Soil Organic Carbon (SOC), and Soil Erosion LEAFs.

## Repository Structure
- [`data/`](data/): Lightweight input datasets (excluding GIS once), starter templates, and sample configurations used by the tutorials and processing pipelines.
- [`documentation/`](documentation/): Theoretical references and step-by-step guides to develope LEAFs, including [SOC documentation](documentation/SOC_Documentation.md), and [Soil Erosion guide](documentation/Soil_Erosion_Documentation.md),as well as data harmonization notes.
- [`examples/`](examples/): Executable notebooks and scripts that demonstrate how to prepare inputs and compute LEAF values.
- [`src/`](src/): Reusable Python modules that implement the harmonization, calculation, and export utilities used across the examples. Shared lookups that power crop and climate calculations are lazily loaded through [`src/sbtn_leaf/data_loader.py`](src/sbtn_leaf/data_loader.py) so that downstream modules can request cached tables or inject test fixtures.
- [`LEAFs/`](LEAFs/): Output directory with precomputed LEAFs tables and maps. New LEAFs generated are stored after here after running the processing workflows. See [LEAFs/README.md](LEAFs/README.md) for guidance on organizing and versioning run artifacts.

## Data and Documentation Resources
### Lightweight datasets
Input samples and reduced-size datasets are stored in [`data/`](data/). Each subfolder includes a `README` or inline metadata to explain processing steps to substitute with new data when generating new or modifying existing LEAFs.

### Theory and methodology documentation
Detailed theory, assumptions, and methodology walkthroughs for each LEAF category live in [`documentation/`](documentation/).

### Example notebooks and scripts
The [`examples/`](examples/) folder contains Jupyter notebooks and Python scripts that mirror the workflows described in the documentation. Each example references the lightweight datasets in [`data/`](data/) and saves results into [`LEAFs/`](LEAFs/).

### LEAF outputs
Precomputed LEAFs for SOC, Soil Erosion and Acidification are available in datatables inside [`LEAFs/`](LEAFs/) folder. 
New LEAFs outputs should be stored here inlcuidng category name, commodity, and  date to keep multiple iterations side by side. When creating new LEAFs, consult [LEAFs/README.md](LEAFs/README.md) for recommendations on which outputs to version (for example lightweight CSV summaries such as [`sample_soc_summary.csv`](LEAFs/sample_soc_summary.csv)) and how to keep large rasters or proprietary
datasets out of Git.

## Quick Start
1. Install the project dependencies:
   ```bash
   pip install -e .
   ```
2. Explore the theory and data preparation guidance in [`documentation/`](documentation/).
3. Open an example notebook from [`examples/`](examples/) and update the configuration cells to point at the desired inputs inside [`data/`](data/).
4. Execute the notebook or script. Generated LEAFs will appear under [`LEAFs/`](LEAFs/).

## Running the Examples
To work through the notebooks:
1. Clone repository into your local drive and launch Jupyter Lab or Notebook from the project root:
   ```bash
   jupyter lab
   ```
2. Open one of the example notebooks in [`examples/`](examples/) and follow the embedded instructions. Each example demonstrates how to load data, preprocess it, and generate LEAF using utilities from [`src/`](src/), and export results.
3. After the run completes, inspect the generated outputs in [`LEAFs/`](LEAFs/) to validate results or share them with your team.

## Versioning
The package version is managed from the single source declared in [`pyproject.toml`](pyproject.toml). The `sbtn_leaf`
package exposes this value at runtime by reading the installed distribution metadata and, during local development,
falling back to the `pyproject.toml` entry. Update the version in `pyproject.toml` when preparing a release to ensure
that both the package metadata and the `sbtn_leaf.__version__` attribute report the same number.

## Contributing
We welcome improvements to the data pipelines, documentation, examples, and new LEAFs. To contribute:
1. Fork the repository and create a feature branch.
2. Make your changes, ensuring new notebooks or datasets are placed in the appropriate top-level folders.
3. Add or update documentation in [`documentation/`](documentation/) and reference new assets from the relevant sections above.
4. Run the provided examples or tests to confirm your changes, then open a pull request describing the updates and any expected LEAF outputs that reviewers should verify in [`LEAFs/`](LEAFs/).
