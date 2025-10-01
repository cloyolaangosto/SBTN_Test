"""
sbtn_leaf: Environmental modeling utilities for SBTN-LEAF project.

This package provides tools for:
- Running environmental and land-use models.
- Working with geospatial raster/vector data.
- Managing data preprocessing and result analysis.

Typical usage example:
    import sbtn_leaf as sl

    result = sl.run_model(config="config.yml")
    sl.plot_results(result)
"""

# Re-export key functions/classes from submodules
# from .core import run_model, ModelConfig
# from .helpers import (
#     load_raster,
#     save_raster,
#     calculate_statistics,
# )

# # Optional: define what gets imported by `from sbtn_leaf import *`
# __all__ = [
#     "run_model",
#     "ModelConfig",
#     "load_raster",
#     "save_raster",
#     "calculate_statistics",
# ]

# Package metadata
from importlib import metadata as _metadata
from pathlib import Path


try:
    __version__ = _metadata.version("sbtn_leaf")
except _metadata.PackageNotFoundError:
    try:  # Python 3.11+
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python <3.11
        tomllib = None

    _pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    if tomllib is not None and _pyproject.exists():
        with _pyproject.open("rb") as _fp:
            __version__ = tomllib.load(_fp)["project"]["version"]
    else:
        __version__ = "0.0.0.dev1"

__author__ = "Cristóbal Loyola"
