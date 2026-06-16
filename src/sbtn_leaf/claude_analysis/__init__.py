"""Claude-generated spatial-aggregation analyses of SBTN-Land LEAFs.

This subpackage groups the aggregation-comparison work:

* :mod:`~sbtn_leaf.claude_analysis.se_aggregation_analysis` -- the original
  soil-erosion-only comparison (country / subcountry / ecoregion).
* :mod:`~sbtn_leaf.claude_analysis.indicator_aggregation` -- a generic engine that
  harmonises any indicator onto one schema and provides the shared statistics and
  figures.
* :mod:`~sbtn_leaf.claude_analysis.indicators` -- the :class:`IndicatorConfig`
  definitions for SOC, soil erosion and terrestrial acidification.
* :mod:`~sbtn_leaf.claude_analysis.cross_indicator` -- the cross-indicator
  comparison (how ecoregions reframe averages differently than country / admin-1
  units) and the one-shot ``run_cross_indicator_analysis`` pipeline.
* :mod:`~sbtn_leaf.claude_analysis.manuscript_support` -- extra statistics and
  figures backing specific claims in the LEAF manuscript (biome significance,
  within-region SD, SOC<->erosion alignment).

The narrative notebooks live under ``paper/claude_analysis/``.
"""

from __future__ import annotations

from sbtn_leaf.claude_analysis import (
    cross_indicator,
    indicator_aggregation,
    indicators,
    manuscript_support,
    se_aggregation_analysis,
)
from sbtn_leaf.claude_analysis.indicator_aggregation import IndicatorConfig
from sbtn_leaf.claude_analysis.indicators import (
    ACIDIFICATION,
    INDICATORS,
    SOC,
    SOIL_EROSION,
    get_indicator,
)
from sbtn_leaf.claude_analysis.cross_indicator import run_cross_indicator_analysis
from sbtn_leaf.claude_analysis.manuscript_support import run_manuscript_support

__all__ = [
    "se_aggregation_analysis",
    "indicator_aggregation",
    "indicators",
    "cross_indicator",
    "manuscript_support",
    "IndicatorConfig",
    "INDICATORS",
    "SOC",
    "SOIL_EROSION",
    "ACIDIFICATION",
    "get_indicator",
    "run_cross_indicator_analysis",
    "run_manuscript_support",
]
