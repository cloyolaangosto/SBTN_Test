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
* :mod:`~sbtn_leaf.claude_analysis.practice_change` -- multi-indicator practice-switch
  co-benefits (SOC gained + erosion avoided) with maps of where to focus, plus the
  :mod:`~sbtn_leaf.claude_analysis.geo` geometry loaders for the choropleths.

The narrative notebooks live under ``paper/claude_analysis/``.
"""

from __future__ import annotations

from sbtn_leaf.claude_analysis import (
    cross_indicator,
    geo,
    indicator_aggregation,
    indicators,
    manuscript_support,
    practice_change,
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
from sbtn_leaf.claude_analysis.practice_change import run_practice_change_analysis

__all__ = [
    "se_aggregation_analysis",
    "indicator_aggregation",
    "indicators",
    "cross_indicator",
    "manuscript_support",
    "practice_change",
    "geo",
    "IndicatorConfig",
    "INDICATORS",
    "SOC",
    "SOIL_EROSION",
    "ACIDIFICATION",
    "get_indicator",
    "run_cross_indicator_analysis",
    "run_manuscript_support",
    "run_practice_change_analysis",
]
