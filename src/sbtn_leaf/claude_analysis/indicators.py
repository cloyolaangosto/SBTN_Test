"""Concrete :class:`IndicatorConfig` definitions for the three LEAF indicators.

Each indicator is reduced to the same harmonised long schema
(:data:`~sbtn_leaf.claude_analysis.indicator_aggregation.STD_COLS`) so the generic
engine can treat them identically:

* **SOC** -- estimated 2030 soil-organic-carbon stock, ``t SOC/ha`` (higher is
  better).  Source: the merged ``SOC_2030_<level>_v1.0.csv`` tables (110
  land-use x management flows).  The ecoregion ``v1.0`` table lacks biome / realm
  columns, so they are joined from the SOC ecoregion crop / forest tables.
* **Soil erosion** -- RUSLE soil loss, ``t soil/ha/yr`` (higher is worse).
  Delegated to :mod:`sbtn_leaf.claude_analysis.se_aggregation_analysis`, whose
  harmonised frame already matches the shared schema.
* **Terrestrial acidification** -- accumulated exceedance characterisation
  factor, ``kg SO2-eq./kg`` emitted (higher is worse).  Three acidifying gases
  (NOx, NH3, SO2); the CF is defined region-wide rather than per commodity.

Canonical-flow keys are shared with the soil-erosion convention
(``Commodity|water|residue|tillage`` / ``BRDC_*`` / ``Grassland``) so the same
commodity lines up across SOC and soil erosion in the cross-indicator view.
"""

from __future__ import annotations

from collections import OrderedDict
from functools import lru_cache

import pandas as pd

from sbtn_leaf.paths import leaf_path
from sbtn_leaf.claude_analysis import se_aggregation_analysis as se
from sbtn_leaf.claude_analysis.indicator_aggregation import (
    LEVELS,
    IndicatorConfig,
    load_long_level,
)

__all__ = ["SOC", "SOIL_EROSION", "ACIDIFICATION", "INDICATORS", "get_indicator"]

_SOC_SUFFIX = "_2030y_SOC"


# --------------------------------------------------------------------------- #
# SOC
# --------------------------------------------------------------------------- #

#: Source ``v1.0`` SOC tables (merged leaf / leaf_median / leaf_std), by level.
_SOC_FILES = {
    "country": "SOC_2030_country_v1.0.csv",
    "subcountry": "SOC_2030_subcountry_v1.0.csv",
    "ecoregion": "SOC_2030_ecoregions_v1.0.csv",
}

#: Ecoregion tables that *do* carry the biome / realm grouping, used to backfill
#: the metadata the ``v1.0`` ecoregion table omits.
_SOC_ECO_META_FILES = (
    "SOC_2030_ecoregions_crops_clipped.csv",
    "SOC_2030_ecoregions_forest_grass.csv",
)


def canonical_flow_soc(raw: str) -> str:
    """Map a SOC ``flow_name`` onto the shared canonical key.

    Strips the ``_2030y_SOC`` suffix, normalises spaces to underscores (so
    ``"Oil palm"`` / ``"BRDC_Boreal dry"`` line up with the soil-erosion
    spellings), handles the ``natural_grassland_*`` herd-management variants, then
    delegates crops / forest to the soil-erosion canonicaliser.  ``Soybean`` is
    aliased to ``Soybeans`` to match the soil-erosion vocabulary.
    """

    s = str(raw)
    if s.endswith(_SOC_SUFFIX):
        s = s[: -len(_SOC_SUFFIX)]
    s = s.replace(" ", "_")

    if s.startswith("natural_grassland"):
        mgmt = s[len("natural_grassland_") :].strip("_")
        return "Grassland" if mgmt == "cattle_avg" else f"Grassland|{mgmt}"

    key = se.canonical_flow(s)
    if key.startswith("Soybean|"):
        key = "Soybeans|" + key.split("|", 1)[1]
    return key


@lru_cache(maxsize=1)
def _soc_eco_meta() -> pd.DataFrame:
    frames = [
        pd.read_csv(leaf_path("SOC", f))[["ECO_ID", "ECO_NAME", "BIOME_NAME", "REALM"]]
        for f in _SOC_ECO_META_FILES
    ]
    return pd.concat(frames, ignore_index=True).drop_duplicates("ECO_ID")


@lru_cache(maxsize=1)
def _load_soc() -> pd.DataFrame:
    eco_meta = _soc_eco_meta()
    parts = []
    for level in LEVELS:
        parts.append(
            load_long_level(
                leaf_path("SOC", _SOC_FILES[level]),
                level,
                pivot_col="variable",
                value_map={"leaf": "leaf", "leaf_median": "leaf_median", "leaf_std": "leaf_std"},
                canonical=canonical_flow_soc,
                eco_meta=eco_meta if level == "ecoregion" else None,
            )
        )
    return pd.concat(parts, ignore_index=True)


#: Representative focal flows, sharing the soil-erosion keys so the same commodity
#: lines up across indicators.  Labels come from the soil-erosion labeller.
_SOC_FOCAL = OrderedDict(
    (k, se.flow_label(k))
    for k in (
        "Wheat|rf|roff|ct",
        "Maize|rf|roff|ct",
        "Soybeans|rf|na|ct",
        "Oil_palm|rf|na|ct",
        "Coffee|rf|na|ct",
        "Sugarcane|rf|na|ct",
        "Cotton|rf|na|ct",
        "Grassland",
        "BRDC_Tropical",
    )
)

SOC = IndicatorConfig(
    key="soc",
    name="SOC stock",
    unit="t SOC/ha",
    flow_kind="commodity",
    higher_is_better=True,
    focal_flows=_SOC_FOCAL,
    loader=_load_soc,
    labeller=se.flow_label,
)


# --------------------------------------------------------------------------- #
# Soil erosion (delegate loading to the soil-erosion module)
# --------------------------------------------------------------------------- #


@lru_cache(maxsize=1)
def _load_se() -> pd.DataFrame:
    # se.load_harmonized already returns the shared schema (drop_na=False keeps the
    # no-data cells so coverage is comparable).
    return se.load_harmonized(drop_na=False)


SOIL_EROSION = IndicatorConfig(
    key="soil_erosion",
    name="Soil erosion",
    unit="t soil/ha/yr",
    flow_kind="commodity",
    higher_is_better=False,
    focal_flows=se.FOCAL_FLOWS,
    loader=_load_se,
    labeller=se.flow_label,
)


# --------------------------------------------------------------------------- #
# Terrestrial acidification
# --------------------------------------------------------------------------- #

_ACID_FILES = {
    "country": "acidification_country.csv",
    "subcountry": "acidification_subcountry.csv",
    "ecoregion": "acidification_ecoregion.csv",
}

#: The three acidifying gases (the "flows" of this indicator).
_ACID_LABELS = {
    "acid_nox": "NOₓ (acidification)",
    "acid_nh3": "NH₃ (acidification)",
    "acid_so2": "SO₂ (acidification)",
}


def canonical_flow_acid(raw: str) -> str:
    """Acidification flows are already canonical gas codes (``acid_nox`` ...)."""

    return str(raw)


def label_acid(flow: str) -> str:
    return _ACID_LABELS.get(flow, flow)


@lru_cache(maxsize=1)
def _load_acid() -> pd.DataFrame:
    parts = []
    for level in LEVELS:
        parts.append(
            load_long_level(
                leaf_path("acidification", _ACID_FILES[level]),
                level,
                pivot_col="metric",
                value_map={"cf_mean": "leaf", "cf_median": "leaf_median", "cf_std": "leaf_std"},
                canonical=canonical_flow_acid,
            )
        )
    return pd.concat(parts, ignore_index=True)


_ACID_FOCAL = OrderedDict((k, label_acid(k)) for k in ("acid_nox", "acid_nh3", "acid_so2"))

ACIDIFICATION = IndicatorConfig(
    key="acidification",
    name="Acidification CF",
    unit="kg SO2-eq./kg",
    flow_kind="pollutant",
    higher_is_better=False,
    focal_flows=_ACID_FOCAL,
    loader=_load_acid,
    labeller=label_acid,
)


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #

#: All indicators, ordered for display.
INDICATORS = OrderedDict(
    [
        ("soc", SOC),
        ("soil_erosion", SOIL_EROSION),
        ("acidification", ACIDIFICATION),
    ]
)


def get_indicator(key: str) -> IndicatorConfig:
    """Look up an :class:`IndicatorConfig` by its short ``key``."""

    try:
        return INDICATORS[key]
    except KeyError:
        raise ValueError(f"Unknown indicator {key!r}; expected one of {list(INDICATORS)}.") from None
