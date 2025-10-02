"""Utility helpers for lazily loading shared crop and climate tables.

The PET and crop calculation modules both rely on the crop coefficient
(`K_Crop_Data.csv`), absolute day (`AbsoluteDayTable.csv`), and climate
lookup tables. Historically each module read these files independently at
import time which duplicated IO and made it awkward to inject custom test
data.

This module centralises the logic for locating those files, loads them on
demand, and caches the parsed ``polars`` objects.  Callers can request a
fresh clone of each table whenever needed or pass the tables around
explicitly for tests.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import calendar
from typing import Dict, Iterable, Mapping, Tuple

import geopandas as gpd
import polars as pl

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"


def _ensure_exists(path: Path, *, description: str) -> Path:
    """Ensure ``path`` exists before attempting to read it."""

    if not path.exists():
        raise FileNotFoundError(f"Expected {description} at '{path}'.")
    return path


@lru_cache(maxsize=None)
def _load_crop_coefficients_table() -> pl.DataFrame:
    """Read the crop coefficient CSV from disk."""

    path = _ensure_exists(
        DATA_DIR / "crops" / "K_Crop_Data.csv",
        description="crop coefficient lookup table",
    )

    return pl.read_csv(path)


def get_crop_coefficients_table() -> pl.DataFrame:
    """Return a cached copy of the crop coefficient table."""

    return _load_crop_coefficients_table().clone()


@lru_cache(maxsize=None)
def _load_absolute_day_table() -> pl.DataFrame:
    """Read the absolute day lookup table from disk."""

    path = _ensure_exists(
        DATA_DIR / "crops" / "AbsoluteDayTable.csv",
        description="absolute day lookup table",
    )

    return pl.read_csv(path)


def get_absolute_day_table() -> pl.DataFrame:
    """Return a cached copy of the absolute day lookup table."""

    return _load_absolute_day_table().clone()


@lru_cache(maxsize=None)
def _build_days_in_month_table(year: int = 2023) -> pl.DataFrame:
    """Construct a Polars table containing the number of days per month."""

    return pl.DataFrame(
        {
            "Month": list(range(1, 13)),
            "Days_in_Month": [calendar.monthrange(year, month)[1] for month in range(1, 13)],
        }
    )


def get_days_in_month_table(year: int = 2023) -> pl.DataFrame:
    """Return a cached copy of the days-in-month table for ``year``."""

    return _build_days_in_month_table(year).clone()


@lru_cache(maxsize=None)
def _load_crop_naming_index_table() -> pl.DataFrame:
    """Read the crop naming index CSV from disk."""

    path = _ensure_exists(
        DATA_DIR / "crops" / "crop_naming_index.csv",
        description="crop naming index table",
    )

    return pl.read_csv(path)


def get_crop_naming_index_table() -> pl.DataFrame:
    """Return a cached copy of the crop naming index table."""

    return _load_crop_naming_index_table().clone()


@lru_cache(maxsize=None)
def _load_fao_statistics_table() -> pl.DataFrame:
    """Read the FAO production statistics CSV from disk."""

    path = _ensure_exists(
        DATA_DIR / "crops" / "Production_Crops_Livestock_E_All_Data.csv",
        description="FAO production statistics table",
    )

    return pl.read_csv(path)


def get_fao_statistics_table() -> pl.DataFrame:
    """Return a cached copy of the FAO production statistics table."""

    return _load_fao_statistics_table().clone()


@lru_cache(maxsize=None)
def _load_fao_crop_yields_table() -> pl.DataFrame:
    """Read the FAO crop yields CSV from disk."""

    path = _ensure_exists(
        DATA_DIR / "crops" / "fao_crop_yields_1423.csv",
        description="FAO crop yields table",
    )

    return pl.read_csv(path, separator=";")


def get_fao_crop_yields_table() -> pl.DataFrame:
    """Return a cached copy of the FAO crop yields table."""

    return _load_fao_crop_yields_table().clone()


@lru_cache(maxsize=None)
def _load_country_boundaries() -> gpd.GeoDataFrame:
    """Read the country boundary shapefile from disk."""

    path = _ensure_exists(
        DATA_DIR / "CountryLayers" / "Country_Level0" / "g2015_2014_0.shp",
        description="country boundary shapefile",
    )

    return gpd.read_file(path)


def get_country_boundaries() -> gpd.GeoDataFrame:
    """Return a cached copy of the country boundary shapefile."""

    return _load_country_boundaries().copy()


@lru_cache(maxsize=None)
def _load_ecoregions_shapefile() -> gpd.GeoDataFrame:
    """Read the ecoregions shapefile from disk."""

    path = _ensure_exists(
        DATA_DIR / "Ecoregions2017" / "Ecoregions2017.shp",
        description="ecoregions shapefile",
    )

    return gpd.read_file(path)


def get_ecoregions_shapefile() -> gpd.GeoDataFrame:
    """Return a cached copy of the ecoregions shapefile."""

    return _load_ecoregions_shapefile().copy()


@lru_cache(maxsize=None)
def _load_crop_ag_residue_table() -> pl.DataFrame:
    """Read the crop above-ground residue Excel sheet from disk."""

    path = _ensure_exists(
        DATA_DIR / "crops" / "crop_residue_data.xlsx",
        description="crop residue workbook",
    )

    return pl.read_excel(path, sheet_name="crop_ABG_Res")


def get_crop_ag_residue_table() -> pl.DataFrame:
    """Return a cached copy of the crop above-ground residue table."""

    return _load_crop_ag_residue_table().clone()


@lru_cache(maxsize=None)
def _load_crop_residue_ratio_table() -> pl.DataFrame:
    """Read the crop residue ratio Excel sheet from disk."""

    path = _ensure_exists(
        DATA_DIR / "crops" / "crop_residue_data.xlsx",
        description="crop residue workbook",
    )

    return pl.read_excel(path, sheet_name="crop_res_ratios")


def get_crop_residue_ratio_table() -> pl.DataFrame:
    """Return a cached copy of the crop residue ratio table."""

    return _load_crop_residue_ratio_table().clone()


THERMAL_CLIMATE_ROWS = [
    (1, "Tropics, lowland", "Tropics"),
    (2, "Tropics, highland", "Tropics"),
    (3, "Subtropics, summer rainfall", "Subtropics summer rainfall"),
    (4, "Subtropics, winter rainfall", "Subtropics winter rainfall"),
    (5, "Subtropics, low rainfall", "Subtropics winter rainfall"),
    (6, "Temperate, oceanic", "Oceanic temperate"),
    (7, "Temperate, sub-continental", "Sub-continental temperate and continental temperate"),
    (8, "Temperate, continental", "Sub-continental temperate and continental temperate"),
    (9, "Boreal, oceanic", "Sub-continental boreal, continental boreal and polar/arctic"),
    (10, "Boreal, sub-continental", "Sub-continental boreal, continental boreal and polar/arctic"),
    (11, "Boreal, continental", "Sub-continental boreal, continental boreal and polar/arctic"),
    (12, "Arctic", "Sub-continental boreal, continental boreal and polar/arctic"),
]


@lru_cache(maxsize=None)
def _build_thermal_climate_tables() -> Tuple[pl.DataFrame, Dict[int, str], Dict[str, Tuple[int, ...]]]:
    """Create the thermal climate lookup table and related dictionaries."""

    table = pl.DataFrame(THERMAL_CLIMATE_ROWS, schema=["id", "TC_Name", "TC_Group"])

    climate_zone_lookup: Dict[int, str] = dict(zip(table["id"].to_list(), table["TC_Group"].to_list()))

    zone_ids: Dict[str, list[int]] = {}
    for zone_id, group in climate_zone_lookup.items():
        zone_ids.setdefault(group, []).append(zone_id)

    zone_ids_by_group = {group: tuple(ids) for group, ids in zone_ids.items()}

    return table, climate_zone_lookup, zone_ids_by_group


def get_thermal_climate_tables(
    *,
    include_lookup: bool = True,
    include_zone_ids: bool = True,
) -> Tuple[pl.DataFrame, Mapping[int, str], Mapping[str, Iterable[int]]]:
    """Return cached climate tables and associated helper mappings."""

    table, lookup, zone_ids = _build_thermal_climate_tables()

    output_table = table.clone()
    output_lookup: Mapping[int, str]
    output_zone_ids: Mapping[str, Iterable[int]]

    output_lookup = dict(lookup) if include_lookup else {}
    output_zone_ids = dict(zone_ids) if include_zone_ids else {}

    return output_table, output_lookup, output_zone_ids
