"""Geometry loaders for choropleth maps of country / ecoregion regions.

No boundary polygons ship with this checkout (they are DVC-tracked and usually
absent), so this module resolves geometry through a fallback chain:

1. the DVC shapefiles via :mod:`sbtn_leaf.data_loader`, when present;
2. otherwise a one-time download of public geometry — Natural Earth admin-0 for
   countries, RESOLVE/WWF *Ecoregions 2017* for ecoregions — cached under
   ``~/.cache/sbtn_leaf_geo`` (override with ``$SBTN_GEO_CACHE``);
3. otherwise ``None`` with a warning, so every map degrades gracefully to a skip
   (mirroring :func:`sbtn_leaf.claude_analysis.se_aggregation_analysis.load_boundaries`).

Country values join by name (FAO ADM0 → Natural Earth, via :func:`country_key`
normalisation plus a small alias table); ecoregion values join by ``ECO_ID``.
"""

from __future__ import annotations

import os
import re
import warnings
from pathlib import Path

import pandas as pd

__all__ = [
    "geo_cache_dir",
    "country_geometry",
    "ecoregion_geometry",
    "country_key",
    "choropleth",
]

_NE_COUNTRIES_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_110m_admin_0_countries.geojson"
)
_ECOREGIONS_URL = "https://storage.googleapis.com/teow2016/Ecoregions2017.zip"

#: FAO ADM0 spellings that do not normalise onto the Natural Earth name; keyed on
#: the normalised FAO name, mapped to the normalised Natural Earth name.
_COUNTRY_ALIAS = {
    "russian federation": "russia",
    "lao peoples democratic republic": "laos",
    "dem peoples rep of korea": "north korea",
    "republic of korea": "south korea",
    "moldova republic of": "moldova",
    "republic of moldova": "moldova",
    "united republic of tanzania": "tanzania",
    "syrian arab republic": "syria",
    "viet nam": "vietnam",
    "brunei darussalam": "brunei",
    "libyan arab jamahiriya": "libya",
    "czech republic": "czechia",
    "the former yugoslav republic of macedonia": "north macedonia",
    "republic of north macedonia": "north macedonia",
    "united kingdom of great britain and northern ireland": "united kingdom",
    "congo": "republic of the congo",
    "cote divoire": "ivory coast",
    "swaziland": "eswatini",
    "burma": "myanmar",
    "united states of america": "united states of america",
    "bolivia plurinational state of": "bolivia",
    "venezuela bolivarian republic of": "venezuela",
    "iran islamic republic of": "iran",
    "tanzania united republic of": "tanzania",
    "the bahamas": "bahamas",
    "gambia the": "gambia",
}


def geo_cache_dir() -> Path:
    """Directory for cached downloaded geometry (``$SBTN_GEO_CACHE`` or ``~/.cache``)."""

    d = Path(os.environ.get("SBTN_GEO_CACHE", Path.home() / ".cache" / "sbtn_leaf_geo"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download(url: str, dest: Path, *, timeout: int = 120) -> Path:
    import urllib.request

    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=timeout) as r, open(tmp, "wb") as f:
        f.write(r.read())
    tmp.replace(dest)
    return dest


def country_key(name: str) -> str:
    """Normalise a country name for joining (drop parentheticals/punctuation, alias)."""

    s = re.sub(r"\(.*?\)", "", str(name))
    s = s.lower().replace("&", "and")
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return _COUNTRY_ALIAS.get(s, s)


# --------------------------------------------------------------------------- #
# Country geometry
# --------------------------------------------------------------------------- #


def country_geometry(*, allow_download: bool = True):
    """Country polygons with a normalised ``key`` column, or ``None`` if unavailable."""

    try:
        import geopandas as gpd
    except Exception as exc:  # pragma: no cover
        warnings.warn(f"geopandas unavailable ({exc}); country maps skipped.")
        return None

    # 1) DVC shapefile, if present.
    try:
        from sbtn_leaf.data_loader import get_country_boundaries

        gdf = get_country_boundaries()
        name_col = next((c for c in ("ADM0_NAME", "NAME", "ADMIN") if c in gdf.columns), None)
        if name_col is not None:
            gdf = gdf[[name_col, "geometry"]].rename(columns={name_col: "name"})
            gdf["key"] = gdf["name"].map(country_key)
            return gdf
    except Exception:
        pass

    # 2) Cached Natural Earth file (always reusable), else download it.
    cache = geo_cache_dir() / "ne_110m_admin_0_countries.geojson"
    if not (cache.exists() and cache.stat().st_size > 0):
        if not allow_download:
            warnings.warn("Country geometry not local and downloads disabled; skipping country maps.")
            return None
        try:
            _download(_NE_COUNTRIES_URL, cache)
        except Exception as exc:  # pragma: no cover - network dependent
            warnings.warn(f"Country geometry download failed ({type(exc).__name__}: {exc}); skipping.")
            return None
    try:
        gdf = gpd.read_file(cache)
        name_col = next((c for c in ("NAME", "ADMIN", "NAME_LONG", "SOVEREIGNT") if c in gdf.columns), None)
        gdf = gdf[[name_col, "geometry"]].rename(columns={name_col: "name"})
        gdf["key"] = gdf["name"].map(country_key)
        return gdf
    except Exception as exc:  # pragma: no cover
        warnings.warn(f"Country geometry read failed ({type(exc).__name__}: {exc}); skipping.")
        return None


# --------------------------------------------------------------------------- #
# Ecoregion geometry
# --------------------------------------------------------------------------- #


def ecoregion_geometry(*, allow_download: bool = True, simplify_tol: float = 0.08):
    """WWF/RESOLVE Ecoregions-2017 polygons keyed by ``ECO_ID``, or ``None``.

    The 149 MB source is downloaded once, simplified and cached as a lightweight
    GeoPackage so subsequent loads are fast.
    """

    try:
        import geopandas as gpd
    except Exception as exc:  # pragma: no cover
        warnings.warn(f"geopandas unavailable ({exc}); ecoregion maps skipped.")
        return None

    cache = geo_cache_dir() / "ecoregions2017_simplified.gpkg"
    if cache.exists():
        try:
            return gpd.read_file(cache)
        except Exception:
            pass

    # 1) DVC shapefile, if present.
    gdf = None
    try:
        from sbtn_leaf.data_loader import get_ecoregions_shapefile

        gdf = get_ecoregions_shapefile()
    except Exception:
        gdf = None

    # 2) Download the RESOLVE archive.
    if gdf is None:
        if not allow_download:
            warnings.warn("Ecoregion geometry not local and downloads disabled; skipping ecoregion maps.")
            return None
        try:
            zpath = _download(_ECOREGIONS_URL, geo_cache_dir() / "Ecoregions2017.zip", timeout=600)
            gdf = gpd.read_file(f"zip://{zpath}")
        except Exception as exc:  # pragma: no cover - network dependent
            warnings.warn(f"Ecoregion geometry download failed ({type(exc).__name__}: {exc}); skipping.")
            return None

    if "ECO_ID" not in gdf.columns:
        warnings.warn("Ecoregion geometry lacks ECO_ID; skipping ecoregion maps.")
        return None

    gdf = gdf[["ECO_ID", "geometry"]].copy()
    gdf["ECO_ID"] = pd.to_numeric(gdf["ECO_ID"], errors="coerce")
    gdf = gdf[gdf["ECO_ID"].notna()]
    gdf["ECO_ID"] = gdf["ECO_ID"].astype(int)
    try:
        gdf["geometry"] = gdf.geometry.simplify(simplify_tol)
        gdf.to_file(cache, driver="GPKG")
    except Exception:  # pragma: no cover - simplification/caching is best-effort
        pass
    return gdf


# --------------------------------------------------------------------------- #
# Generic choropleth
# --------------------------------------------------------------------------- #


def choropleth(
    gdf,
    value_col: str = "value",
    *,
    ax=None,
    cmap: str = "viridis",
    norm=None,
    title: str = "",
    cbar_label: str = "",
    basemap=None,
    missing_color: str = "whitesmoke",
):
    """Plot one ``value_col`` choropleth; returns ``(fig, ax)``."""

    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 5.5))
    else:
        fig = ax.figure
    if basemap is not None:
        basemap.plot(ax=ax, color="0.92", edgecolor="white", linewidth=0.2)
    gdf.plot(
        ax=ax,
        column=value_col,
        cmap=cmap,
        norm=norm,
        legend=False,
        missing_kwds={"color": missing_color},
    )
    ax.set_title(title)
    ax.set_axis_off()
    if norm is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, shrink=0.6, label=cbar_label)
    return fig, ax
