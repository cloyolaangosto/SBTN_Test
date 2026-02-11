# MAP FUNCTIONS #

# MODULES
# import polars as pl
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm, Normalize, LinearSegmentedColormap
import matplotlib.cm as cm
import plotly.express as px
import rioxarray
import rasterio
from rasterio.enums import Resampling
from rasterio.coords import BoundingBox
from pathlib import Path
from typing import Dict, Optional, Tuple, Union, Any, Sequence

# World shapefile (loaded lazily)
from sbtn_leaf.paths import data_path

_WORLD_MAP_PATHS = {
    "lr": data_path("world_maps", "low_res", "ne_110m_admin_0_countries.shp"),
    "hr": data_path("world_maps", "high_res", "ne_10m_admin_0_countries.shp"),
}

_world_map_cache: Dict[str, gpd.GeoDataFrame] = {}


def _safe_read_file(path: str) -> gpd.GeoDataFrame:
    try:
        return gpd.read_file(path)
    except Exception:  # pragma: no cover - optional dependency for plotting
        return gpd.GeoDataFrame()


def _get_world_map(resolution: str = "hr") -> gpd.GeoDataFrame:
    """Return a cached world basemap for the requested ``resolution``.

    Parameters
    ----------
    resolution:
        Key identifying the desired basemap resolution. ``"hr"`` (high
        resolution) remains the default to keep existing behaviour. ``"lr"``
        (low resolution) is also available out of the box.

    Notes
    -----
    Advanced callers can pre-load their own basemap and register it with
    :func:`_register_world_map` to avoid redundant disk reads when plotting
    repeatedly.
    """

    resolution_key = resolution.lower()
    try:
        shapefile_path = _WORLD_MAP_PATHS[resolution_key]
    except KeyError as exc:  # pragma: no cover - defensive guard
        raise ValueError(f"Unknown world map resolution: {resolution!r}") from exc

    if resolution_key not in _world_map_cache:
        _world_map_cache[resolution_key] = _safe_read_file(shapefile_path)

    return _world_map_cache[resolution_key]


def _register_world_map(resolution: str, basemap: gpd.GeoDataFrame) -> None:
    """Register ``basemap`` for ``resolution`` to reuse without disk I/O.

    This is primarily intended for advanced and performance-sensitive callers
    that want to supply a custom GeoDataFrame (perhaps already filtered) before
    invoking the plotting utilities in this module.
    """

    _world_map_cache[resolution.lower()] = basemap


# FUNCTIONS
def _preprocess_raster_data_eliminate_nodata(
    raster_data: np.ndarray,
    nodata_value: Optional[float] = None,
    eliminate_zeros: bool = True
) -> np.ndarray:
    """
    Replace nodata (e.g. -32000) and ±inf with np.nan.
    """
    out = raster_data.astype("float32", copy=True)

    if nodata_value is not None:
        out[out == nodata_value] = np.nan

    if eliminate_zeros:
        out[out == 0] = np.nan

    # any inf / -inf -> NaN
    out[~np.isfinite(out)] = np.nan

    return out


def _preprocess_raster_data_percentiles(
    raster_data: np.ndarray,
    nodata_value: Optional[float] = None,
    p_min: Union[int, float] = 1,
    p_max: Union[int, float] = 99,
    hard_min: Optional[float] = None,
    hard_max: Optional[float] = None,
    eliminate_zeros: Optional[bool] = None
) -> np.ndarray:
    """
    1. Convert nodata and inf to NaN.
    2. Clip to percentiles to reduce extreme spikes.
    3. Optionally clip to explicit physical bounds (hard_min / hard_max).
    """

    # Step 1: remove nodata, inf
    raster_data = _preprocess_raster_data_eliminate_nodata(raster_data, nodata_value, eliminate_zeros)

    # Step 2: percentile clipping on finite pixels
    finite_mask = np.isfinite(raster_data)
    finite_vals = raster_data[finite_mask]
    if finite_vals.size > 0:
        lower_p, upper_p = np.percentile(finite_vals, [p_min, p_max])
        raster_data = np.clip(raster_data, lower_p, upper_p)

    # Step 3: optional hard clip (physical plausibility)
    if hard_min is not None:
        raster_data = np.where(raster_data < hard_min, hard_min, raster_data)
    if hard_max is not None:
        raster_data = np.where(raster_data > hard_max, hard_max, raster_data)

    return raster_data


def _preprocess_raster_data_eliminate_low_values(raster_data, nodata_value=None, threshold = None):
    """
    Preprocess raster data to handle invalid and extreme values.
    
    Parameters:
        raster_data (ndarray): The raster data array.
        nodata_value (float): Value representing no-data in the raster.
    
    Returns:
        ndarray: Cleaned raster data array.
    """
    # Eliminate no_data and Inf values
    raster_data = _preprocess_raster_data_eliminate_nodata(raster_data, nodata_value)

    # Determine minimum value for as threshold if none is given
    if threshold is None:
        threshold = np.min(raster_data)
    
    # Elimiante minimum values
    finite_data = raster_data[np.isfinite(raster_data)]
    raster_data = np.where(finite_data <= threshold, np.nan, finite_data)
    
    return raster_data


def _prepare_raster_plot_input(
    source: Union[str, Path, xr.DataArray],
    *,
    raster_band: int = 1,
    band: Optional[int] = None,
    perc_cutoff: Optional[float] = None,
    p_min: Optional[float] = None,
    p_max: Optional[float] = None,
    hard_min: Optional[float] = None,
    hard_max: Optional[float] = None,
    eliminate_zeros: bool = False,
) -> Tuple[np.ndarray, BoundingBox, Optional[float], Any, Dict[str, Optional[float]]]:
    """Normalise raster-like inputs for plotting.

    Parameters
    ----------
    source:
        Either a path to a raster on disk or an :class:`xarray.DataArray`.
    raster_band:
        Band to read when ``source`` is a file path. 1-indexed to match
        :mod:`rasterio` conventions.
    band:
        Band to select when ``source`` is a DataArray with a ``"band"``
        dimension. Ignored otherwise.
    perc_cutoff:
        Symmetric percentile cutoff used when ``p_min``/``p_max`` are not
        provided. Mirrors the legacy ``perc_cutoff``/``alpha`` parameters.
    p_min, p_max:
        Explicit percentile cut-offs. When omitted they default to
        ``(perc_cutoff, 100 - perc_cutoff)`` if ``perc_cutoff`` is supplied,
        or ``(1, 99)`` otherwise.
    hard_min, hard_max:
        Optional hard clipping bounds applied after percentile trimming.
    eliminate_zeros:
        Whether zeros should be treated as nodata values during pre-processing.

    Returns
    -------
    tuple
        ``(array, bounds, nodata, crs, options)`` where ``array`` is the
        processed ``float32`` raster ready for plotting, ``bounds`` is a
        :class:`rasterio.coords.BoundingBox`, ``nodata`` holds the nodata
        marker (if any), ``crs`` is the raster CRS, and ``options`` contains
        the resolved percentile configuration used by the helper.
    """

    resolved_p_min = p_min
    resolved_p_max = p_max

    if resolved_p_min is None or resolved_p_max is None:
        default_cutoff = perc_cutoff if perc_cutoff is not None else 1.0
        if resolved_p_min is None:
            resolved_p_min = default_cutoff
        if resolved_p_max is None:
            resolved_p_max = 100 - default_cutoff

    if isinstance(source, (str, Path)):
        with rasterio.open(source) as src:
            data = src.read(raster_band).astype("float32", copy=False)
            bounds = src.bounds
            nodata = src.nodata
            crs = src.crs
    elif isinstance(source, xr.DataArray):
        if "band" in source.dims:
            band_index = band if band is not None else 0
            arr = source.isel(band=band_index)
        else:
            arr = source

        data = arr.values.astype("float32", copy=False)
        if hasattr(arr, "rio"):
            nodata = arr.rio.nodata
            minx, miny, maxx, maxy = arr.rio.bounds()
            bounds = BoundingBox(left=minx, bottom=miny, right=maxx, top=maxy)
            crs = arr.rio.crs
        else:  # pragma: no cover - defensive guard for unexpected input
            nodata = None
            bounds = BoundingBox(left=0.0, bottom=0.0, right=float(data.shape[1]), top=float(data.shape[0]))
            crs = None
    else:  # pragma: no cover - defensive guard
        raise TypeError(
            "source must be a path or xarray.DataArray"
        )

    processed = _preprocess_raster_data_percentiles(
        data,
        nodata_value=nodata,
        p_min=resolved_p_min,
        p_max=resolved_p_max,
        hard_min=hard_min,
        hard_max=hard_max,
        eliminate_zeros=eliminate_zeros,
    )

    options = {
        "p_min": resolved_p_min,
        "p_max": resolved_p_max,
        "hard_min": hard_min,
        "hard_max": hard_max,
    }

    return processed, bounds, nodata, crs, options


def _create_plt_choropleth(
    raster_data: np.ndarray,
    bounds,
    title: str,
    region: str = '',
    label_title: str = 'Raster Values',
    quantiles: Optional[Union[int, Sequence[float]]] = None,
    cmap: str = "viridis",
    n_categories: int = 20,
    x_size = 14,
    y_size = 8,
    base_shp=None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    divergence_center: Optional[float] = None,
    raster_crs=None,
    plt_show: bool = True,
    truncate_one_sided: bool = False
):
    """
    Plot the raster in its native CRS.
    - NaNs are fully transparent.
    - base_shp (world boundaries) is reprojected to the raster CRS before overlay.
    - Axes are labeled as projected coordinates, not "Longitude/Latitude".
    - When ``truncate_one_sided`` is ``True``, one-sided data (all positive or all negative) uses the corresponding half of the colormap. When ``False`` (default), the full colormap is used while keeping the same normalization.
    - ``quantiles`` can be an int (number of bins) or an explicit sequence of bin edges. If all provided edges fall between 0 and 100, they are interpreted as percentiles of the data.
    """

    # 1) Optionally filter base_shp by region (still in its own CRS at this point)
    shp_to_plot = None
    if base_shp is not None:
        shp_to_plot = base_shp
        if region:
            mask = (
                shp_to_plot['NAME'].str.contains(region, case=False, na=False)
                | shp_to_plot['CONTINENT'].str.contains(region, case=False, na=False)
            )
            shp_to_plot = shp_to_plot[mask]
            if shp_to_plot.empty:
                raise ValueError(f"Region '{region}' not found in world shapefile.")

    # 2) Determine categorical vs continuous
    valid_mask = np.isfinite(raster_data)
    values = raster_data[valid_mask]

    if values.size == 0:
        raise ValueError("Raster is empty or fully NaN after preprocessing.")

    unique_vals = np.unique(values)
    is_categorical = (len(unique_vals) <= n_categories)

    print(
        f"Raster has {unique_vals.size:,} unique values. "
        f"Min: {values.min():,.2f}. Max: {values.max():,.2f}"
    )

    # 3) Build color normalization
    divergence_midpoint = 0.0 if divergence_center is None else float(divergence_center)
    quantile_edges = None
    if quantiles is not None and not isinstance(quantiles, (int, np.integer)):
        edges = np.asarray(quantiles, dtype=float)
        if edges.ndim != 1 or edges.size < 2:
            raise ValueError("quantiles must be a 1D sequence with at least two entries.")
        if np.all((edges >= 0) & (edges <= 100)):
            edges = np.percentile(values, edges)
        quantile_edges = edges

    if is_categorical:
        ncat = len(unique_vals)

        # --- build boundaries (len = ncat+1) ---
        if ncat > 1:
            diffs = np.diff(unique_vals)
            boundaries = np.concatenate([
                [unique_vals[0] - diffs[0] / 2],
                (unique_vals[:-1] + unique_vals[1:]) / 2,
                [unique_vals[-1] + diffs[-1] / 2],
            ])
        else:
            boundaries = np.array([unique_vals[0] - 0.5, unique_vals[0] + 0.5])

        # --- choose colormap: USE the one the user passed ---
        # assume `cmap` is your input argument (string or Colormap)
        if cmap is None:
            cmap_obj = ListedColormap(plt.cm.tab20.colors[:ncat])
        else:
            cmap_obj = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap

            # If it's a continuous cmap, sample it into ncat discrete colors
            if not isinstance(cmap_obj, ListedColormap):
                cmap_obj = ListedColormap(cmap_obj(np.linspace(0, 1, ncat)))

            # If it's categorical but too short, resample/cycle to ncat
            if getattr(cmap_obj, "N", ncat) < ncat:
                cmap_obj = ListedColormap(cmap_obj(np.linspace(0, 1, ncat)))

        # IMPORTANT: ncolors should align with cmap_obj
        norm = BoundaryNorm(boundaries, ncolors=cmap_obj.N, clip=True)

    else:
        # Continuous
        provided_bounds = (vmin is not None) or (vmax is not None)
        if (vmin is None) and (vmax is None):
            vmin = float(np.nanmin(values))
            vmax = float(np.nanmax(values))

        neg = values[values < 0]
        pos = values[values > 0]

        if quantiles is not None:
            if quantile_edges is not None:
                print("Using quantile edges")
                edges = quantile_edges
                n_int = len(edges) - 1
                colors = plt.get_cmap(cmap)(np.linspace(0, 1, n_int))
                cmap_obj = ListedColormap(colors)
                norm = BoundaryNorm(edges, ncolors=n_int)
            else:
                print("Using quantiles")
                if (len(neg) > 0) and (len(pos) > 0):
                    print("2-sided route (quantiles)")
                    edges = np.linspace(vmin, vmax, quantiles + 1)
                    n_int = len(edges) - 1
                    colors = plt.get_cmap(cmap)(np.linspace(0, 1, n_int))
                    cmap_obj = ListedColormap(colors)
                    norm = BoundaryNorm(edges, ncolors=n_int)

                elif len(pos) == 0:
                    print("All negatives route (quantiles)")
                    if truncate_one_sided:
                        cmap_obj = _truncate_colormap(cmap, 0.0, 0.5, n=quantiles)
                    else:
                        cmap_obj = _truncate_colormap(cmap, 0.0, 1.0, n=quantiles)
                    norm_vmin = vmin if (vmin is not None) else 0.0
                    norm_vmax = vmax if provided_bounds else 0.0
                    if provided_bounds and vmax is None:
                        norm_vmax = 0.0
                    norm = Normalize(vmin=norm_vmin, vmax=norm_vmax)

                else:
                    print("All positives route (quantiles)")
                    if truncate_one_sided:
                        cmap_obj = _truncate_colormap(cmap, 0.5, 1.0, n=quantiles)
                    else:
                        cmap_obj = _truncate_colormap(cmap, 0.0, 1.0, n=quantiles)
                    norm_vmin = vmin if provided_bounds else 0.0
                    if provided_bounds and vmin is None:
                        norm_vmin = 0.0
                    norm_vmax = vmax if (vmax is not None) else 0.0
                    norm = Normalize(vmin=norm_vmin, vmax=norm_vmax)
        else:
            print("Not using quantiles")
            if (len(pos) > 0) and (len(neg) > 0):
                print("2-sided route (continuous)")
                cmap_obj = _truncate_colormap(cmap, 0.0, 1.0)
                norm = TwoSlopeNorm(vmin=vmin, vcenter=divergence_midpoint, vmax=vmax)
            elif len(pos) == 0:
                print("All negatives route (continuous)")
                if truncate_one_sided:
                    cmap_obj = _truncate_colormap(cmap, 0.0, 0.5)
                else:
                    cmap_obj = plt.get_cmap(cmap)
                norm_vmin = vmin if (vmin is not None) else 0.0
                norm_vmax = vmax if provided_bounds else 0.0
                if provided_bounds and vmax is None:
                    norm_vmax = 0.0
                norm = Normalize(vmin=norm_vmin, vmax=norm_vmax)
            else:
                print("All positives route (continuous)")
                if truncate_one_sided:
                    cmap_obj = _truncate_colormap(cmap, 0.5, 1.0)
                else:
                    cmap_obj = plt.get_cmap(cmap)
                norm_vmin = vmin if provided_bounds else 0.0
                if provided_bounds and vmin is None:
                    norm_vmin = 0.0
                norm_vmax = vmax if (vmax is not None) else 0.0
                norm = Normalize(vmin=norm_vmin, vmax=norm_vmax)

    # 4) Start figure/axes
    fig, ax = plt.subplots(figsize=(x_size,y_size))
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

    img = ax.imshow(
        raster_data,
        extent=extent,
        origin="upper",
        cmap=cmap_obj,
        norm=norm,
        # no alpha here; NaN pixels are automatically transparent
    )

    # 5) Plot colorbar
    cbar = fig.colorbar(img, ax=ax, label=label_title)
    if is_categorical:
        cbar.set_ticks(unique_vals)
        cbar.set_ticklabels(unique_vals)

    # 6) Reproject and overlay shapefile in raster CRS
    if shp_to_plot is not None and raster_crs is not None:
        try:
            # Reproject shapefile to raster CRS
            shp_proj = shp_to_plot.to_crs(raster_crs)
            shp_proj.boundary.plot(ax=ax, edgecolor='grey', linewidth=0.5)
        except Exception as e:
            print("Warning: could not reproject/plot base_shp:", e)

    # 7) Title and axis labels
    ax.set_title(title)

    if raster_crs is not None:
        # We know this isn't lon/lat, it's the raster's projected space (meters-ish)
        ax.set_xlabel("Projected X (map units)")
        ax.set_ylabel("Projected Y (map units)")
    else:
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

    if plt_show:
        plt.show()
    else:
        return fig, ax


def _truncate_colormap(cmap_name, minval=0.0, maxval=1.0, n=256):
    """
    Return a new colormap sampled from [minval, maxval] of the original.
    """
    base = plt.get_cmap(cmap_name)
    colors = base(np.linspace(minval, maxval, n))
    return LinearSegmentedColormap.from_list(
        f"{cmap_name}_trunc_{minval:.2f}_{maxval:.2f}",
        colors,
        N=n
    )

def plot_raster_on_world_extremes_cutoff(
    raster: Union[str, Path, xr.DataArray],
    title: str,
    label_title: str = 'Raster Values',
    raster_band: int = 1,
    band: Optional[int] = None,
    perc_cutoff: Optional[float] = 1,
    p_min: Optional[float] = None,
    p_max: Optional[float] = None,
    quantiles: Optional[Union[int, Sequence[float]]] = None,
    region: Optional[str] = None,
    cmap: str = 'viridis',
    divergence_center: Optional[float] = None,
    n_categories: int = 20,
    base_shp: Optional[gpd.GeoDataFrame] = None,
    plt_show: bool = True,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
    eliminate_zeros: bool = False,
    diverg0: Optional[bool] = None,
    truncate_one_sided: bool = False,
):
    """
    Load raster, clip extremes, and plot in its native projection.
    - Correctly handles nodata (e.g. -32000) by masking to NaN.
    - Reprojects base_shp (assumed EPSG:4326 or anything else) to raster CRS automatically.
    - When ``truncate_one_sided`` is ``True``, one-sided data uses the
      corresponding half of the colormap. When ``False`` (default), the full
      colormap is used while honoring divergence-centered normalization bounds.
    - Explicit ``min_val``/``max_val`` bounds take precedence over divergence
      centering when provided.
    - ``quantiles`` can be an int (number of bins) or an explicit sequence of
      bin edges. If all provided edges fall between 0 and 100, they are
      interpreted as percentiles of the data.
    """

    if base_shp is None:
        base_shp = _get_world_map()

    raster_data, bounds, _, raster_crs, _ = _prepare_raster_plot_input(
        raster,
        raster_band=raster_band,
        band=band,
        perc_cutoff=perc_cutoff,
        p_min=p_min,
        p_max=p_max,
        hard_min=min_val,
        hard_max=max_val,
        eliminate_zeros=eliminate_zeros,
    )

    resolved_divergence = divergence_center
    if diverg0 is not None:
        if divergence_center is not None and diverg0:
            raise ValueError("Specify either 'divergence_center' or 'diverg0', not both.")
        if diverg0:
            resolved_divergence = 0.0
        elif divergence_center is None:
            resolved_divergence = None

    vmin = min_val
    vmax = max_val
    if (vmin is None) and (vmax is None) and (resolved_divergence is not None):
        valid = raster_data[np.isfinite(raster_data)]
        if valid.size:
            max_dev = np.max(np.abs(valid - resolved_divergence))
            vmin = resolved_divergence - max_dev
            vmax = resolved_divergence + max_dev
        else:
            vmin = vmax = resolved_divergence

    fig_ax = _create_plt_choropleth(
        raster_data=raster_data,
        bounds=bounds,
        title=title,
        region=region,
        label_title=label_title,
        quantiles=quantiles,
        cmap=cmap,
        n_categories=n_categories,
        base_shp=base_shp,
        vmin=vmin,
        vmax=vmax,
        divergence_center=resolved_divergence,
        truncate_one_sided=truncate_one_sided,
        raster_crs=raster_crs,
        plt_show=plt_show
    )

    if not plt_show:
        return fig_ax


def plot_da_on_world_extremes_cutoff(
    da: xr.DataArray,
    title: str,
    label_title: str = 'Raster Values',
    band: Optional[int] = None,
    perc_cutoff: float = 1.0,
    p_min: Optional[float] = None,
    p_max: Optional[float] = None,
    quantiles: Optional[int] = None,
    region: Optional[str] = None,
    cmap: str = 'viridis',
    diverg0: bool = False,
    n_categories: int = 20,
    base_shp: Optional[gpd.GeoDataFrame] = None,
    eliminate_zeros: bool = False
):
    """Backward compatible wrapper delegating to
    :func:`plot_raster_on_world_extremes_cutoff`.
    """

    result = plot_raster_on_world_extremes_cutoff(
        da,
        title,
        label_title=label_title,
        band=band,
        perc_cutoff=perc_cutoff,
        p_min=p_min,
        p_max=p_max,
        quantiles=quantiles,
        region=region,
        cmap=cmap,
        divergence_center=0.0 if diverg0 else None,
        n_categories=n_categories,
        base_shp=base_shp,
        plt_show=False,
        eliminate_zeros=eliminate_zeros,
    )

    return result


def plot_all_raster_bands(
    tif_path: str,
    title_prefix: str = "",
    max_cols: int = 3,
    cmap: str = "viridis",
    quantiles: Optional[int] = None,
    n_categories: int = 20,
    base_shp: Optional[gpd.GeoDataFrame] = None,
    x_size: float = 14,
    y_size: float = 8
):
    """
    Plot every band of a multi-band raster in a grid (up to max_cols per row).
    Each subplot gets its own colorbar and uses categorical logic if there are <= n_categories unique values.
    """
    if base_shp is None:
        base_shp = _get_world_map()

    # 1) Load the raster and metadata
    with rasterio.open(tif_path) as src:
        bands  = src.read()         # shape: (count, height, width)
        bounds = src.bounds
        nodata = src.nodata

    count, _, _ = bands.shape
    cols  = min(count, max_cols)
    rows  = int(np.ceil(count / cols))

    # 2) Prepare figure + axes
    fig, axes = plt.subplots(rows, cols,
                             figsize=(cols * x_size, rows * y_size),
                             squeeze=False)
    axes_flat = axes.flatten()

    # 3) Loop over bands
    for idx in range(count):
        ax   = axes_flat[idx]
        data = bands[idx].astype("float32")

        # mask nodata
        data[data == nodata] = np.nan

        # decide categorical vs continuous
        unique = np.unique(data[np.isfinite(data)])
        is_cat = len(unique) <= n_categories

        if is_cat:
            # build discrete norm
            if len(unique) > 1:
                first = unique[0] - (unique[1] - unique[0]) / 2
                mids  = [(unique[i-1] + unique[i]) / 2 for i in range(1, len(unique))]
                last  = unique[-1] + (unique[-1] - unique[-2]) / 2
                bounds_list = [first] + mids + [last]
            else:
                bounds_list = [unique[0] - 0.05, unique[0] + 0.05]

            discrete_cmap = ListedColormap(plt.cm.tab20.colors[:len(unique)])
            norm          = BoundaryNorm(bounds_list, len(unique))
            plot_cmap     = discrete_cmap
        else:
            plot_cmap = plt.get_cmap(cmap)
            if quantiles:
                vals = data[np.isfinite(data)]
                bins = np.quantile(vals, np.linspace(0, 1, quantiles + 1))
                norm = BoundaryNorm(bins, plot_cmap.N, extend="both")
            else:
                norm = None

        # plot the raster
        extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
        im = ax.imshow(data,
                       cmap=plot_cmap,
                       norm=norm,
                       extent=extent,
                       alpha=0.7)

        # overlay boundaries
        base_shp.boundary.plot(ax=ax, color="grey", linewidth=0.5)

        # titles & axes
        ax.set_title(f"{title_prefix} – Band {idx+1}", fontsize=10)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

        # colorbar
        cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        if is_cat:
            cbar.set_ticks(unique)
            cbar.set_ticklabels(unique)

    # 4) turn off any extra subplots
    for j in range(count, len(axes_flat)):
        axes_flat[j].axis("off")

    plt.tight_layout()
    plt.show()


def plot_raster_on_world_no_min(
    tif_path,
    title: str,
    label_title='Raster Values',
    no_data_value=None,
    threshold=None,
    quantiles=None,
    cmap='viridis',
    x_size=14,
    y_size=8,
    base_shp: Optional[gpd.GeoDataFrame] = None,
):
    if base_shp is None:
        base_shp = _get_world_map()

    # Load the raster data
    with rasterio.open(tif_path) as src:
        raster_data = src.read(1)  # Read the first band
        bounds = src.bounds  # Get bounds
        nodata_value = src.nodata  # No-data value

    if no_data_value is None:
        no_data_value = nodata_value

    # Preprocess the raster data
    if threshold is None:
        raster_data = _preprocess_raster_data_eliminate_nodata(raster_data, nodata_value=no_data_value)
    else:
        raster_data = _preprocess_raster_data_eliminate_nodata(raster_data, nodata_value=no_data_value, threshold=threshold)

    # Create the plot
    _create_plt_choropleth(raster_data=raster_data, bounds=bounds, title=title, label_title=label_title, quantiles=quantiles, cmap=cmap,
                           x_size=x_size, y_size=y_size, base_shp=base_shp)

    # Show the plot
    plt.show()


def plot_static_shapefile_on_world(shapefile, color_variable_name: str,title="Shapefile Overlay on World Map", alpha=0.6):
    """
    Plots a shapefile on top of a world map.
    
    Parameters:
        shapefile (shp): Shapefile to be plotted
        color_variable_name (str): Name of the variable to be plotted with different colours
        title (str): Title of the plot.
        alpha (float): Opacity of the shapefile geometries (0.0 to 1.0).
        
    Returns:
        None: Displays the plot.
    """

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    _get_world_map().plot(ax=ax, color="lightgrey", edgecolor="black")  # Plot the world map

    gdf = shapefile

    # Check if the variable is categorical or continuous
    if gdf[color_variable_name].dtype.name == 'category' or gdf[color_variable_name].dtype == 'object':
        # Categorical data
        categories = gdf[color_variable_name].unique()
        color_map = {category: plt.cm.tab20(i / len(categories)) for i, category in enumerate(categories)}
        gdf['color'] = gdf[color_variable_name].map(color_map)

        # Plot with categorical colors
        gdf.plot(ax=ax, color=gdf['color'], legend=True, alpha = alpha)
        handles = [plt.Line2D([0], [0], marker='o', color='w', label=cat, 
                            markerfacecolor=color_map[cat], markersize=10) 
                for cat in categories]
        ax.legend(handles=handles, title=color_variable_name)
        
    else:
        # Continuous data
        norm = mcolors.Normalize(vmin=gdf[color_variable_name].min(), vmax=gdf[color_variable_name].max())
        cmap = cm.viridis  # Choose a colormap
        gdf['color'] = gdf[color_variable_name].apply(lambda x: cmap(norm(x)))

        # Plot with continuous colors
        gdf.plot(ax=ax, color=gdf['color'], alpha = alpha)

        # Add a colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label(color_variable_name)

    # Set the title and labels
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("Longitude", fontsize=12)
    ax.set_ylabel("Latitude", fontsize=12)

    # Adjust aspect ratio and grid
    ax.set_aspect('equal')
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.show()


def plot_raster_over_gdf(raster_path: str,
                         gdf: gpd.GeoDataFrame,
                         band: int = 1,
                         title: str = "Raster over GeoDataFrame",
                         cmap: str = "viridis",
                         alpha: float = 0.7,
                         figsize: tuple = (12, 8),
                         polygon_linewidth = 0.5,
                         polygon_edgecolor = "lightgrey"):
    """
    Clip a raster to the given GeoDataFrame and plot only the overlapping pixels,
    honoring the raster's own nodata value by rendering those pixels transparent.

    Parameters:
    -----------
    raster_path : str
        Path to the raster file (e.g., GeoTIFF).
    gdf : GeoDataFrame
        Polygons to clip by and overlay.
    band : int
        1-based index of the raster band to plot.
    title : str
        Plot title.
    cmap : str
        Matplotlib colormap name.
    alpha : float
        Alpha transparency for the raster.
    figsize : tuple
        Figure size.
    """
    # Load and check CRS
    raster = rioxarray.open_rasterio(raster_path)
    if raster.rio.crs is None:
        raise ValueError("Raster must have a CRS.")

    # Reproject the GeoDataFrame
    gdf_proj = gdf.to_crs(raster.rio.crs)

    # Clip raster by geometry
    clipped = raster.rio.clip(gdf_proj.geometry, drop=True)

    # Extract the desired band
    try:
        da = clipped.sel(band=band)
    except ValueError:
        da = clipped.isel(band=band - 1)

    data = da.values
    nodata_val = clipped.rio.nodata

    # Mask nodata pixels → turn them into NaN
    if nodata_val is not None:
        data = np.where(data == nodata_val, np.nan, data)

    # Prepare colormap so that NaNs are transparent
    cmap_obj = cm.get_cmap(cmap)
    cmap_obj.set_bad(color="none")

    # Compute plotting extent
    minx, miny, maxx, maxy = clipped.rio.bounds()
    extent = [minx, maxx, miny, maxy]

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    img = ax.imshow(
        data,
        extent=extent,
        origin="upper",
        cmap=cmap_obj,
        alpha=alpha
    )

    # Overlay polygon boundaries
    gdf_proj.boundary.plot(ax=ax, edgecolor=polygon_edgecolor, linewidth=polygon_linewidth)

    # Colorbar and labels
    cbar = fig.colorbar(img, ax=ax, label="Value")
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.show()

def plot_raster_over_gdf_showpolygonvalues(
    raster_path: str,
    gdf: gpd.GeoDataFrame,
    band: int = 1,
    title: str = "Raster over GeoDataFrame",
    cmap: str = "viridis",
    alpha: float = 0.7,
    figsize: tuple = (12, 8),
    polygon_linewidth: float = 0.5,
    polygon_edgecolor: str = "lightgrey",
    # polygons
    polygon_value_col: str | None = None,
    polygon_alpha: float = 0.6,
    polygon_cmap: str | None = None,
    # scale/linking
    match_color_scale: bool = True,
    add_polygon_colorbar_when_unmatched: bool = True,
    # labels
    raster_colorbar_label: str = "Raster value",
    polygon_colorbar_label: str = "Polygon value",
    # NEW: quantiles
    n_quantiles: int = 0,   # 0 = continuous, >=1 = quantile bins
):
    """
    Plot a raster clipped to gdf and (optionally) fill polygons by an attribute.
    Supports quantile-based color binning via `n_quantiles`:
      - n_quantiles = 0  → continuous color scale (Normalize)
      - n_quantiles >= 1 → quantile bins (BoundaryNorm)
    If `match_color_scale=True`, polygons reuse the raster's colormap AND scale
    (continuous or quantile bins), so colors mean the same across layers.
    """
    # Load raster & CRS
    raster = rioxarray.open_rasterio(raster_path)
    if raster.rio.crs is None:
        raise ValueError("Raster must have a CRS.")
    gdf_proj = gdf.to_crs(raster.rio.crs)

    # Clip and extract band
    clipped = raster.rio.clip(gdf_proj.geometry, drop=True)
    try:
        da = clipped.sel(band=band)
    except ValueError:
        da = clipped.isel(band=band - 1)

    data = da.values
    nodata_val = clipped.rio.nodata
    if nodata_val is not None:
        data = np.where(data == nodata_val, np.nan, data)

    # Base cmap with NaNs transparent
    base_cmap = cm.get_cmap(cmap).copy()
    base_cmap.set_bad(color="none")

    # Raster scale: continuous or quantile
    valid_data = data[~np.isnan(data)]
    if valid_data.size == 0:
        raise ValueError("No valid raster pixels after clipping.")

    if n_quantiles and n_quantiles >= 1:
        # Build quantile bin edges from raster
        probs = np.linspace(0, 1, n_quantiles + 1)
        boundaries = np.quantile(valid_data, probs)
        # Ensure strictly increasing boundaries (may collapse if data are constant)
        boundaries = np.unique(boundaries)
        if boundaries.size < 2:
            # fallback to min/max
            vmin, vmax = float(np.nanmin(valid_data)), float(np.nanmax(valid_data))
            boundaries = np.array([vmin, vmax])
        raster_norm = mcolors.BoundaryNorm(boundaries, ncolors=base_cmap.N, clip=False)
        raster_boundaries = boundaries  # keep for colorbar and polygons
        poly_scale_kind = "quantiles"
    else:
        vmin, vmax = float(np.nanmin(valid_data)), float(np.nanmax(valid_data))
        raster_norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        raster_boundaries = None
        poly_scale_kind = "continuous"

    # Extent
    minx, miny, maxx, maxy = clipped.rio.bounds()
    extent = [minx, maxx, miny, maxy]

    fig, ax = plt.subplots(figsize=figsize)

    # 1) Polygons beneath raster
    drew_polys = False
    if polygon_value_col is not None:
        if polygon_value_col not in gdf_proj.columns:
            raise ValueError(f"Column '{polygon_value_col}' not found in GeoDataFrame.")
        series = gdf_proj[polygon_value_col].astype(float)
        finite_vals = series[np.isfinite(series)]

        if match_color_scale:
            # Reuse raster scale/cmap
            poly_cmap = base_cmap
            if poly_scale_kind == "quantiles":
                poly_norm = mcolors.BoundaryNorm(raster_boundaries, ncolors=poly_cmap.N, clip=False)
            else:
                poly_norm = raster_norm
        else:
            # Independent polygon scale
            poly_cmap = cm.get_cmap(polygon_cmap or cmap)
            if finite_vals.size == 0:
                poly_norm = None
            else:
                if n_quantiles and n_quantiles >= 1:
                    probs = np.linspace(0, 1, n_quantiles + 1)
                    poly_bounds = np.quantile(finite_vals, probs)
                    poly_bounds = np.unique(poly_bounds)
                    if poly_bounds.size < 2:
                        poly_bounds = np.array([finite_vals.min(), finite_vals.max()])
                    poly_norm = mcolors.BoundaryNorm(poly_bounds, ncolors=poly_cmap.N, clip=False)
                else:
                    poly_norm = mcolors.Normalize(vmin=float(finite_vals.min()),
                                                  vmax=float(finite_vals.max()))

        gdf_proj.plot(
            column=polygon_value_col,
            ax=ax,
            cmap=poly_cmap,
            norm=poly_norm,
            alpha=polygon_alpha,
            linewidth=0,
            legend=False
        )
        drew_polys = True

    # 2) Raster on top
    img = ax.imshow(
        data,
        extent=extent,
        origin="upper",
        cmap=base_cmap,
        norm=raster_norm,
        alpha=alpha
    )

    # 3) Polygon boundaries
    gdf_proj.boundary.plot(ax=ax, edgecolor=polygon_edgecolor, linewidth=polygon_linewidth)

    # 4) Colorbars
    # Raster colorbar
    if n_quantiles and n_quantiles >= 1 and raster_boundaries is not None:
        cbar = fig.colorbar(img, ax=ax, label=raster_colorbar_label, boundaries=raster_boundaries, spacing="proportional")
    else:
        cbar = fig.colorbar(img, ax=ax, label=raster_colorbar_label)

    # Optional separate polygon colorbar (only when not matching)
    if drew_polys and (not match_color_scale) and add_polygon_colorbar_when_unmatched:
        poly_mappable = cm.ScalarMappable(norm=poly_norm, cmap=poly_cmap)  # type: ignore
        poly_mappable.set_array([])
        cax = fig.add_axes([ax.get_position().x0 - 0.035, ax.get_position().y0, 0.02, ax.get_position().height])
        if isinstance(poly_norm, mcolors.BoundaryNorm):
            fig.colorbar(poly_mappable, cax=cax, label=polygon_colorbar_label, boundaries=poly_norm.boundaries, spacing="proportional")
        else:
            fig.colorbar(poly_mappable, cax=cax, label=polygon_colorbar_label)

    # Labels
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    plt.show()
    return fig, ax

###############################
### SHAPE PLOTLY FORMATTING ###
###############################
# Helper function to compute map center
def _calculate_map_center(shapefile: gpd.GeoDataFrame):
    """
    Calculate the geometric center of a group of polygons in a GeoDataFrame.

    Parameters:
    - shapefile (GeoDataFrame): GeoDataFrame containing polygon geometries.

    Returns:
    - center (dict): A dictionary with 'lat' and 'lon' representing the center coordinates.
    """
    if shapefile is None or shapefile.empty:
        raise ValueError("Shapefile is empty or None. Cannot calculate map center.")

    # Ensure CRS compatibility
    shapefile = shapefile.to_crs("EPSG:4326")

    # Calculate the centroid of all geometries combined
    combined_geom = shapefile.union_all()  # Union of all geometries
    centroid = combined_geom.centroid

    # Extract lat/lon from the centroid
    center = {"lat": centroid.y, "lon": centroid.x}
    return center

# Preprocess the GeoDataFrame: reproject, optionally apply log-scale transformation or quantile binning.
def _preprocess_gdf(gdf, category_column, log_scale=False, n_quantile=None):
    gdf = gdf.copy().to_crs("EPSG:4326")
    if log_scale:
        # Calculate an offset to avoid log issues with zeros/negatives
        offset = abs(np.min(gdf[category_column])) / 100
        gdf[category_column] = np.log10(gdf[category_column] + offset)
    elif n_quantile is not None:
        # Compute quantile bins and capture the bin edges.
        quantile_bins, bins = pd.qcut(
            gdf[category_column], q=n_quantile, retbins=True, labels=False, duplicates='drop'
        )
        # Create a new column that holds the bin midpoint for each observation.
        gdf['quantile_mid'] = quantile_bins.map(lambda x: (bins[x] + bins[x+1]) / 2)
        # Save the bin edges as an attribute (used later to update the colorbar legend)
        gdf.attrs['quantile_bins'] = bins
    return gdf

# Filter the GeoDataFrame to a specific country (using a global boundaries dataset)
def _filter_by_country(gdf, world, country_zoom):
    if country_zoom:
        country_gdf = world[world["NAME_EN"] == country_zoom]
        if country_gdf.empty:
            raise ValueError(f"Country '{country_zoom}' not found in the boundaries dataset.")
        filtered_gdf = gpd.sjoin(gdf, country_gdf, how="inner", predicate="intersects")
        if filtered_gdf.empty:
            raise ValueError(f"No geometries found in the shapefile within country '{country_zoom}'.")
        return filtered_gdf
    return gdf

# Format hover data to display the continuous variable in a friendly format.
def _format_hover_data(gdf, category_column):
    hover_data = {}
    formatted_name = f'formatted_{category_column}'
    if np.nanmedian(gdf[category_column]) < 1e-2:
        gdf[formatted_name] = gdf[category_column].apply(lambda x: f"{x:.2e}")
    else:
        gdf[formatted_name] = gdf[category_column].apply(lambda x: f"{x:.2f}")
    hover_data[formatted_name] = True
    if 'ECO_NAME' in gdf.columns:
        hover_data['ECO_NAME'] = True
    return hover_data

# Determine the color range for continuous mapping
def _determine_range_color(gdf, category_column, log_scale, min_color_value, max_color_value):
    if log_scale:
        offset = abs(np.min(gdf[category_column])) / 100
        if min_color_value is not None and max_color_value is not None:
            return [np.log10(min_color_value + offset), np.log10(max_color_value + offset)]
        else:
            return None
    else:
        if min_color_value is not None and max_color_value is not None:
            return [min_color_value, max_color_value]
        else:
            return None


####################
### RASTER TOOLS ###
####################
def downsample_raster(tif_path, scale_factor, output_path):
    with rasterio.open(tif_path) as src:
        # Calculate new dimensions
        new_width = int(src.width * scale_factor)
        new_height = int(src.height * scale_factor)

        # Resample raster data
        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=Resampling.average  # Use average to aggregate pixel values
        )
        
        # Update transform to match the new dimensions
        transform = src.transform * src.transform.scale(
            src.width / new_width,
            src.height / new_height
        )

        # Save the downsampled raster
        profile = src.profile
        profile.update({
            'height': new_height,
            'width': new_width,
            'transform': transform
        })
        
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(data)

def inspect_raster(file_path):
    """
    Inspect and print metadata and band information for a raster file.
    
    Args:
        file_path (str): Path to the raster file.
    """
    try:
        with rasterio.open(file_path) as src:
                # General metadata
                print(f"File: {file_path}")
                print(f"Driver: {src.driver}")
                print(f"Width, Height: {src.width}, {src.height}")
                print(f"Number of Bands: {src.count}")
                print(f"CRS: {src.crs}")
                print(f"Bounds: {src.bounds}")
                print(f"Pixel Size: {src.res}")
                print(f"No-data Value: {src.nodata}")

                print("\n--- Raster Metadata ---")
                tags = src.tags()
                if tags:
                    for key, value in tags.items():
                        print(f"{key}: {value}")
                else:
                    print("No additional metadata found.")

                print("\n--- Band Information ---")
                nodata = src.nodata

                for b in range(1, src.count + 1):
                    # Read the raw data (no auto-masking)
                    data = src.read(b)

                    # Build a mask of invalid data
                    mask_nodata = (data == nodata) if nodata is not None else np.zeros_like(data, dtype=bool)
                    mask_nan    = np.isnan(data)
                    mask        = mask_nodata | mask_nan

                    # Extract only the valid pixels
                    valid = data[~mask]

                    # Compute stats safely (handle empty valid array)
                    if valid.size:
                        vmin = valid.min()
                        vmax = valid.max()
                        vmean = valid.mean()
                        vstd = valid.std()
                    else:
                        vmin = vmax = vmean = vstd = float('nan')

                    print(f"\nBand {b}:")
                    print(f"  Data Type: {src.dtypes[b-1]}")
                    print(f"  Min Value: {vmin}")
                    print(f"  Max Value: {vmax}")
                    print(f"  Mean Value: {vmean}")
                    print(f"  Standard Deviation: {vstd}")
                
                    # Check for band-specific metadata
                    band_tags = src.tags(b)
                    if band_tags:
                        print("  Band Metadata:")
                        for key, value in band_tags.items():
                            print(f"    {key}: {value}")
        
    except Exception as e:
        print(f"Error reading raster file: {e}")

def get_raster_band_count(file_path):
    """
    Get the number of bands in a raster file.

    Parameters:
    - file_path (str): Path to the raster file.

    Returns:
    - int: The number of bands in the raster file.
    """
    with rasterio.open(file_path) as src:
        return src.count  

def extract_raster_values_from_tif(file_path, mask_invalid=True):
    """
    Extracts all raster values from a .tif file.

    Parameters:
    - file_path (str): Path to the raster .tif file.
    - mask_invalid (bool): Whether to exclude invalid values like NaN or masked values.

    Returns:
    - numpy.ndarray: A 1D array of valid raster values.
    """
    with rasterio.open(file_path) as src:
        raster_data = src.read(1)  # Read the first band
        if mask_invalid:
            return raster_data[np.isfinite(raster_data)].flatten()
        return raster_data.flatten()

def plot_raster_data_histogram(values, bins=50, title="Raster Value Histogram", xlabel="Value", ylabel="Frequency"):
    """
    Plots a histogram of raster values.

    Parameters:
    - values (numpy.ndarray): A 1D array of raster values.
    - bins (int): Number of bins for the histogram.
    - title (str): Title of the histogram.
    - xlabel (str): Label for the x-axis.
    - ylabel (str): Label for the y-axis.
    """
    plt.hist(values, bins=bins, color='blue', alpha=0.7, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()

def plot_raster_histogram(
    raster_path,
    band: int = 1,
    bins: int = 50,
    title: str = "Raster Value Histogram",
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    log_scale: bool = False,
    min_positive: float = 1e-6,
    mask_zeros: bool = False
):
    """
    Plots a histogram of raster values with optional logarithmic scale.

    Parameters
    ----------
    raster_path : str
        Path to the raster file.
    band : int
        Raster band number to read.
    bins : int
        Number of bins for the histogram.
    title : str
        Title of the histogram.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    log_scale : bool
        If True, plot histogram in log10 scale (x-axis).
    min_positive : float
        Minimum positive threshold to avoid log(0) if log_scale=True.
    """

    # --- 1. Read raster
    with rasterio.open(raster_path) as src:
        data = src.read(band).astype("float32")
        nodata = src.nodata

    # --- 2. Replace nodata with NaN and flatten
    mask_nodata = (data == nodata) if nodata is not None else np.zeros_like(data, dtype=bool)
    mask_nan    = np.isnan(data)
    zero_mask  = data == 0
    mask        = (mask_nodata | mask_nan | zero_mask) if mask_zeros else (mask_nodata | mask_nan) 

    # Extract only the valid pixels
    valid = data[~mask]
    
    #if nodata is not None:
    #    data[data == nodata] = np.nan

    values = valid[np.isfinite(valid)].flatten()

    # --- 3. Handle log scale
    if log_scale:
        # remove non-positive values
        valid_mask = values > 0
        n_removed = np.count_nonzero(~valid_mask)
        if n_removed > 0:
            print(f"⚠️ Removed {n_removed:,} non-positive values before log-transform.")
        values = values[valid_mask]
        # avoid zeros / negatives
        values = np.log10(np.clip(values, min_positive, None))
        xlabel = f"log₁₀({xlabel})"

    # --- 4. Plot
    plt.figure(figsize=(8, 5))
    plt.hist(values, bins=bins, color='steelblue', alpha=0.75, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

def plot_overlapping_histograms(raster1_path, raster2_path, title, x_label, label1, label2, bins=50, std_dev_filter=0, filter_quantiles=0.0, quantiles_tails = 'both'):
    """
    Plot overlapping histograms of two rasters, filtering values beyond a specified number of standard deviations from the mean.
    
    Parameters:
    -----------
    raster1_path : str
        Path to first raster file
    raster2_path : str  
        Path to second raster file
    title : str
        Plot title
    x_label : str
        X-axis label
    label1 : str
        Legend label for first raster
    label2 : str
        Legend label for second raster
    bins : int
        Number of histogram bins
    std_dev_filter : float
        Number of standard deviations from mean to filter data
    """
    # Check if both standard deviation and filter quantiles are applied
    if (std_dev_filter > 0)  and filter_quantiles > 0:
        raise ValueError("Cannot apply both standard deviation and quantile filters at the same time. Choose 1")
    if quantiles_tails not in ['both', 'left', 'right']:
        raise ValueError("Quantile tails filtering must be either both, left, or right")

    with rasterio.open(raster1_path) as src1, rasterio.open(raster2_path) as src2:
        arr1 = src1.read(1, masked=True)
        arr2 = src2.read(1, masked=True)

    # Also mask any NaN/Inf explicitly
    a1 = np.ma.masked_invalid(arr1)
    a2 = np.ma.masked_invalid(arr2)

    # Flatten and remove masked/nan values
    data1 = a1.compressed()
    data2 = a2.compressed()

    if std_dev_filter >0:
        # Filter outliers based on standard deviations
        mean1, std1 = np.mean(data1), np.std(data1)
        mean2, std2 = np.mean(data2), np.std(data2)

        print(f"Raster 1 mean: {mean1} and std_dev: {std1}")
        print(f"Raster 2 mean: {mean2} and std_dev: {std2}")

        mask1 = (data1 >= (mean1 - std_dev_filter*std1)) & (data1 <=( mean1 + std_dev_filter*std1))
        mask2 = (data2 >= (mean2 - std_dev_filter*std2)) & (data2 <=( mean2 + std_dev_filter*std2))

        data1_f = data1[mask1]
        data2_f = data2[mask2]
    elif filter_quantiles > 0:
        # Compute thresholds
        if quantiles_tails == 'both':
            qlow1, qhigh1 = np.quantile(data1, [filter_quantiles, 1 - filter_quantiles])
            qlow2, qhigh2 = np.quantile(data2, [filter_quantiles, 1 - filter_quantiles])
            mask1 = (data1 >= qlow1) & (data1 <= qhigh1)
            mask2 = (data2 >= qlow2) & (data2 <= qhigh2)
        elif quantiles_tails == 'left':   # drop lower tail only
            qlow1 = np.quantile(data1, filter_quantiles)
            qlow2 = np.quantile(data2, filter_quantiles)
            mask1 = data1 >= qlow1
            mask2 = data2 >= qlow2
        else:  # 'right' -> drop upper tail only
            qhigh1 = np.quantile(data1, 1 - filter_quantiles)
            qhigh2 = np.quantile(data2, 1 - filter_quantiles)
            mask1 = data1 <= qhigh1
            mask2 = data2 <= qhigh2

        data1_f = data1[mask1]
        data2_f = data2[mask2]
    else:
        data1_f = data1
        data2_f = data2

    # Defining bins
    all_vals = np.concatenate([data1_f, data2_f])
    if all_vals.min() == all_vals.max():
        bins = bins  # or keep the user-provided bins
    else:
        bins = np.linspace(all_vals.min(), all_vals.max(), bins+1)

    # Plotting
    plt.figure(figsize=(10,6))
    plt.hist(data1_f, bins=bins, alpha=0.5, label=label1, color='blue')
    plt.hist(data2_f, bins=bins, alpha=0.5, label=label2, color='orange')
    plt.xlabel(x_label)
    plt.ylabel("Frequency")
    plt.title(title)
    plt.legend()
    plt.show()

#########################
#### ROTHC FUNCTIONS ####
#########################

def plot_multiband_raster_timesires(raster_path: str, title: str, xlabel: str = "Year", ylabel: str = 'Mean SOC (t C/ha)', show_iq = False):
    """
    Plot global mean SOC and cumulative CO2 emissions time series.
    
    Parameters
    ----------
    soc_da : xarray.DataArray
        Annual SOC with dims ('year','y','x').
    co2_da : xarray.DataArray
        Annual CO₂ emissions with dims ('year','y','x').
    """
    # Compute spatial means (ignoring NaNs)
    da = rioxarray.open_rasterio(raster_path, masked=True)

    # 3) rename the band‐axis to “time”
    da = da.rename({da.dims[0]: "time"})
    
    # 4) assign integer time steps 1…N (or real dates)
    duration = da.sizes["time"]
    da = da.assign_coords(time = ("time", np.arange(1, duration+1)))

    # Global mean
    ts_mean  = da.mean(dim=["x","y"], skipna=True)              # shape = (time,)
    if show_iq:
        q1 = da.quantile(q=0.25, dim=["x","y"], skipna=True)
        q3 = da.quantile(q=0.75, dim=["x","y"], skipna=True)                   

    # Plot SOC
    plt.figure()
    plt.plot(ts_mean["time"].values, ts_mean.values, marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, "major")
    
    if show_iq:
        plt.fill_between(ts_mean["time"].values, q1.values, q3.values, color='blue', alpha=0.3, label='IQR')

    plt.show()

def plot_average_band_values(raster_paths, labels, title: str = "Average Value per Year for Maize Scenarios", xlabel: str = "Year (Band)", ylabel: str = 'Mean Value', show_iq=False, x_size=14, y_size=8):
    """
    Plot the average value per band (year) for multiple rasters, with optional interquartile shading.

    Parameters
    ----------
    raster_paths : list of str
        List of raster file paths.
    labels : list of str
        List of labels for each raster.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    show_iq : bool
        If True, plot interquartile range shading for each raster.
    x_size, y_size : float
        Figure size.
    """
    plt.figure(figsize=(x_size, y_size))

    for path, label in zip(raster_paths, labels):
        with rasterio.open(path) as src:
            arr = src.read()  # shape: (bands, height, width)
            band_means = np.nanmean(arr, axis=(1, 2))
            plt.plot(np.arange(1, arr.shape[0] + 1), band_means, label=label)
            if show_iq:
                q1 = np.nanpercentile(arr, 25, axis=(1, 2))
                q3 = np.nanpercentile(arr, 75, axis=(1, 2))
                plt.fill_between(np.arange(1, arr.shape[0] + 1), q1, q3, alpha=0.2)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(axis="y", visible=True, color="#f0f0f0")
    plt.tight_layout()
    plt.show()

def plot_soc_distribution(soc_da, year):
    """
    Plot histogram of SOC distribution across all pixels for a given year.
    
    Parameters
    ----------
    soc_da : xarray.DataArray
        Annual SOC with dims ('year','y','x').
    year : int
        Year number (1-based index) to plot histogram for.
    """
    # Select the requested year
    soc_year = soc_da.sel(year=year)
    # Flatten and remove NaNs
    values = soc_year.values.flatten()
    values = values[np.isfinite(values)]

    plt.figure()
    plt.hist(values, bins=50)
    plt.xlabel('SOC (t C/ha)')
    plt.ylabel('Pixel Count')
    plt.title(f'SOC Distribution in Year {year}')
    plt.grid(True)
    plt.show()
