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
from typing import Optional, Union

# World shapefile (loaded outside the function)
world_global_lr = gpd.read_file("../data/world_maps/low_res/ne_110m_admin_0_countries.shp")
world_global_hr = gpd.read_file("../data/world_maps/high_res/ne_10m_admin_0_countries.shp")


# FUNCTIONS
def preprocess_raster_data_eliminate_nodata(raster_data, nodata_value=None):
    """
    Preprocess raster data to handle invalid and extreme values.
    
    Parameters:
        raster_data (ndarray): The raster data array.
        nodata_value (float): Value representing no-data in the raster.
    
    Returns:
        ndarray: Cleaned raster data array.
    """
    # Replace no-data values with NaN
    if nodata_value is not None:
        raster_data = np.where(raster_data == nodata_value, np.nan, raster_data)
    
    # Replace inf/-inf values with NaN
    raster_data = np.where(np.isfinite(raster_data), raster_data, np.nan)
    
    return raster_data




def preprocess_raster_data_percentiles(
    raster_data, 
    nodata_value=None, 
    p_min: Union[int, float] = 1, 
    p_max: Union[int, float] = 99
):
    """
    Preprocess raster data to handle invalid and extreme values.
    
    Parameters:
        raster_data (ndarray): The raster data array.
        nodata_value (float): Value representing no-data in the raster.
        p_min (int or float): Lower percentile cutoff.
        p_max (int or float): Upper percentile cutoff.
    
    Returns:
        ndarray: Cleaned raster data array.
    """
    # Eliminate no_data and Inf values
    raster_data = preprocess_raster_data_eliminate_nodata(raster_data, nodata_value)
    
    # Clip raster data to a reasonable range (e.g., 1st and 99th percentiles)
    finite_data = raster_data[np.isfinite(raster_data)]
    lower, upper = np.percentile(finite_data, [p_min, p_max])  # 1st and 99th percentiles
    raster_data = np.clip(raster_data, lower, upper)
    
    return raster_data


def preprocess_raster_data_eliminate_low_values(raster_data, nodata_value=None, threshold = None):
    """
    Preprocess raster data to handle invalid and extreme values.
    
    Parameters:
        raster_data (ndarray): The raster data array.
        nodata_value (float): Value representing no-data in the raster.
    
    Returns:
        ndarray: Cleaned raster data array.
    """
    # Eliminate no_data and Inf values
    raster_data = preprocess_raster_data_eliminate_nodata(raster_data, nodata_value)

    # Determine minimum value for as threshold if none is given
    if threshold is None:
        threshold = np.min(raster_data)
    
    # Elimiante minimum values
    finite_data = raster_data[np.isfinite(raster_data)]
    raster_data = np.where(finite_data <= threshold, np.nan, finite_data)
    
    return raster_data


def create_plt_choropleth(
    raster_data,
    bounds,
    title: str,
    region: str = '',
    label_title: str = 'Raster Values',
    quantiles: Optional[int] = None,
    cmap: str = "viridis",
    n_categories: int = 20,
    figsize=(14, 8),
    base_shp=None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    plt_show = True
):
    """
    Plot a raster as a choropleth, with categorical or continuous (optionally diverging) scales.
    """
    # 1) Region filter
    if region and base_shp is not None:
        mask = (
            base_shp['NAME'].str.contains(region, case=False, na=False)
            | base_shp['CONTINENT'].str.contains(region, case=False, na=False)
        )
        base_shp = base_shp[mask]
        if base_shp.empty:
            raise ValueError(f"Region '{region}' not found in world shapefile.")

    # 2) Detect categorical
    valid_mask = np.isfinite(raster_data)
    values = raster_data[valid_mask]
    unique = np.unique(values)
    is_categorical = len(unique) <= n_categories

    # Print some basic data from the raster
    print(f'Raster has {unique.size:,} different values. Min: {values.min():,.2f}. Max: {values.max():,.2f}')

    # Prepare colormap and norm
    if is_categorical:
        # Discrete bins around each unique value
        if len(unique) > 1:
            diffs = np.diff(unique)
            boundaries = np.concatenate([
                [unique[0] - diffs[0] / 2],
                (unique[:-1] + unique[1:]) / 2,
                [unique[-1] + diffs[-1] / 2],
            ])
        else:
            boundaries = [unique[0] - 0.5, unique[0] + 0.5]

        cmap_obj = ListedColormap(plt.cm.tab20.colors[: len(unique)])
        norm = BoundaryNorm(boundaries, ncolors=len(unique))

    else:
        # Continuous data
        values = raster_data[np.isfinite(raster_data)]
        if (vmin is None) and (vmax is None):
            vmin, vmax = np.min(values), np.max(values)
        neg = values[values <  0]
        pos = values[values >  0]

        if quantiles is not None:
            print("Using quantiles")
            # 1) two‑sided split
            if len(neg) > 0 and len(pos) > 0:
                print("2-sided route")
                edges = np.linspace(vmin, vmax, quantiles + 1)

                n_int = len(edges) - 1
                # sample the *diverging* cmap at n_int points
                colors   = plt.get_cmap(cmap)(np.linspace(0, 1, n_int))
                cmap_obj = ListedColormap(colors)
                norm     = BoundaryNorm(edges, ncolors=n_int)

            # 2) all negative → lower half of diverging cmap
            elif len(pos) == 0:
                print("All negatives route")
                cmap_obj = truncate_colormap(cmap, 0.0, 0.5, n=quantiles)
                norm     = Normalize(vmin=vmin, vmax=0)

            # 3) all positive → upper half of diverging cmap
            else:  # len(neg) == 0
                print("All positives route")
                cmap_obj = truncate_colormap(cmap, 0.5, 1.0, n=quantiles)
                norm     = Normalize(vmin=0, vmax=vmax)
        else:
            print("Not using quantiles")
            # continuous scale, whether or not 0 lies inside
            if (len(pos) > 0) & (len(neg) >0):
                print("2-sided route")
                cmap_obj = truncate_colormap(cmap, 0.0, 1.0)
                norm     = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
            elif len(pos) == 0:
                print("All negatives route")
                cmap_obj = truncate_colormap(cmap, 0.0, 0.5)
                norm     = Normalize(vmin=vmin, vmax=0)
            else:
                print("All positives route")
                cmap_obj = truncate_colormap(cmap, 0.5, 1.0)
                norm     = Normalize(vmin=0, vmax=vmax)

    # 3) Plot
    fig, ax = plt.subplots(figsize=figsize)
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

    img = ax.imshow(
        raster_data,
        extent=extent,
        cmap=cmap_obj,
        norm=norm,
        alpha=0.7,
    )

    cbar = fig.colorbar(img, ax=ax, label=label_title)
    if is_categorical:
        cbar.set_ticks(unique)
        cbar.set_ticklabels(unique)

    if base_shp is not None:
        base_shp.boundary.plot(ax=ax, edgecolor='grey', linewidth=0.5)

    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    if plt_show:
        plt.show()
    else:
        return fig, ax


def truncate_colormap(cmap_name, minval=0.0, maxval=1.0, n=256):
    """
    Return a new colormap sampled from [minval, maxval] of the original.
    """
    base = plt.get_cmap(cmap_name)
    colors = base(np.linspace(minval, maxval, n))
    return LinearSegmentedColormap.from_list(
        f"{cmap_name}_trunc_{minval:.2f}_{maxval:.2f}", colors, N=n
    )

def plot_raster_on_world_extremes_cutoff(tif_path, title: str, label_title='Raster Values', raster_band = 1, 
                                         alpha=1, p_min=None, p_max=None, quantiles=None, 
                                         region: Optional[str] = None, cmap='viridis', n_categories = 20, base_shp = world_global_hr):
    """
    Plot a raster file with optional preprocessing and region filtering.
    """
    # Load the raster data
    with rasterio.open(tif_path) as src:
        raster_data = src.read(raster_band)  # Read the first band
        bounds = src.bounds
        nodata_value = src.nodata
    
    # Calculate percentile cutoffs
    if p_min is None:
        p_min = alpha
    if p_max is None:
        p_max = 100 - alpha

    # Preprocess the raster data
    raster_data = preprocess_raster_data_percentiles(raster_data, nodata_value, p_min, p_max)

    # Create the plot
    create_plt_choropleth(
        raster_data=raster_data, 
        bounds=bounds, 
        title=title, 
        region=region, 
        label_title=label_title, 
        quantiles=quantiles, 
        cmap=cmap,
        base_shp=base_shp
    )


def plot_da_on_world_extremes_cutoff(
    da: xr.DataArray,
    title: str,
    label_title: str = 'Raster Values',
    band: Optional[int] = None,
    alpha: float = 1.0,
    p_min: Optional[float] = None,
    p_max: Optional[float] = None,
    quantiles: Optional[int] = None,
    region: Optional[str] = None,
    cmap: str = 'viridis',
    diverg0 = False,
    n_categories: int = 20,
    base_shp = world_global_hr,
    plt_show = False
):
    """
    Plot a single‑band xarray DataArray on a world basemap, with
    percentile cutoff pre‑processing and optional region filter.
    
    Parameters
    ----------
    da : xarray.DataArray
        A 2D or 3D DataArray with a 'band' dimension (or already 2D).
    title : str
        Plot title.
    band : int, optional
        If `da` is 3D, which band index to plot (0‑based). Otherwise ignored.
    p_min, p_max : float, optional
        Percentile cutoff (1–99). If None, defaults to (α, 100–α).
    alpha : float
        “alpha” used to default p_min/p_max when they are None.
    quantiles : int, optional
        Number of quantile bins (for color breaks).
    region : str, optional
        Country or continent name to zoom into.
    cmap : str
        Matplotlib colormap name.
    n_categories : int
        Max unique values to treat as categorical.
    base_shp : GeoDataFrame
        World map for context.
    """
    # 1) extract the raw numpy array + metadata
    if 'band' in da.dims:
        arr = da.isel(band=band or 0).values.astype(float)
    else:
        arr = da.values.astype(float)
    
    nodata = da.rio.nodata
    raw_bounds = da.rio.bounds()

    minx, miny, maxx, maxy = raw_bounds
    bounds = BoundingBox(left=minx,
                         bottom=miny,
                         right=maxx,
                         top=maxy)
    
    # 2) default percentile cutoffs
    if p_min is None: p_min = alpha
    if p_max is None: p_max = 100 - alpha
    
    # 3) run your existing preprocessors
    raster_data = preprocess_raster_data_percentiles(arr, nodata, p_min, p_max)
    
    # 3.5) Optional - Sets cmap to divergence point 0
    if diverg0:
        valid_data = raster_data[np.isfinite(raster_data)]
        max_abs = np.max(np.abs(valid_data))
        vmin, vmax = -max_abs, max_abs
        cmap = cmap  # or any diverging cmap you prefer
    else:
        vmin = vmax = None  # allow normal behavior
    
    # 4) call your choropleth plotter
    fig, ax = create_plt_choropleth(
        raster_data=raster_data,
        bounds=bounds,
        title=title,
        region=region,
        label_title=label_title,
        quantiles=quantiles,
        cmap=cmap,
        n_categories=n_categories,
        base_shp=base_shp,
        vmin = vmin,
        vmax = vmax,
        plt_show = plt_show
    )

    return fig, ax


def plot_all_raster_bands(
    tif_path: str,
    title_prefix: str = "",
    max_cols: int = 3,
    cmap: str = "viridis",
    quantiles: Optional[int] = None,
    n_categories: int = 20,
    base_shp = world_global_hr,
    x_size: float = 14,
    y_size: float = 8
):
    """
    Plot every band of a multi-band raster in a grid (up to max_cols per row).
    Each subplot gets its own colorbar and uses categorical logic if there are <= n_categories unique values.
    """
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


def plot_raster_on_world_no_min(tif_path, title: str, label_title = 'Raster Values', no_data_value = None, threshold = None, quantiles = None, cmap='viridis', x_size = 14, y_size = 8, base_shp = world_global_hr):
    # Load the raster data
    with rasterio.open(tif_path) as src:
        raster_data = src.read(1)  # Read the first band
        bounds = src.bounds  # Get bounds
        nodata_value = src.nodata  # No-data value

    if no_data_value is None:
        no_data_value = nodata_value

    # Preprocess the raster data
    if threshold is None:
        raster_data = preprocess_raster_data_eliminate_nodata(raster_data, nodata_value=no_data_value)
    else:
        raster_data = preprocess_raster_data_eliminate_nodata(raster_data, nodata_value=no_data_value, threshold=threshold)

    # Create the plot
    create_plt_choropleth(raster_data=raster_data, bounds=bounds, title=title, label_title=label_title, quantiles=quantiles, cmap=cmap,
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
    world_global_hr.plot(ax=ax, color="lightgrey", edgecolor="black")  # Plot the world map

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


def plotly_shapefile_continuous(shapefile, category_column, title=None, subtitle=None, legend_title=None, country_zoom=None, log_scale=False, n_quantile=None, color_palette='Viridis', min_color_value=None, max_color_value=None):
    """
    Plots a shapefile using Plotly, showing only data for a specific country based on spatial location.
    
    Parameters:
        shapefile (GeoDataFrame): The shapefile or GeoJSON to be plotted.
        category_column (str): Column with continuous values for coloring.
        title (str, optional): Title of the plot.
        subtitle (str, optional): Subtitle of the plot.
        legend_title (str, optional): Title of the legend (colorbar).
        country_zoom (str, optional): Name of the country to zoom into.
        log_scale (bool, optional): If True, apply a log10 transformation to the data.
        n_quantile (int, optional): Number of quantiles to bin the data into (for discrete coloring).
        color_palette (str, optional): Color scale to use (default is 'Viridis').
        min_color_value (float, optional): Manual minimum value for the color scale.
        max_color_value (float, optional): Manual maximum value for the color scale.
    
    Returns:
        None: Displays the interactive map.
    """
    # Loading shapefile
    gdf = shapefile

    # Ensure geometries are valid
    gdf = gdf[gdf.geometry.notnull()]

    # Check if category_column is in the shapefile
    if category_column not in gdf.columns:
        print('Variable to plot is not in shapefile.')
        return

    # Load country boundaries (default to Natural Earth dataset if not provided)
    # if country_boundaries is not None:
    #     world = gpd.read_file(country_boundaries)

    # Ensure CRS compatibility
    print('Ensuring same projection is used')
    gdf = gdf.to_crs("EPSG:4326")
    world = world_global_hr.to_crs("EPSG:4326")

    # Ensure only one transformation is applied.
    if log_scale and n_quantile is not None:
        raise ValueError("Please choose either log_scale or n_quantiles, not both.")


    # Apply log-scale transformation if requested
    if log_scale:
        print('Transforming into log10')
        # Calculate a small offset to avoid issues with zeros/negatives
        gdf_min = abs(np.min(gdf[category_column]))
        gdf[category_column] = np.log10(gdf[category_column] + gdf_min / 100)


    # Apply quantile binning if requested (only if not using log_scale)
    if n_quantile is not None and not log_scale:
        print(f'Using {n_quantile} quantiles')
        gdf['quantile_bins'], q_bins = pd.qcut(gdf[category_column], q=n_quantile, labels=False, retbins=True)


    # Filter for a specific country if country_zoom is provided
    if country_zoom:
        print("Filtering geometries by country")
        country_gdf = world[world["NAME_EN"] == country_zoom]
        if country_gdf.empty:
            raise ValueError(f"Country '{country_zoom}' not found in the boundaries dataset.")
        # Spatial join to filter geometries within the specified country
        filtered_gdf = gpd.sjoin(gdf, country_gdf, how="inner", predicate="intersects")
        if filtered_gdf.empty:
            raise ValueError(f"No geometries found in the shapefile within country '{country_zoom}'.")
    else:
        filtered_gdf = gdf

    # Convert the filtered GeoDataFrame to GeoJSON format
    print('Converting to geojson')
    geojson_data = filtered_gdf.__geo_interface__

    # Format values in scientific notation for hover data
    print('Formatting values')
    if category_column in filtered_gdf.columns:
        formated_value_name = 'formated_' + category_column
        if np.nanmedian(filtered_gdf[category_column]) < 1e-2:
            filtered_gdf[formated_value_name] = filtered_gdf[category_column].apply(lambda x: f"{x:.2e}")
        else:
            filtered_gdf[formated_value_name] = filtered_gdf[category_column].apply(lambda x: f"{x:.2f}")
        hover_data = {formated_value_name: True}
    else:
        hover_data = {}

    # Adding ecoregions to the hover data if present
    if 'ECO_NAME' in filtered_gdf.columns:
        hover_data['ECO_NAME'] = True

    # Set default title if not provided
    if not title and country_zoom:
        title = f"{category_column}"
        if country_zoom:
            title = f"{category_column} data for {country_zoom}"

    # Determine the range_color parameter for continuous data (only if n_quantile is not used)
    if n_quantile is None:
        if log_scale:
            if min_color_value is not None and max_color_value is not None:
                # Transform the provided min and max values to log-scale using the same offset
                range_color = [np.log10(min_color_value + gdf_min / 100), 
                               np.log10(max_color_value + gdf_min / 100)]
            else:
                range_color = None
        else:
            range_color = [min_color_value, max_color_value] if (min_color_value is not None and max_color_value is not None) else None
    else:
        range_color = None    

    print('Plotting the choropleth map')
    if n_quantile is None:
        fig = px.choropleth_map(
            filtered_gdf,
            geojson=geojson_data,
            locations=filtered_gdf.index,  # Plotly requires an identifier for geometries
            color=category_column,
            hover_data=hover_data,
            opacity=0.6,
            color_continuous_scale=color_palette,
            range_color=range_color,  # Set the color range if provided
            map_style="carto-positron",
            title=title,
            center={
                "lat": 0,  # You can modify these values to recenter the map if needed
                "lon": 0,
            },
            zoom=1,
        )
    else:
        fig = px.choropleth_map(
            filtered_gdf,
            geojson=geojson_data,
            locations=filtered_gdf.index,
            color='quantile_bins',
            hover_data=hover_data,
            opacity=0.6,
            color_continuous_scale=color_palette,
            map_style="carto-positron",
            title=title,
            center={
                "lat": 0,
                "lon": 0,
            },
            zoom=1,
        )
        # If quantile binning is used, update the colorbar to show quantile break values.
        if q_bins.shape[0] > 0:
            # For each bin, show the midpoint as the tick label.
            tickvals = list(range(n_quantile))
            ticktext = [f"{(q_bins[i] + q_bins[i+1]) / 2:.2f}" for i in range(len(q_bins) - 1)]
            fig.update_coloraxes(colorbar=dict(tickvals=tickvals, ticktext=ticktext))

    # Add a subtitle as an annotation if provided
    if n_quantile:
        subtitle = f"Quantiles: {n_quantile}"

    if subtitle:
        fig.update_layout(
            title=dict(
                text=title,
                subtitle=dict(
                    text=subtitle,
                    font=dict(
                        color="gray",
                        style = "italic",
                        size = 16),
                ),
            )
        )

    # Update layout margins
    fig.update_layout(margin={"r": 0, "t": 50, "l": 0, "b": 0})

    # Update legend (colorbar) title if provided
    if legend_title:
        fig.update_coloraxes(colorbar_title=legend_title)

    # Further number formatting in the colorbar ticks
    if np.median(filtered_gdf[category_column]) < 1e-2 and n_quantile is None:
        fig.update_coloraxes(colorbar_tickformat=".2e")  # Scientific notation
    else:
        fig.update_coloraxes(colorbar_tickformat=".2f")  # Two decimal places

    fig.show()


def plotly_shapefile_categorical(shapefile: gpd.GeoDataFrame, categorical_variable: str, title: str, calculate_center = False, in_notebook = False):
    """
    Create a categorical choropleth map using Plotly.

    Parameters:
    - shapefile (GeoDataFrame): GeoDataFrame containing polygons and categorical values.
    - categorical_variable (str): Name of the column containing the categorical values.
    - title (str): Title of the plot.
    """

    # Load the shapefile
    if not shapefile.empty:
        gdf = shapefile
    else:
        print('No shapefile entered')
        return

    # Ensure geometries are valid
    gdf = gdf[gdf.geometry.notnull()]

    # Check if category_column is in the shapefile
    if categorical_variable not in gdf.columns:
        print('Variable to plot is not in shapefile.')
        return

    # Ensure CRS compatibility
    print('Ensuring same projection is used')
    gdf = gdf.to_crs("EPSG:4326")
    world = world_global_hr.to_crs("EPSG:4326")

    print('Converting to geojson')
    fig = px.choropleth_map(
        gdf,
        geojson=gdf.geometry.__geo_interface__,  # GeoJSON representation
        locations=gdf.index,  # Unique identifier for each geometry
        color=categorical_variable,  # Categorical column in your data
        hover_name=categorical_variable,  # Display category in hover
        center = calculate_map_center(gdf) if calculate_center else {"lat":0, "lon": 0},  # Center maps on the middle of all polygons
        zoom = 0,  # Set the zoom level
        # color_discrete_map='Set1',  # Assign specific colors to categories
        title=title,  # Title of the plot
        opacity=0.6,
        map_style="carto-positron"  # Choose map style
    )

    # Update layout to remove the margin and adjust the map zoom
    fig.update_layout(
        margin={"r":0,"t":50,"l":0,"b":0}, 
        geo=dict(showcoastlines=True, 
                 coastlinecolor="Black")
                 )
    
    if in_notebook:
        return fig
    else:
        fig.show()


def plot_raster_over_gdf(raster_path: str,
                         gdf: gpd.GeoDataFrame,
                         band: int = 1,
                         title: str = "Raster over GeoDataFrame",
                         cmap: str = "viridis",
                         alpha: float = 0.7,
                         figsize: tuple = (12, 8)):
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
    gdf_proj.boundary.plot(ax=ax, edgecolor="black", linewidth=1)

    # Colorbar and labels
    cbar = fig.colorbar(img, ax=ax, label="Value")
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.show()



###############################
### SHAPE PLOTLY FORMATTING ###
###############################
# Helper function to compute map center
def calculate_map_center(shapefile: gpd.GeoDataFrame):
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
def preprocess_gdf(gdf, category_column, log_scale=False, n_quantile=None):
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
def filter_by_country(gdf, world, country_zoom):
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
def format_hover_data(gdf, category_column):
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
def determine_range_color(gdf, category_column, log_scale, min_color_value, max_color_value):
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

def plot_raster_histogram(raster_path, band=1, bins=50, title="Raster Value Histogram", xlabel="Value", ylabel="Frequency"):
    """
    Plots a histogram of raster values.

    Parameters:
    - values (numpy.ndarray): A 1D array of raster values.
    - bins (int): Number of bins for the histogram.
    - title (str): Title of the histogram.
    - xlabel (str): Label for the x-axis.
    - ylabel (str): Label for the y-axis.
    """
    # Opens data
    with rasterio.open(raster_path) as src:
        data = src.read(band)
        nodata_value = src.nodata

    # Filters nodata value
    values_wo_nd = np.where(data != nodata_value, data, np.nan)
    valid_values = values_wo_nd[~np.isnan(values_wo_nd)]
    values = valid_values.flatten()

    plt.hist(values, bins=bins, color='blue', alpha=0.7, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()

def plot_overlapping_histograms(raster1_path, raster2_path, title, x_label, label1, label2, bins=50):
    with rasterio.open(raster1_path) as src1, rasterio.open(raster2_path) as src2:
        arr1 = src1.read(1, masked=True)
        arr2 = src2.read(1, masked=True)

    # Flatten and remove masked/nan values
    data1 = arr1.compressed()
    data2 = arr2.compressed()

    plt.figure(figsize=(10,6))
    plt.hist(data1, bins=bins, alpha=0.5, label=label1, color='blue')
    plt.hist(data2, bins=bins, alpha=0.5, label=label2, color='orange')
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