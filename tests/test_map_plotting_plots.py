import numpy as np
import pytest
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from rasterio.transform import from_origin
import rasterio
import rioxarray  # noqa: F401 - ensures the rio accessor is registered

from sbtn_leaf.map_plotting import (
    _prepare_raster_plot_input,
    plot_raster_on_world_extremes_cutoff,
    plot_da_on_world_extremes_cutoff,
    plot_n_rasters_on_world_extremes_cutoff,
)


def _create_test_raster(tmp_path):
    data = np.array([[1, 2], [3, -9999]], dtype="float32")
    transform = from_origin(0, 2, 1, 1)
    path = tmp_path / "test.tif"

    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999,
    ) as dst:
        dst.write(data, 1)

    return path, data, transform


def _create_dataarray(data, transform):
    y_size, x_size = data.shape
    y = np.arange(y_size, dtype=float)
    x = np.arange(x_size, dtype=float)
    da = xr.DataArray(data, dims=("y", "x"), coords={"y": y, "x": x})
    da = da.rio.write_transform(transform)
    da = da.rio.write_crs("EPSG:4326")
    da = da.rio.write_nodata(-9999)
    return da


def test_prepare_raster_input_from_path(tmp_path):
    raster_path, _, _ = _create_test_raster(tmp_path)

    array, bounds, nodata, crs, options = _prepare_raster_plot_input(
        raster_path,
        perc_cutoff=0,
    )

    assert np.isnan(array[-1, -1])
    assert bounds.left == pytest.approx(0)
    assert bounds.top == pytest.approx(2)
    assert nodata == -9999
    assert str(crs) == "EPSG:4326"
    assert options["p_min"] == 0
    assert options["p_max"] == 100


def test_prepare_raster_input_from_dataarray(tmp_path):
    raster_path, data, transform = _create_test_raster(tmp_path)
    da = _create_dataarray(data, transform)

    array, bounds, nodata, crs, options = _prepare_raster_plot_input(
        da,
        perc_cutoff=0,
    )

    assert np.isnan(array[-1, -1])
    expected_left, _, _, expected_top = da.rio.bounds()
    assert bounds.left == pytest.approx(expected_left)
    assert bounds.top == pytest.approx(expected_top)
    assert nodata == -9999
    assert str(crs) == "EPSG:4326"
    assert options["p_min"] == 0
    assert options["p_max"] == 100


def test_plot_raster_accepts_path(tmp_path):
    raster_path, _, _ = _create_test_raster(tmp_path)

    fig, ax = plot_raster_on_world_extremes_cutoff(
        raster_path,
        "Test",
        perc_cutoff=0,
        plt_show=False,
        base_shp=gpd.GeoDataFrame(),
    )

    assert fig is not None and ax is not None
    plt.close(fig)


def test_plot_raster_accepts_dataarray(tmp_path):
    raster_path, data, transform = _create_test_raster(tmp_path)
    da = _create_dataarray(data, transform)

    fig, ax = plot_raster_on_world_extremes_cutoff(
        da,
        "Test",
        perc_cutoff=0,
        plt_show=False,
        base_shp=gpd.GeoDataFrame(),
        divergence_center=0,
    )

    assert fig is not None and ax is not None
    plt.close(fig)


def test_plot_da_wrapper_respects_diverg0(tmp_path):
    raster_path, data, transform = _create_test_raster(tmp_path)
    da = _create_dataarray(data, transform)

    fig, ax = plot_da_on_world_extremes_cutoff(
        da,
        "Test",
        perc_cutoff=0,
        quantiles=5,
        diverg0=True,
        base_shp=gpd.GeoDataFrame(),
    )

    assert fig is not None and ax is not None
    plt.close(fig)


def test_plot_raster_divergence_guard(tmp_path):
    raster_path, _, _ = _create_test_raster(tmp_path)

    with pytest.raises(ValueError):
        plot_raster_on_world_extremes_cutoff(
            raster_path,
            "Test",
            perc_cutoff=0,
            plt_show=False,
            base_shp=gpd.GeoDataFrame(),
            divergence_center=0,
            diverg0=True,
        )


def test_plot_n_rasters_returns_figures(tmp_path):
    raster_path, _, _ = _create_test_raster(tmp_path)

    outputs = plot_n_rasters_on_world_extremes_cutoff(
        [raster_path, raster_path],
        titles=["A", "B"],
        perc_cutoff=0,
        base_shp=gpd.GeoDataFrame(),
        plt_show=False,
    )

    assert len(outputs) == 2
    for fig, ax in outputs:
        assert fig is not None and ax is not None
        plt.close(fig)


def test_plot_n_rasters_title_length_guard(tmp_path):
    raster_path, _, _ = _create_test_raster(tmp_path)

    with pytest.raises(ValueError):
        plot_n_rasters_on_world_extremes_cutoff(
            [raster_path, raster_path],
            titles=["Only one"],
            perc_cutoff=0,
            base_shp=gpd.GeoDataFrame(),
            plt_show=False,
        )
