"""sbtn_leaf: Environmental modelling utilities for SBTN-LEAF projects.

The top-level package intentionally re-exports the most commonly used
workflows so downstream notebooks can depend on a stable, well-documented
surface area.  The curated groups are:

* Potential evapotranspiration helpers (:mod:`sbtn_leaf.PET`)
  - :func:`calculate_PET_location_based`
  - :func:`calculate_PET_crop_based`
  - :func:`calculate_crop_based_PET_raster_optimized`
  - :func:`calculate_crop_based_PET_raster_vPipeline`
  - :func:`create_KC_Curve`
  - :func:`monthly_KC_curve`
* RothC core routines (:mod:`sbtn_leaf.RothC_Core`)
  - :class:`CarbonPools`
  - :func:`initialize_pools`
  - :func:`onemonth_step_rothc`
  - :func:`run_equilibrium`
  - :func:`run_simulation`
  - :func:`RMF_Tmp`
  - :func:`RMF_Moist`
  - :func:`RMF_PC`
  - :func:`RMF_TRM`
* Raster RothC pipelines (:mod:`sbtn_leaf.RothC_Raster`)
  - :func:`load_single_band`, :func:`load_multiband`
  - :func:`align_and_resample`, :func:`mask_by_landuse`, :func:`stack_time_series`
  - :func:`build_pc_mask`
  - :func:`write_single_band_tif`, :func:`write_multiband_tif`
  - :func:`raster_rothc_annual_only`
  - :func:`raster_rothc_annual_results`
  - :func:`save_annual_results`
  - :func:`run_RothC_crops`, :func:`run_RothC_forest`, :func:`run_RothC_grassland`
  - :func:`run_rothC_sceneraios_from_csv`
  - :func:`calcuate_annual_perc_changes`, :func:`calcuate_practice_change_benefit`
* Data lookups (:mod:`sbtn_leaf.data_loader`)
  - :func:`get_crop_coefficients_table`, :func:`get_absolute_day_table`
  - :func:`get_days_in_month_table`, :func:`get_crop_naming_index_table`
  - :func:`get_fao_statistics_table`, :func:`get_fao_crop_yields_table`
  - :func:`get_country_boundaries`, :func:`get_ecoregions_shapefile`
  - :func:`get_crop_ag_residue_table`, :func:`get_crop_residue_ratio_table`
  - :func:`get_thermal_climate_tables`
* Crop calculations (:mod:`sbtn_leaf.cropcalcs`)
  - :func:`index_files`
  - :func:`create_crop_yield_shapefile`, :func:`create_crop_yield_raster`
  - :func:`create_crop_yield_raster_withIrrigationPracticeScaling`
  - :func:`calculate_SPAM_yield_modifiers`
  - :func:`calculate_average_yield_by_ecoregion_and_biome`
  - :func:`calculate_crop_residues`
  - :func:`apply_residues_to_raster_flexible`, :func:`create_residue_raster_rasterops`
  - :func:`create_plant_cover_monthly_curve`, :func:`create_plant_cover_monthly_raster`
  - :func:`convert_K2C_raster`
  - :func:`compute_residue_raster`, :func:`compute_monthly_residue_raster`
  - :func:`compute_monthly_residue_raster_fromAnnualRaster`
  - :func:`calculate_irrigation_fromArray`, :func:`calculate_irrigation_fromTif`
  - :func:`prepare_crop_data`, :func:`prepare_crop_scenarios`
  - :func:`binarize_raster_pipeline`
  - :func:`calculate_irrigation_vPipeline`
  - :func:`create_crop_yield_raster_withIrrigationPracticeScaling_vPipeline`
  - :func:`create_monthly_residue_vPipeline`
  - :func:`get_forest_litter_rate`, :func:`get_forest_litter_monthlyrate_fromda`
  - :func:`generate_grassland_residue_map`, :func:`calculate_carbon_dung`
* Map calculations (:mod:`sbtn_leaf.map_calculations`)
  - :func:`calculate_area_weighted_cfs_from_shp_with_std_and_median`
  - :func:`calculate_area_weighted_cfs_from_raster_with_std_and_median`
  - :func:`calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers`
  - :func:`run_diagnostic`
  - :func:`build_cfs_gpkg_from_rasters`
  - :func:`rasterize_shapefile_to_target_raster`
  - :func:`calculate_polygon_means_from_raster`
  - :func:`multiply_rasters`, :func:`create_binary_mask`
  - :func:`resample_raster_to_match`, :func:`resample_to_match`, :func:`resample_to_match_multiband`, :func:`resample_to_match_noSaving`
  - :func:`subtract_rasters_union`
  - :func:`mask_raster1_by_overlap_with_raster2`
  - :func:`calculate_average_from_raster_timeseries`
  - :func:`calculate_average_from_raster_timeseries_byblocks`
  - :func:`binarize_raster`
  - :func:`extract_coordinates_values_from_shapefile`
  - :func:`extract_coordinates_values_from_raster`
* Plotting helpers (:mod:`sbtn_leaf.map_plotting`)
  - :func:`plot_raster_on_world_extremes_cutoff`
  - :func:`plot_da_on_world_extremes_cutoff`
  - :func:`plot_all_raster_bands`
  - :func:`plot_raster_on_world_no_min`
  - :func:`plot_static_shapefile_on_world`
  - :func:`plotly_shapefile_continuous`, :func:`plotly_shapefile_categorical`
  - :func:`plot_raster_over_gdf`, :func:`plot_raster_over_gdf_showpolygonvalues`
  - :func:`calculate_map_center`, :func:`preprocess_gdf`, :func:`filter_by_country`
  - :func:`format_hover_data`, :func:`determine_range_color`
  - :func:`downsample_raster`, :func:`inspect_raster`, :func:`get_raster_band_count`
  - :func:`extract_raster_values_from_tif`
  - :func:`plot_raster_data_histogram`, :func:`plot_raster_histogram`
  - :func:`plot_overlapping_histograms`
  - :func:`plot_multiband_raster_timesires`
  - :func:`plot_average_band_values`
  - :func:`plot_soc_distribution`

Refer to the module documentation for detailed usage patterns and expected
inputs for each routine.
"""

from .PET import (
    calculate_PET_crop_based,
    calculate_PET_location_based,
    calculate_crop_based_PET_raster_optimized,
    calculate_crop_based_PET_raster_vPipeline,
    create_KC_Curve,
    monthly_KC_curve,
)
from .RothC_Core import (
    CarbonPools,
    RMF_Moist,
    RMF_PC,
    RMF_TRM,
    RMF_Tmp,
    initialize_pools,
    onemonth_step_rothc,
    run_equilibrium,
    run_simulation,
)
from .RothC_Raster import (
    align_and_resample,
    build_pc_mask,
    calcuate_annual_perc_changes,
    calcuate_practice_change_benefit,
    load_multiband,
    load_single_band,
    mask_by_landuse,
    raster_rothc_annual_only,
    raster_rothc_annual_results,
    run_RothC_crops,
    run_RothC_forest,
    run_RothC_grassland,
    DEPRECATED_run_rothC_crop_scenarios_from_csv,
    save_annual_results,
    stack_time_series,
    write_multiband_tif,
    write_single_band_tif,
)
from .cropcalcs import (
    apply_residues_to_raster_flexible,
    _binarize_raster_pipeline,
    _calculate_SPAM_yield_modifiers,
    calculate_average_yield_by_ecoregion_and_biome,
    calculate_carbon_dung,
    calculate_crop_residues,
    calculate_irrigation_fromArray,
    calculate_irrigation_fromTif,
    calculate_irrigation_vPipeline,
    compute_monthly_residue_raster,
    compute_monthly_residue_raster_fromAnnualRaster,
    compute_residue_raster,
    convert_K2C_raster,
    create_crop_yield_raster,
    create_crop_yield_raster_withIrrigationPracticeScaling,
    create_crop_yield_raster_with_irrigation_scaling_pipeline,
    create_crop_yield_shapefile,
    create_monthly_residue_vPipeline,
    create_plant_cover_monthly_curve,
    create_plant_cover_monthly_raster,
    create_residue_raster_rasterops,
    generate_grassland_residue_map,
    get_forest_litter_monthlyrate_fromda,
    get_forest_litter_rate,
    index_files,
    prepare_crop_data,
    prepare_crop_scenarios,
)
from .data_loader import (
    get_absolute_day_table,
    get_country_boundaries,
    get_crop_ag_residue_table,
    get_crop_coefficients_table,
    get_crop_naming_index_table,
    get_crop_residue_ratio_table,
    get_days_in_month_table,
    get_ecoregions_shapefile,
    get_fao_crop_yields_table,
    get_fao_statistics_table,
    get_thermal_climate_tables,
)
from .map_calculations import (
    binarize_raster,
    build_cfs_gpkg_from_rasters,
    calculate_area_weighted_cfs_from_raster_with_std_and_median,
    calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers,
    calculate_area_weighted_cfs_from_shp_with_std_and_median,
    calculate_average_from_raster_timeseries,
    calculate_average_from_raster_timeseries_byblocks,
    calculate_polygon_means_from_raster,
    create_binary_mask,
    extract_coordinates_values_from_raster,
    extract_coordinates_values_from_shapefile,
    mask_raster1_by_overlap_with_raster2,
    multiply_rasters,
    rasterize_shapefile_to_target_raster,
    resample_raster_to_match,
    resample_to_match,
    resample_to_match_multiband,
    resample_to_match_noSaving,
    run_diagnostic,
    subtract_rasters_union,
)
from .map_plotting import (
    _calculate_map_center,
    _determine_range_color,
    downsample_raster,
    extract_raster_values_from_tif,
    _filter_by_country,
    _format_hover_data,
    get_raster_band_count,
    inspect_raster,
    plot_all_raster_bands,
    plot_average_band_values,
    plot_da_on_world_extremes_cutoff,
    plot_multiband_raster_timesires,
    plot_overlapping_histograms,
    plot_raster_data_histogram,
    plot_raster_histogram,
    plot_raster_on_world_extremes_cutoff,
    plot_raster_on_world_no_min,
    plot_raster_over_gdf,
    plot_raster_over_gdf_showpolygonvalues,
    plot_soc_distribution,
    plot_static_shapefile_on_world,
    _preprocess_gdf,
)


_PET_EXPORTS = [
    "calculate_PET_crop_based",
    "calculate_PET_location_based",
    "calculate_crop_based_PET_raster_optimized",
    "calculate_crop_based_PET_raster_vPipeline",
    "create_KC_Curve",
    "monthly_KC_curve",
]
_ROTHC_CORE_EXPORTS = [
    "CarbonPools",
    "RMF_Moist",
    "RMF_PC",
    "RMF_TRM",
    "RMF_Tmp",
    "initialize_pools",
    "onemonth_step_rothc",
    "run_equilibrium",
    "run_simulation",
]
_ROTHC_RASTER_EXPORTS = [
    "align_and_resample",
    "build_pc_mask",
    "calcuate_annual_perc_changes",
    "calcuate_practice_change_benefit",
    "load_multiband",
    "load_single_band",
    "mask_by_landuse",
    "raster_rothc_annual_only",
    "raster_rothc_annual_results",
    "run_RothC_crops",
    "run_RothC_forest",
    "run_RothC_grassland",
    "run_rothC_sceneraios_from_csv",
    "save_annual_results",
    "stack_time_series",
    "write_multiband_tif",
    "write_single_band_tif",
]
_CROPCALCS_EXPORTS = [
    "apply_residues_to_raster_flexible",
    "binarize_raster_pipeline",
    "calculate_SPAM_yield_modifiers",
    "calculate_average_yield_by_ecoregion_and_biome",
    "calculate_carbon_dung",
    "calculate_crop_residues",
    "calculate_irrigation_fromArray",
    "calculate_irrigation_fromTif",
    "calculate_irrigation_vPipeline",
    "compute_monthly_residue_raster",
    "compute_monthly_residue_raster_fromAnnualRaster",
    "compute_residue_raster",
    "convert_K2C_raster",
    "create_crop_yield_raster",
    "create_crop_yield_raster_withIrrigationPracticeScaling",
    "create_crop_yield_raster_withIrrigationPracticeScaling_vPipeline",
    "create_crop_yield_shapefile",
    "create_monthly_residue_vPipeline",
    "create_plant_cover_monthly_curve",
    "create_plant_cover_monthly_raster",
    "create_residue_raster_rasterops",
    "generate_grassland_residue_map",
    "get_forest_litter_monthlyrate_fromda",
    "get_forest_litter_rate",
    "index_files",
    "prepare_crop_data",
    "prepare_crop_scenarios",
]
_DATA_LOADER_EXPORTS = [
    "get_absolute_day_table",
    "get_country_boundaries",
    "get_crop_ag_residue_table",
    "get_crop_coefficients_table",
    "get_crop_naming_index_table",
    "get_crop_residue_ratio_table",
    "get_days_in_month_table",
    "get_ecoregions_shapefile",
    "get_fao_crop_yields_table",
    "get_fao_statistics_table",
    "get_thermal_climate_tables",
]
_MAP_CALCULATIONS_EXPORTS = [
    "binarize_raster",
    "build_cfs_gpkg_from_rasters",
    "calculate_area_weighted_cfs_from_raster_with_std_and_median",
    "calculate_area_weighted_cfs_from_raster_with_std_and_median_vOutliers",
    "calculate_area_weighted_cfs_from_shp_with_std_and_median",
    "calculate_average_from_raster_timeseries",
    "calculate_average_from_raster_timeseries_byblocks",
    "calculate_polygon_means_from_raster",
    "create_binary_mask",
    "extract_coordinates_values_from_raster",
    "extract_coordinates_values_from_shapefile",
    "mask_raster1_by_overlap_with_raster2",
    "multiply_rasters",
    "rasterize_shapefile_to_target_raster",
    "resample_raster_to_match",
    "resample_to_match",
    "resample_to_match_multiband",
    "resample_to_match_noSaving",
    "run_diagnostic",
    "subtract_rasters_union",
]
_MAP_PLOTTING_EXPORTS = [
    "calculate_map_center",
    "determine_range_color",
    "downsample_raster",
    "extract_raster_values_from_tif",
    "filter_by_country",
    "format_hover_data",
    "get_raster_band_count",
    "inspect_raster",
    "plot_all_raster_bands",
    "plot_average_band_values",
    "plot_da_on_world_extremes_cutoff",
    "plot_multiband_raster_timesires",
    "plot_overlapping_histograms",
    "plot_raster_data_histogram",
    "plot_raster_histogram",
    "plot_raster_on_world_extremes_cutoff",
    "plot_raster_on_world_no_min",
    "plot_raster_over_gdf",
    "plot_raster_over_gdf_showpolygonvalues",
    "plot_soc_distribution",
    "plot_static_shapefile_on_world",
    "plotly_shapefile_categorical",
    "plotly_shapefile_continuous",
    "preprocess_gdf",
]


__all__ = (
    _PET_EXPORTS
    + _ROTHC_CORE_EXPORTS
    + _ROTHC_RASTER_EXPORTS
    + _CROPCALCS_EXPORTS
    + _DATA_LOADER_EXPORTS
    + _MAP_CALCULATIONS_EXPORTS
    + _MAP_PLOTTING_EXPORTS
    + ["__version__", "__author__"]
)

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
