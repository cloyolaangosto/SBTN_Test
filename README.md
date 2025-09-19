# SBTN Land's Land Environmental Assessment Factors (LEAFs)
Land Environmental Assessment Factors (LEAF) are numerical factors used in SBTN Land v2 methods to help companies assess their direct operations impacts on land ecosystems (link). LEAF cover 3 categories: Soil Organic Carbon (SOC), Soil Erosion, and Terrestrial Acidification.

This repository contains 1) readily available LEAF for the most common land uses, and 2) instructions and basic code to derive new LEAF for each of the 3 categories.

## Available LEAF
LEAF are available in different formats and geographic levels of aggregation. Country, sub-country, and ecoregion average LEAF are available in Excel and Shapefiles. Unaggregated LEAF are available in their own original resolution as GeoTIFF files.

They can be found under the folder Available_LEAF in their respective category, containing excel, shp, and tiff folders. They have been generated using the method described below (AND MAYBE PAPER INSERT HERE).

## New LEAF Development
Documentation is provided to create LEAF for new practice changes is available only for SOC and Soil Erosion, as practice changes that affect terrestrial acidification affect emissions which already have availalbe LEAF.

For SOC and Soil Erosion, new LEAF generation is organized in 3 steps: data gathering and processing, data harmonization (same for both), and LEAF calculation. 

SOC documentation can be found [here](SOC_Documentation.md) and Soil Erosion [here](Soil_Erosion_Documentation.md).