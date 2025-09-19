### SOC
SOC LEAF are generated following the methods described in [Morais, Teixeria & Domingos (2019)](https://doi.org/10.1371/journal.pone.0222604) and [Teixeira, Morais & Domingos (2021)](https://doi.org/10.1038/s41597-021-01018-2), these in turn use the Rothmadel C (RothC) model to simulate SOC changes through time. RothC simulates "the turnover of organic carbon in non-waterlogged top-soils that allows for the effects of soil type, temperature, moisture content and plant cover on the turnover process" (reference). Complete documentation can be found on LINK.

The model has been adaptaed into GIS format, using as a base the one found on XXX.

#### Data Gathering and Processing
The model needs several inputs, despite being the most simple, that ara availalbe globally. A data download support script is availalbe in XXX.

##### Data Download
Weather related data has been obtained from NASA projects, while soil related data has been sourced from SoilGrids. 

###### 1) Monthly Precipitation (mm/month)
For the newly generated practice change LEAFs, [NASA's GPM_3IMERGM](https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGM_07/summary?keywords=3IMERG) data, between 2013 and 2023 and given at 0.1°, has been used.

###### 2) Temperature (mm)
For the newly generated practice change LEAFs, [GLDAS Catchment Land Surface Model L4](https://disc.gsfc.nasa.gov/datasets/GLDAS_CLSM10_M_2.1/summary), monthly 1.0 x 1.0 degree V2.1 (GLDAS_CLSM10_M), at a 1° resolution has been used.

###### 3) Soil Organic Carbon (t C/ha)
Organic Carbon Stock (0-30 cm) mean has been used. [link](https://files.isric.org/soilgrids/latest/data/ocs/) 

###### 4) Clay Content (percentage)
Average clay content (15-30 cm) has beeb used. [link](https://files.isric.org/soilgrids/latest/data/clay/) 

###### 5) Sand Content, only needed for reduced tillage (percentage)
Average sand content (15-30 cm) has beeb used. [link](https://files.isric.org/soilgrids/latest/data/sand/)

###### 6) Clay Content (percentage)

