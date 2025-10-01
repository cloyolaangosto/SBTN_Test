# Soil Erosion Documentation
Soil Erosion LEAFs are factors that estimate the annual rate of soil loss (t soil/year) for different land uses and land management practices. Despite the existence of other models, due to global applicability and data sources, they are currently based on the Revised Universal Soil Loss Equation (RUSLE, **INSERT CITATION HERE**). Guidance on how to apply LEAFs to different land uses has been covered in detail on SBTN Land's AGILED document [link](https://sciencebasedtargetsnetwork.org/wp-content/uploads/2025/04/SBTN-Land-Accounting-Guidelines-Draft-for-Public-Consultation.pdf).

LEAFs are available in [`LEAFs/soil_erosion`](../LEAFs) in csv format and shapefile for ecoregions, countries and subcountries. Native resolution (25km) raster files can be found in **XXX - Replace when we figure this out**. How to generate new LEAFs can be found on the notebook [New Leafs Soil Erosion Example](../examples/New_LEAFs_Soil_Erosion_Example.ipynb).

## RUSLE
The RUSLE equation estimates soil loss simply by:

$RE = R\cdot K_{st}\cdot LS \cdot C \cdot P$

Where,
- *$R$*: Rainfall erosivity factor $[MJ*mm/(ha*h*a)]$
- *$K_{st}$*: Soil erodibility factor taking into account the soil skeletal fraction $[t*ha*h/(MJ*ha*mm)]$
- *$LS$:* Slope length and slope steepness factor [ ]
- *$C$*: Ground cover and tillage factor [ ]
- *$P$*: Erosion protection factor [ ]

The variables can then be divided in soil and weather ($R$, $K_{st}$, and $LS$), and crop and land management relaed ($C$ and $P$). As soil and weather remain constant independent of the land use and land management option, a single $R\cdot K_{st}\cdot LS$ raster layer has been generated and used as a base to simplify future calculations. $C$ is then adjusted depending on crops and tillage factor, while $P$ has been assumed as 1, a common assumption on LCA characterization factor generation (**CITATION CITATION HERE**).

## Step 1 - Data Gathering and Processing
### Soil and Weather
The Global Soil Erosion Modelling platform ([GloSEM](https://esdac.jrc.ec.europa.eu/content/global-soil-erosion)) has released several estimations of global soil erosion, the most recent being v1.3, along with publications that are the base of SBTN LEAFs ([Borelli et al., 2017](https://www.nature.com/articles/s41467-017-02142-7), [Borelli et al., 2022](https://www.nature.com/articles/s41467-017-02142-7), [Panagos et al., 2015](https://www.sciencedirect.com/science/article/pii/S0264837715001611?via%3Dihub)).


v1.1 released publicly the $R$, $K_{st}$, and $LS$ layers at 25km resolution ([download link](https://esdac.jrc.ec.europa.eu/public_path/shared_folder/dataset/58/GloSEM_25km.zip)) and they serve as the base of these calculations. Those 3 layers were downloaded and multiplied. 

(![RLKSK Factors](../archives/se_rlsk.png))

### Crop Data
#### Land use
Land use maps for each commodity was sourced from [Morais, Teixeria & Domingos (2019)](https://doi.org/10.1371/journal.pone.0222604), which divided the globe into UHTH zones of X by Y resolution following FAO GAEZ v2 suitability index maps.

**[TODO - GET RESOLUTION OF UHTH]**

#### C-Factor - Commodities
Crop C-Factors were obtained from [Borelli et al., 2017](https://www.nature.com/articles/s41467-017-02142-7) according to the following tables:

Table 1 - C Values for IGBP Land-cover type (Source: Borelli et al., 2017)
|Class|IGBP Land-cover type|C-Factor|
|:---:|:---|:-----:|
|0| Water | No data |
|1| Evergreen Needleaf forest | 0.0001 - 0.003 |
|2| Evergreen Broadleaf forest | 0.0001 - 0.003 |
|3| Decidious Needleaf forest | 0.0001 - 0.003 |
|4| Decidious Broadleaf forest | 0.0001 - 0.003 |
|5| Mixed forest | 0.0001 - 0.003 |
|6| Closed shrublands | 0.01 - 0.15 |
|7| Open shrublands | 0.01 - 0.15 |
|8| Woody savannas | 0.01 - 0.15 |
|9| Savannas | 0.01 - 0.15 |
|10| Grasslands | 0.01 - 0.15 |
|11| Permanent wetlands | No data |
|13| Urban and built-up | No data |
|15| Snow and ice | No data |
|16| Barren or sparsely vegetated | 0.1 - 0.5 |
|21| Transitional woodland-shrub | 0.01 - 0.15|

Table 2 - C Values for IGBP Land-cover type  (Source: Borelli et al., 2017)
| Crop Group | Crop Sub-group | C-Factor |
|:---|:---|:-----:|
| Cereal Grains | Various | 0.20 |
| Cereal Grains | Maize | 0.20 |
| Cereal Grains | Rice | 0.20 |
| Legume Vegetables | Various | 0.32 |
| Root and Tuber Vegetables | Various | 0.34 |
| Fruiting Vegetables| Various | 0.25 |
| Cucurbit Vegetables | Various | 0.25 |
| Bulby Vegetables | Various | 0.30 |
| Leafy Vegetables | Various | 0.25 |
| Leafy Vegetables | Tobacco | 0.50 |
| Forage, Fodder, and Straw of Cereal Grain Group | Mixed-legumes | 0.15 |
| Forage, Fodder, and Straw of Cereal Grain Group | Mixed-grasses | 0.10 |
| Grain and Hops | Grains| 0.35 |
| Grain and Hops | Hops | 0.42 |
| Oilseed Group | Various | 0.25 |
| Oilseed Group | Cotton | 0.40 |
| Fibre Crops | Fibre Crops | 0.28 |
| Berries Group | Various | 0.15 |
| Berries Group | Strawberries | 0.20 |
| Shrubs Herbs and Spices Group | Shrubs Herbs and Spices | 0.15 |
| Shrubs Herbs and Spices Group | Coffee | 0.20 |
| Trees/Fruit Trees| Various | 0.15 |

Each land use was mapped to a C-Factor from both tables according to FAO commodity classification, and can be found at the [`lu_c_factor_inventory.xlsx`](../data/crops/lu_c_factor_inventory.xlsx) file in `data/crops/`.

After assigning C-Factors to each commodity, a C-Factor raster was created for each by simply assigning that value for each pixel where the commodity can be produced. This factors can **BE OBTAINED HERE WHEN WE DECIDE WHERE TO UPLOAD THEM**

#### C-Factors - Land Management: Tillage, Residue Management and Cover Crops
According to Panagos et al. (2015), C-Factors for crops can be further dissagregted using the following equations:

$C_{arable} = C_{crop} \times C_{mgmt}$

Where $C_{mgmg}$ represents the C-Factor of land management options including tillage, residues, and cover crops. This is also expressed by a simple multiplication:

$C_{mgmt} = C_{tillage} \times C_{residues} \times C_{cover}$

For **tillage**, a value of 1, 0.35, or 0.25 is assigned respectively to convetional tillage, reduced/conservation tillage, and No-Till practices.

For **residue management**, it can be calculated as:

$C_{residues} = 1 - 0.12 \cdot F_{residues}$

Where $F_{residues} is the fraction of land that leaves residues on the field.

Finally, for **cover crops**, $C_{cover}$ can be calculated as:

$C_{cover} = 1 - 0.2 \cdot F_{cover}$

Where $F_{cover} is the fraction of land that uses cover crops during winter or spring.

## Step 2 - Data Harmonization
All crops C-Factors raster were downsampled into 25km values using the *RLSK* raster as base, using a nearest value approach.

No other data harmonization was needed.

## Step 3 - LEAFs Calculations

## Step 4 - LEAFs Averages for Ecoregions, Countries and Subcountries
