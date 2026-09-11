# DATA
ASSESSING EXPOSURE OF EASTERN NORTH AMERICAN LANDBIRDS TO CLIMATE AND LAND USE CHANGE

All input and output data for this study is stored in separate respositories on Figshare due to storage and file number constraints. The full repository can be found here: 

_____

## Species List Data (input)
The full list of 218 species with >0.01% of maximum relative abundance globally is provided in the _species_list.csv_ file within this folder. For each species the data includes: common name, 6-letter species code, family, habitat, ecological guild, and the proportion in the east. 

## Abundance Data (input)
Abundance data is sourced from eBird Status and Trends, cropped to the bounds of the study area, and then normalized. Cropped rasters of maximum relative abundance for each species can be found on the Figshare repository under the **species_max_abundance_maps** folder. 

## Climate Data (input)
Climate velocity data is sourced from Carrol et al., 2023. We used four data files for forward-projected climate velocity, one for each time step and emissions scenario (SSP2 2041-2070, SSP5 2041-2080, SSP2 2071-2100, & SSP5 2071-2100). These can be found in the **climate_velocity_data** folder. 

## Land Use Data (input)
Future projects of land use are sourced from Chen et al., 2022. We used five data files for forward-projected climate velocity, one for each time step and emissions scenario (SSP2 2041-2070, SSP5 2041-2080, SSP2 2071-2100, & SSP5 2071-2100), and one basline file to compare it to (2015 data). These can be found in the **LULC_change_data** folder. 

## Community-level Exposure Maps (output)
Community-level (inclusive of all landbird species) exposure maps were calculated for each emissions scenario and time step combination, and can be found on Figshare in the **community_maps** folder. This includes 8 maps: 2 threats (LEI, CEI) x 2 socio-economic pathways x 2 time steps. 

## Species-level Exposure Maps (output)
Species-level exposure maps for each specieswere calculated for each emissions scenario and time step combination. CEI exposure maps can be found in **species_climate_impact_maps.zip** and LEI exposure maps can be found in **species_LULC_impact_maps.zip**. These zip files each include 872 maps: 218 species x 2 socio-economic pathways x 2 time steps. 

## Table S1 (output)
Table S1 in supplemental information is included in this folder at TableS1.csv. It includes all the information from the species_list.csv as well as the calculated values for CEI and LEI for each time step and emissions scenario. 
