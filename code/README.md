
SCRIPTS

1_Max_abd_raster.R
This script downloads the 2022 version of eBird Data Products and produces a map of year-round maximum relative abundance for each species for just the Eastern United States. All the species maps produced by this workflow can be found on our figshare repository in the species_max_abundance_maps folder. The second part of this script runs through a list of all species on ebirdst and calculates the proportion of their global maximum relative abundance which occurs in the Eastern U.S.

    INPUT: 
    None, but eBird S&T Data Access Key Needed

    OUTPUT: 
    Maximum year-round relative abundance maps for each species
    CSV file of the percent of global max abundance which occurs in the Eastern U.S. for each  
    species

2_Community_raster_builder.R
This script collects the species maximum abundance maps, aligns them with climate velocity and/or changes in LULC, calculates the CEI/LEI for each species, creates and saves species-specific impact maps for each species, and creates a community-wide exposure map for the selected scenario. The loop also outputs a CSV file with each species and its CEI and LEI for that scenario.

    INPUT: 
    Climate velocity rasters for selected scenario(s)
    Land-cover forecast raster for selected scenario(s)
    Current land cover raster
    Maximum year-round relative abundance maps for each species

    OUTPUT: 
    Raster of exposure to climate change across all selected species
    Raster of exposure to LULC change across all selected species
    Raster of species specific exposure map for climate for selected scenario
    Raster of species specific exposure map for LULC change for selected scenario
    CSV file with the CEI and LEI for each species

3_Impact_mapper.R
This script loads files calculated by the community_raster_builder.r to map total effect of climate and LULC change exposure on Eastern Landbirds, as well as map the species-specific effect of climate change.

    INPUT:
    Community-level Exposure Rasters for selected scenario
    Species-level exposure Rasters for selected scenario and species

    OUTPUT:
    Figure 1 (Species-level exposure map for CEI and LEI for Wood Thrush)
    Figure 3 (Community-level exposure map for CEI and LEI)

4_glm_model.R
Uses data frames calculated by the community_raster_builder.r script to analyze patterns in CEI between groups of species using a log-adjusted generalized linear model (GLM).

    INPUT: 
    CSV files of CEI and LEI for each species in each time step

    OUTPUT: 
    Figure 4 (model outputs for CEI/LEI by habitat type)
    Figure 5 (model outputs for CEI/LEI by latitude and longitude)
    Supplemental model diagnostic figures
    Table S2, S3, S4

4.1_diagnostic_tests.R
Compares 3 model families (gaussian linear, log-transformed linear, and gamma) for CEI and 2 model families (gaussian linear, beta) for LEI using DHARMa model diagnostics.

    INPUT: 
    CSV files of CEI and LEI for each species in each time step

    OUTPUT: 
    Plots of diagnostics figures (Figures S3-S10)


4.2_alt_model_w_order.R
Uses data frames calculated by the community_raster_builder.r script to analyze differences in climate/LULC change exposure between groups of species using a linear model. This shows results for an alternative model with delta_AIC<2 but not the best model overall. This model includes lat, long, habitat, and order.

    INPUT: 
    CSV files of CEI and LEI for each species in each time step

    OUTPUT: 
     Figures S12-S13: Predictive plots showing modeled relationships for alternate linear model


4.3_alt_model_w_guild.R
Uses data frames calculated by the community_raster_builder.r script to analyze differences in climate/LULC change exposure between groups of species using a linear model. This shows results for an alternative model with delta_AIC<2 but not the best model overall. This model includes lat, long, and guild.

    INPUT: 
    CSV files of CEI and LEI for each species in each time step

    OUTPUT: 
    Figures S14-S15: Predictive plots showing modeled relationships for alternate linear model

4.4_coefficient_estimates.R
Uses data frames calculated by the community_raster_builder.r script to calculate coefficient estimates for all models with substantial support. 

    INPUT: 
    CSV files of CEI and LEI for each species in each time step

    OUTPUT: 
    Table S5

5_appendix_figures.R
This script creates a handful of miscellaneous figures, as well as calculates some of the information relevant for each plot (i.e. R^2 values).

    INPUT:  
    CSV files of CEI and LEI for each species in each time step
    Community-level CEI/LEI exposure maps for all scenarios
    Species-level CEI/LEI exposure maps for selected species for all scenarios
    Climate velocity projections for each scenario
    LULC projections for each scenario

    OUTPUT: 
    Figure 2
    Table 1
    Figures S1-S2, S11, S16-S17


