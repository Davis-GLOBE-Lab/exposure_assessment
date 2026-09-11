
# ------------------------------------------------------------------------------

#               SCRIPT FOR ANALYZING COMMUNITY EXPOSURE

# ------------------------------------------------------------------------------

# This script collects the species maximum abundance maps, aligns them with climate 
# velocity and/or changes in LULC, calculates the CEI/LEI for each species, creates 
# and saves species-specific impact maps for each species, and creates a community-
# wide exposure map for the selected scenario. The loop also outputs a CSV file with 
# each species and it's CEI and LEI for that scenario. 


# Last modified by: Matthew Gilbert
# Last modified on: 15 Mar 2025


# --------------------------        SETUP        -------------------------------


library(terra)
library(ebirdst)
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(lubridate)
library(sf)
library(patchwork)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
# if struggling with install for naturalearthhires download, use this command:
# install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")


unitedstates <- rnaturalearth::ne_states(country = "United States of America", 
                                         returnclass = "sf")
east <- unitedstates[unitedstates$region == "Northeast" |
                       unitedstates$name == "Delaware" |
                       unitedstates$name == "Kentucky" | 
                       unitedstates$name == "West Virginia" | 
                       unitedstates$name == "Virginia" |
                       unitedstates$name == "Ohio" |
                       unitedstates$name == "North Carolina" |
                       unitedstates$name == "South Carolina" |
                       unitedstates$name == "Georgia" |
                       unitedstates$name == "Florida" |
                       unitedstates$name == "Alabama" |
                       unitedstates$name == "Mississippi" |
                       unitedstates$name == "Indiana" |
                       unitedstates$name == "Michigan" |
                       unitedstates$name == "Illinois" |
                       unitedstates$name == "Wisconsin" |
                       unitedstates$name == "Tennessee" |
                       unitedstates$name ==  "Maryland" ,]
crs <- "EPSG:3857"

setwd("~/Desktop/vulnerability_analysis")


# -----------   FILTERING, THRESHOLDS, AND CHOOSING DATE/SSP     ---------------


#  "species_list" must be a list with the species codes for all of the 
#  mapped bird species which have ANY modeled population in the Eastern United 
#  states (as defined by the states above) during ANY time of year

# load in species data frame
species_data <- read.csv("species_list.csv", header=T)

# filter to only eBird products from 2022
year <- 2022
species_data_current <- species_data[species_data$year == year, ]

# filter to only species with >0.1% of their global max abundance in the study area
#       (this value is pre-calculated using the Max_abd_raster.r script)
filtered_species_data <- species_data_current %>%
  filter(proportion_in_east >= 0.001)

# set as current data frame
species_df_1 <- filtered_species_data

# SET TIME STEP as "2050's" or "2080's"
time_step <- "2050's"

# SET CLIMATE SCENARIO as "ssp2" or ssp5"
climate_scenario <- "ssp5"


# --------  FUNCTION FOR CALCULTING Exposure and MAKING MAPS  -------------


CEI_weighted_map <- function(species_df, time_step, climate_scenario, east) {
  # Set climate scenario variables
  if (climate_scenario == "ssp5") {
    climate_scenario_num <- "585"
    climate_scenario_string <- "SSP5_RCP85"
  } else if (climate_scenario == "ssp2") {
    climate_scenario_num <- "245"
    climate_scenario_string <- "SSP2_RCP45"
  } else {
    stop("Invalid climate scenario: choose ssp2 or ssp5")
  }
  
  # Set time step variables
  if (time_step == "2050's") {
    time_step_num <- "2055"
    time_step_range <- "2041_2070"
  } else if (time_step == "2080's") {
    time_step_num <- "2085"
    time_step_range <- "2071_2100"
  } else {
    stop("Invalid time step: choose 2050's or 2080's")
  }
  
  # Helper function: Crop raster to the extent of 'east'
  crop_to_east <- function(raster) {
    east_bad <- st_transform(east, crs(raster))
    raster_east <- mask(crop(raster, vect(east_bad)),east_bad)
    raster <- project(raster_east, crs) |> trim()
    east <- st_transform(east_bad, crs)
    return(raster)
  }
  
  # Load and crop current LULC data
  current_LULC <- rast("LULC_change_data/global_LULC_2015.tif")
  current_LULC <- crop_to_east(current_LULC)
  
  # Load, crop, and align climate data
  climate_file <- paste0(
    "climate_velocity_data/fwvel731_ensemble_8gcm_",
    climate_scenario_num, "_", time_step_range, ".tif"
  )
  climate_velocity <- rast(climate_file)
  climate_velocity <- crop_to_east(climate_velocity)
  
  # Load and crop projected LULC data
  lulc_file <- paste0(
    "LULC_change_data/global_",
    climate_scenario_string, "_", time_step_num, ".tif"
  )
  projected_LULC <- rast(lulc_file)
  projected_LULC <- crop_to_east(projected_LULC)
  
  # Calculate LULC change raster
  LULC_change <- ifel(projected_LULC != current_LULC, 1, 0)
  
  # Initialize community maps
  community_map_climate <- NULL
  community_map_LULC <- NULL
  
  # Loop through species data frame
  for (i in seq_len(nrow(species_df))) {
    species_code <- species_df$species_code[i]
    species_name <- species_df$species_name[i]
    
    # Load, crop, and align species abundance map
    abundance_file <- paste0(
      "species_max_abundance_maps/", 
      species_code, "_abundance.tif"
    )
    max_abundance <- rast(abundance_file)
    max_abundance <- crop_to_east(max_abundance)
    
    # Convert to percent max abundance
    sum_max_abd <- global(max_abundance, fun = "sum", na.rm = TRUE)[1, 1]
    percent_max_abd <- max_abundance / sum_max_abd
    
    # Calculate climate product map
    climate_velocity <- resample(climate_velocity, percent_max_abd)
    climate_product_map <- percent_max_abd * climate_velocity
    CEI <- global(climate_product_map, fun = "sum", na.rm = TRUE)[1, 1]
    
    # Calculate LULC product map
    LULC_change <- resample(LULC_change, percent_max_abd)
    LULC_product_map <- percent_max_abd * LULC_change
    prop_LULC_change <- global(LULC_product_map, fun = "sum", na.rm = TRUE)[1, 1]
    
    # Update species data frame with results
    species_df$CEI[i] <- CEI
    species_df$prop_LULC_change[i] <- prop_LULC_change
    
    # Load the terra package for raster operations
    library(terra)
    
    # If the CRS is projected (e.g., UTM), reproject to a geographic CRS (e.g., WGS84)
    percent_max_abd_geo <- project(percent_max_abd, "EPSG:4326")  # EPSG:4326 corresponds to WGS84
    
    # Use the reprojected raster for further calculations
    longitude_raster <- init(percent_max_abd_geo, "x")
    latitude_raster <- init(percent_max_abd_geo, "y")
    
    # Calculate mean longitude
    weighted_sum_longitude <- as.numeric(global(longitude_raster * percent_max_abd_geo, fun = "sum", na.rm = TRUE)[1, 1])
    total_abundance <- as.numeric(global(percent_max_abd_geo, fun = "sum", na.rm = TRUE)[1, 1])
    weighted_avg_longitude <- weighted_sum_longitude / total_abundance
    species_df$mean_longitude[i] <- weighted_avg_longitude
    
    # Calculate mean latitude
    weighted_sum_latitude <- as.numeric(global(latitude_raster * percent_max_abd_geo, fun = "sum", na.rm = TRUE)[1, 1])
    weighted_avg_latitude <- weighted_sum_latitude / total_abundance
    species_df$mean_latitude[i] <- weighted_avg_latitude
    
    # Update community maps
    climate_impact_map <- percent_max_abd * CEI
    community_map_climate <- if (is.null(community_map_climate)) {
      climate_impact_map
    } else {
      community_map_climate + climate_impact_map
    }
    
    LULC_impact_map <- percent_max_abd * prop_LULC_change
    community_map_LULC <- if (is.null(community_map_LULC)) {
      LULC_impact_map
    } else {
      community_map_LULC + LULC_impact_map
    }
    
    # Save individual species impact maps
    climate_output <- paste0(
      "species_climate_impact_maps/", 
      species_code, ".", time_step, ".", climate_scenario, ".tif"
    )
    writeRaster(climate_product_map, climate_output, overwrite = TRUE)
    
    LULC_output <- paste0(
      "species_LULC_impact_maps/", 
      species_code, ".", time_step, ".", climate_scenario, ".tif"
    )
    writeRaster(LULC_product_map, LULC_output, overwrite = TRUE)
    
    message("Calculation Completed for ", species_name)
  }
  
  # Save final community map for Climate
  community_climate_output <- paste0("community_maps/", 
                                     "community_map_climate", "_", 
                                     time_step, "_", climate_scenario, ".tif"
  )
  writeRaster(community_map_climate, community_climate_output, overwrite = TRUE)
  
  # Save final community map for LULC
  community_LULC_output <- paste0("community_maps/",
                                  "community_map_LULC", "_", 
                                  time_step, "_", climate_scenario, ".tif"
  )
  writeRaster(community_map_LULC, community_LULC_output, overwrite = TRUE)
  
  # Save final community spreadsheet for CEI and prop_LULC_change
  df_output_name <- paste0("output_data_frames/", 
                           "output_df", "_", 
                           time_step, "_", climate_scenario, ".csv"
  )
  
  write.csv(species_df, df_output_name)
  
  message("Function completed successfully. Files saved to folder")
  
  # Return community maps
  return(list(
    community_map_climate <- community_map_climate,
    community_map_LULC <- community_map_LULC, 
    output_df <- species_df
  ))
}


# --------------------    FUNCTION CALL and OUTPUT    --------------------------


community_rasters <- CEI_weighted_map(species_df_1, time_step, climate_scenario, east)

community_map_climate <- community_rasters[[1]]  # First raster
community_map_LULC <- community_rasters[[2]]     # Second raster
output_df <- community_rasters[[3]]