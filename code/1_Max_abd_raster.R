# ------------------------------------------------------------------------------

#             R SCRIPT FOR CALCULATING MAX WEEKLY PROPORTION 
#                   OF EASTERN U.S. POPULATION

# ------------------------------------------------------------------------------

# Created by Matthew Gilbert on 22 Oct, 2024
# Last modified by Matthew Gilbert on 22 Oct, 2024

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


# ------------------------------------------------------------------------------

#  "species_list" must be a list with the species codes for all of the 
#  mapped bird species which have ANY modeled population in the Eastern United 
#  states (as defined by the states above) during ANY time of year

example_list <- list("norpar", "magwar", "babwar", "bkbwar", "yelwar")

species_list <- example_list

# ------------------------------------------------------------------------------

for (species in species_list) {
  
  # processing message
  common_name <- ebirdst_runs$common_name[ebirdst_runs$species_code == species]
  cat("Processing:", common_name, "\n")
  
  # download weekly abundance data
  set_ebirdst_access_key("vhav03jqtha1", overwrite = TRUE)
  ebirdst_download_status(species = species, download_abundance = TRUE, pattern = "_3km_")
  abd_weekly_median <- load_raster(species, product = "abundance", 
                                   period = "weekly", metric = "median", 
                                   resolution = "3km")
  
  # crop data to the Eastern US
  east_rp <- st_transform(east, crs(abd_weekly_median))
  abd_weekly_median_east <- mask(crop(abd_weekly_median, vect(east_rp)),east_rp)
  abd_weekly_median_merc <- project(abd_weekly_median_east, crs) |> trim()
  east_merc <- st_transform(east_rp, crs)

  # calculate percent of max abundance
  max_abd <- app(abd_weekly_median_merc, max, na.rm = TRUE)
  sum_max_abd <- global(max_abd, fun = "sum", na.rm = TRUE)[1, 1]
  percent_max_abd <- max_abd / sum_max_abd

  # save raster as a file
  output_file <- paste0(species, "_percent_max_abd.tif")
  writeRaster(percent_max_abd, output_file, overwrite = TRUE)
  cat("Saved:", output_file, "\n")
}

cat("All species processed \n")








# ------------------------------------------------------------------------------

#       SCRIPT TO CALCULATE PROPORTION OF GLOBAL POPULATION IN EASTERN U.S.

# ------------------------------------------------------------------------------

#  "species_list" must be a LIST with the species codes for all of the 
#  mapped bird species which have ANY modeled population in the Eastern United 
#  states (as defined by the states above) during ANY time of year

data <- read.csv("species_list_10.29.24.csv", header = T)

data_list <- as.list(data$species_code)

example_list <- list("antnig", "monpar", "tenwar", "bkbwar", "yelwar")

species_list <- data_list


# ------------------------------------------------------------------------------

df <- data.frame(species_code = unlist(species_list))

for (i in seq_along(species_list)) {
  
  # download weekly abundance data
  species <- as.character(species_list[i])
  set_ebirdst_access_key("vhav03jqtha1", overwrite = TRUE)
  ebirdst_download_status(species = species, download_abundance = TRUE, pattern = "_27km_")
  abd_weekly_median <- load_raster(species, product = "abundance", 
                                   period = "weekly", metric = "median", 
                                   resolution = "27km")
  
  # calculate percent of max abundance
  max_abd <- app(abd_weekly_median, max, na.rm = TRUE)
  sum_max_abd <- global(max_abd, fun = "sum", na.rm = TRUE)[1, 1]
  percent_max_abd <- max_abd / sum_max_abd
  
  # crop data and calculate proportion
  east_rp <- st_transform(east, crs(percent_max_abd))
  percent_max_abd_east <- mask(crop(percent_max_abd, vect(east_rp)),east_rp)
  proportion <- global(percent_max_abd_east, fun = "sum", na.rm = TRUE)[1, 1]
  
  df$proportion_in_east[i] <- proportion

}

write.csv(df, "percent_max_abd_in_east.csv", row.names = FALSE)













