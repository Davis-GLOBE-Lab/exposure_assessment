# ------------------------------------------------------------------------------

#                       SCRIPT for PLOTTING IMPACT MAPS

# this script loads files calculated by the community_raster_builder.r to map 
# total effect of climate and LULC change vulnerability on Eastern Landbirds, 
# as well as maps the species-specific effect of climate change

# Last modified by: Matthew Gilbert
# Last modified on: 5 Aug 2026

library(cowplot)
library(viridis)
library(ggplotify)
library(terra)
library(png)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(gridGraphics)
library(rnaturalearthhires)
library(grid)
library(gridExtra)

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
east <- st_transform(east, crs)

setwd("~/Desktop/vulnerability_analysis")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


par(mar = c(2, 2, 2, 2))
plot(
  x = seq_along(pal),
  y = rep(2, length(pal)),
  col = pal,
  pch = 15,
  cex = 8,
  axes = FALSE,
  xlab = "",
  ylab = ""
)


# ----------------------    COMMUNITY SETUP    ---------------------------------


# SET TIME STEP as "2050's" or "2080's"
time_step <- "2080's"

# SET CLIMATE SCENARIO as "ssp2" or ssp5"
climate_scenario <- "ssp5"

# create file paths
climate_file_path <- file.path("community_maps", paste0("community_map_climate_", time_step, "_", climate_scenario, ".tif"))
LULC_file_path <- file.path("community_maps", paste0("community_map_LULC_", time_step, "_", climate_scenario, ".tif"))

# load in files
community_map_climate <- rast(climate_file_path)
community_map_LULC <- rast(LULC_file_path)




# ---------------    MAPPING COMMUNITY-WIDE CLIMATE IMPACT   -------------------

raster_values <- values(community_map_climate, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community Climate Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_CEI_map.png", width = 4000, height = 4000, res = 1000)

plot(community_map_climate, 
     #main = title, 
     cex.main = 1.5, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)

dev.off()  # Closes the device




# ---------------    MAPPING COMMUNITY-WIDE LULC IMPACT   ----------------------

raster_values <- values(community_map_LULC, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community LULC Change Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_LEI_map.png", width = 4000, height = 4000, res = 1000)

plot(community_map_LULC, 
     #main = title, 
     cex.main = 1.5, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)

dev.off()  # Closes the device






# ------------------------------------------------------------------------------

# Load images
CEI_map_image <- rasterGrob(readPNG("plots/community_CEI_map.png"), interpolate = TRUE)
LEI_map_image <- rasterGrob(readPNG("plots/community_LEI_map.png"), interpolate = TRUE)
legend_image <- rasterGrob(readPNG("plots/vulnerability_map_legend.png"), interpolate = TRUE)

# Create labels "A" and "B"
CEI_with_label <- ggdraw() +
  draw_grob(CEI_map_image) +
  draw_label("a", x = 0.1, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate", x = 0.4, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold")

LEI_with_label <- ggdraw() +
  draw_grob(LEI_map_image) +
  draw_label("b", x = 0.1, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Land-Use", x = 0.4, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold")

# Combine the CEI and LEI maps side by side
side_by_side <- plot_grid(CEI_with_label, LEI_with_label, nrow = 1, rel_widths = c(1, 1))

# Overlay the legend on top of the side-by-side maps
final_plot <- ggdraw() +
  draw_plot(side_by_side, x = 0, y = 0, width = 1, height = 1) + # Base maps
  draw_grob(legend_image, x = 0.35, y = 0.28, width = 0.22, height = 0.25) # Legend overlay

# Display the plot
final_plot

# Save the final combined plot
ggsave("plots/Figure_3.png", final_plot, width = 12, height = 6)
















# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# --------------------      SPECIES-SPECIFIC SETUP       -----------------------


# SET TIME STEP as "2050's" or "2080's"
time_step <- "2080's"

# SET CLIMATE SCENARIO as "ssp2" or ssp5"
climate_scenario <- "ssp5"

# SET THE SPECIES CODE (using eBird standardized 6-letter codes)
species_code <- "woothr"
species_name <- "Wood Thrush"

# create file paths
species_climate_file_path <- file.path("species_climate_impact_maps", paste0(species_code, ".", time_step, ".", climate_scenario, ".tif"))
species_LULC_file_path <- file.path("species_LULC_impact_maps", paste0(species_code, ".", time_step, ".", climate_scenario, ".tif"))

species_map_climate <- rast(species_climate_file_path)
species_map_LULC <- rast(species_LULC_file_path)


# ---------------    MAPPING SPECIES-SPECIFIC CLIMATE IMPACT   -----------------



raster_values <- values(species_map_climate, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0(species_name, " Climate Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/species_CEI_map.png", width = 4000, height = 4000, res = 1000)

plot(species_map_climate, 
     #main = title, 
     #cex.main = 1.5, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)

dev.off()  # Closes the device


# ---------------    MAPPING SPECIES-SPECIFIC LULC IMPACT   -----------------


raster_values <- values(species_map_LULC, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0(species_name, " LULC Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/species_LEI_map.png", width = 4000, height = 4000, res = 1000)

plot(species_map_LULC, 
     #main = title, 
     #cex.main = 1.5, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)

dev.off()  # Closes the device



# ------------------------------------------------------------------------------

# Load images
CEI_map_image <- rasterGrob(readPNG("plots/species_CEI_map.png"), interpolate = TRUE)
LEI_map_image <- rasterGrob(readPNG("plots/species_LEI_map.png"), interpolate = TRUE)
legend_image <- rasterGrob(readPNG("plots/vulnerability_map_legend.png"), interpolate = TRUE)

# Create labels "A" and "B"
CEI_with_label <- ggdraw() +
  draw_grob(CEI_map_image) +
  draw_label("a", x = 0.1, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate", x = 0.4, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold")


LEI_with_label <- ggdraw() +
  draw_grob(LEI_map_image) +
  draw_label("b", x = 0.1, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Land-Use", x = 0.4, y = 0.95, hjust = 0, vjust = 1, size = 20, fontface = "bold")


# Combine the CEI and LEI maps side by side
side_by_side <- plot_grid(CEI_with_label, LEI_with_label, nrow = 1, rel_widths = c(1, 1))

# Overlay the legend on top of the side-by-side maps
final_plot <- ggdraw() +
  draw_plot(side_by_side, x = 0, y = 0, width = 1, height = 1) + # Base maps
  draw_grob(legend_image, x = 0.35, y = 0.28, width = 0.22, height = 0.25) # Legend overlay

# Display the plot
final_plot

# Save the final combined plot
ggsave("plots/Figure_1.png", final_plot, width = 12, height = 6)
