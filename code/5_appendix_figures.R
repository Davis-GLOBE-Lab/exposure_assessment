# ------------------------------------------------------------------------------

# -------------------    LIST of APPENDIX PLOTS   ------------------------------

# ------------------------------------------------------------------------------

# This script creates Table S1 and Figures 2, S1-2, S16-17, as well as calculates some of the information
# relevant for each plot (i.e. correlation values). See the Github for the final plots

# Last modified by: Matthew Gilbert
# Last modified on: 23 Aug 2026

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
library(stringr)
library(grid)
library(dplyr)
library(gridExtra)
library(maps)
library(dplyr)
library(hexbin)
library(scales) 

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

crop_to_east <- function(raster) {
  east_bad <- st_transform(east, crs(raster))
  raster_east <- mask(crop(raster, vect(east_bad)),east_bad)
  raster <- project(raster_east, crs) |> trim()
  east <- st_transform(east_bad, crs)
  return(raster) }

get_CEI <- function(df) {
  if ("CEI" %in% names(df)) {
    return(df$CEI)
  } else if ("CVI" %in% names(df)) {
    return(df$CVI)
  } else {
    stop("Neither CEI nor CVI column found in data frame.")
  }
}

replace_cei_or_cvi <- function(df, newname) {
  if ("CEI" %in% names(df)) {
    names(df)[names(df) == "CEI"] <- newname
  } else if ("CVI" %in% names(df)) {
    names(df)[names(df) == "CVI"] <- newname
  } else {
    stop(paste("Neither CEI nor CVI found in dataframe for", newname))
  }
  return(df)
}


# ------------------------------------------------------------------------------

#        Figure 2 -- Species-level Correlation between LEI and CEI 

## -----------------------------------------------------------------------------

# choose which model colors to use. Figure 4 uses only one color for SSP2 and one 
# for SSP 5 (the second option below). The first option would plot time steps 
# with unique colors as well
model_colors <- c("2050 SSP2" = "#B6D7E8", "2050 SSP5" = "#FDB774", "2080 SSP2" = "#5989BE", "2080 SSP5" = "#D73027")
model_colors <- c("2050 SSP2" = "#5989BE", "2050 SSP5" = "#D73027", "2080 SSP2" = "#5989BE", "2080 SSP5" = "#D73027")


# Retrieve output files from folder
df_2050_ssp2 <- read.csv("output_data_frames/output_df_2050's_ssp2.csv")
df_2050_ssp2 <- replace_cei_or_cvi(df_2050_ssp2, "CEI")
names(df_2050_ssp2)[names(df_2050_ssp2) == "prop_LULC_change"] <- "LEI"

df_2050_ssp5 <- read.csv("output_data_frames/output_df_2050's_ssp5.csv")
df_2050_ssp5 <- replace_cei_or_cvi(df_2050_ssp5, "CEI")
names(df_2050_ssp5)[names(df_2050_ssp5) == "prop_LULC_change"] <- "LEI"

df_2080_ssp2 <- read.csv("output_data_frames/output_df_2080's_ssp2.csv")
df_2080_ssp2 <- replace_cei_or_cvi(df_2080_ssp2, "CEI")
names(df_2080_ssp2)[names(df_2080_ssp2) == "prop_LULC_change"] <- "LEI"

df_2080_ssp5 <- read.csv("output_data_frames/output_df_2080's_ssp5.csv")
df_2080_ssp5 <- replace_cei_or_cvi(df_2080_ssp5, "CEI")
names(df_2080_ssp5)[names(df_2080_ssp5) == "prop_LULC_change"] <- "LEI"

#
df_2050_ssp2$model <- "2050 SSP2"
df_2050_ssp5$model <- "2050 SSP5"
df_2080_ssp2$model <- "2080 SSP2"
df_2080_ssp5$model <- "2080 SSP5"

full_df <- bind_rows(df_2050_ssp2, df_2050_ssp5, df_2080_ssp2, df_2080_ssp5)
just_2080 <- bind_rows(df_2080_ssp2, df_2080_ssp5)
just_2050 <- bind_rows(df_2050_ssp2, df_2050_ssp5)


# Load plots
plot_2050 <- ggplot(just_2050, aes(x = CEI, y = LEI, color = model, fill = model)) +
  geom_point(alpha = 0.7) +
  #geom_smooth(method = "loess", se = TRUE, fullrange = TRUE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(color = "Model", fill = "Model", x = "CEI", y = "LEI") +  # Keep LEI label for 2050
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.margin = margin(20, 10, 10, 10)) +  # Shift graph down
  coord_cartesian(xlim = c(1, 7), ylim = c(0, 0.6))

plot_2080 <- ggplot(just_2080, aes(x = CEI, y = LEI, color = model, fill = model)) +
  geom_point(alpha = 0.7) +
  #geom_smooth(method = "loess", se = TRUE, fullrange = TRUE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(color = "Model", fill = "Model", x = "CEI", y = NULL) +  # Remove LEI label from 2080
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.margin = margin(20, 10, 10, 10)) +  # Shift graph down
  coord_cartesian(xlim = c(1, 7), ylim = c(0, 0.6))

calculate_cor <- function(data) {
  
  data %>%
    group_by(model) %>%
    group_modify(~{
      test <- cor.test(
        .x$CEI,
        .x$LEI,
        method = "spearman",
        exact = FALSE
      )
      
      tibble(
        n = sum(complete.cases(.x$CEI, .x$LEI)),
        rho = unname(test$estimate),
        statistic = unname(test$statistic),
        p.value = test$p.value
      )
    }) %>%
    ungroup()
}

data <- subset(just_2050, just_2050$model=="2050 SSP2")
cor.test(data$CEI, data$LEI, method="spearman")

# Compute correlation for both 2050 and 2080 datasets
cor_2050 <- calculate_cor(just_2050)
cor_2080 <- calculate_cor(just_2080)

# Combine results into a single list
cor_list <- list("2050" = cor_2050, "2080" = cor_2080)


# Load legend image
glm_legend <- rasterGrob(readPNG("plots/glm_legend_reduced.png"), interpolate = TRUE)

# Label the plots
plot_2050_labeled <- ggdraw() +
  draw_plot(plot_2050, y = -0.02) +  # Shift down slightly
  draw_label("a", x = 0.05, y = 0.99, hjust = 0, vjust = 1, size = 18, fontface = "bold") +
  draw_label("2050", x = 0.45, y = 1, hjust = 0, vjust = 1, size = 18, fontface = "bold")

plot_2080_labeled <- ggdraw() +
  draw_plot(plot_2080, y = -0.02) +  # Shift down slightly
  draw_label("b", x = 0.05, y = 0.99, hjust = 0, vjust = 1, size = 18, fontface = "bold") +
  draw_label("2080", x = 0.45, y = 1, hjust = 0, vjust = 1, size = 18, fontface = "bold")

# Arrange plots side by side
correlation_grid <- plot_grid(
  plot_2050_labeled, plot_2080_labeled,
  ncol = 2
)

# Combine with legend
final_plot <- plot_grid(
  correlation_grid,
  ggdraw() + draw_grob(glm_legend),
  ncol = 1,
  rel_heights = c(1, 0.1))  # Adjust legend size

# Save final plot
ggsave("plots/Figure_2.png", final_plot, width = 9, height = 5)



# ------------------------------------------------------------------------------

#  Sup. Data S1 -- Full Data table of CEI and LEI for each species in each scenario

## -----------------------------------------------------------------------------

# read in data files for each time step/emissions scenarios and rename columns
df_2050_ssp2 <- read.csv("output_data_frames/output_df_2050's_ssp2.csv")
df_2050_ssp2 <- replace_cei_or_cvi(df_2050_ssp2, "CEI_ssp2_2050")
names(df_2050_ssp2)[names(df_2050_ssp2) == "CEI"] <- "CEI_ssp2_2050"
names(df_2050_ssp2)[names(df_2050_ssp2) == "prop_LULC_change"] <- "LEI_ssp2_2050"

df_2050_ssp5 <- read.csv("output_data_frames/output_df_2050's_ssp5.csv")
df_2050_ssp5 <- replace_cei_or_cvi(df_2050_ssp5, "CEI_ssp5_2050")
names(df_2050_ssp5)[names(df_2050_ssp5) == "CEI"] <- "CEI_ssp5_2050"
names(df_2050_ssp5)[names(df_2050_ssp5) == "prop_LULC_change"] <- "LEI_ssp5_2050"

df_2080_ssp2 <- read.csv("output_data_frames/output_df_2080's_ssp2.csv")
df_2080_ssp2 <- replace_cei_or_cvi(df_2080_ssp2, "CEI_ssp2_2080")
names(df_2080_ssp2)[names(df_2080_ssp2) == "CEI"] <- "CEI_ssp2_2080"
names(df_2080_ssp2)[names(df_2080_ssp2) == "prop_LULC_change"] <- "LEI_ssp2_2080"

df_2080_ssp5 <- read.csv("output_data_frames/output_df_2080's_ssp5.csv")
df_2080_ssp5 <- replace_cei_or_cvi(df_2080_ssp5, "CEI_ssp5_2080")
names(df_2080_ssp5)[names(df_2080_ssp5) == "CEI"] <- "CEI_ssp5_2080"
names(df_2080_ssp5)[names(df_2080_ssp5) == "prop_LULC_change"] <- "LEI_ssp5_2080"


# add in 2050 ssp2 and 2050 ssp5
combined <- df_2050_ssp2 %>% left_join(df_2050_ssp5, by="species_name")
combined <- combined %>% 
  rename_at(vars(ends_with(".x")), ~str_replace(., "\\..$","")) %>% select_at(vars(-ends_with(".y")))

# add in 2080 ssp2
combined <- combined %>% left_join(df_2080_ssp2, by="species_name")
combined <- combined %>% 
  rename_at(vars(ends_with(".x")), ~str_replace(., "\\..$","")) %>% select_at(vars(-ends_with(".y")))

# addin 2080 ssp5
combined <- combined %>% left_join(df_2080_ssp5, by="species_name")
combined <- combined %>% 
  rename_at(vars(ends_with(".x")), ~str_replace(., "\\..$","")) %>% select_at(vars(-ends_with(".y")))

#reorder column names
combined <- combined %>% 
  dplyr::select(X, species_name, species_code, family, habitat, guild, year, proportion_in_east, mean_longitude, mean_latitude, CEI_ssp2_2050, LEI_ssp2_2050, CEI_ssp5_2050, LEI_ssp5_2050, CEI_ssp2_2080, LEI_ssp2_2080, CEI_ssp5_2080, LEI_ssp5_2080)

# round values to 3 decimel points
rounded <- combined %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# write out CSV
write.csv(rounded, "plots/Table_S1.csv", row.names = FALSE)





# ------------------------------------------------------------------------------

#    Figure S1 -- Maps of Climate Velocity Projections in Different Scenarios

## -----------------------------------------------------------------------------


#load in file for ssp2 2050
time_step <- "2041_2070"
climate_scenario <- "245"
climate_file_path <- file.path(paste0("climate_velocity_data/fwvel731_ensemble_8gcm_", climate_scenario, "_", time_step, ".tif"))
map_climate_ssp2_2050 <- rast(climate_file_path)

#load in file for ssp5 2050
time_step <- "2041_2070"
climate_scenario <- "585"
climate_file_path <- file.path(paste0("climate_velocity_data/fwvel731_ensemble_8gcm_", climate_scenario, "_", time_step, ".tif"))
map_climate_ssp5_2050 <- rast(climate_file_path)

#load in file for ssp2 2080
time_step <- "2071_2100"
climate_scenario <- "245"
climate_file_path <- file.path(paste0("climate_velocity_data/fwvel731_ensemble_8gcm_", climate_scenario, "_", time_step, ".tif"))
map_climate_ssp2_2080 <- rast(climate_file_path)

#load in file for ssp5 2080
time_step <- "2071_2100"
climate_scenario <- "585"
climate_file_path <- file.path(paste0("climate_velocity_data/fwvel731_ensemble_8gcm_", climate_scenario, "_", time_step, ".tif"))
map_climate_ssp5_2080 <- rast(climate_file_path)

#crop files to Eastern US
map_climate_ssp2_2050 <- crop_to_east(map_climate_ssp2_2050)
map_climate_ssp5_2050 <- crop_to_east(map_climate_ssp5_2050)
map_climate_ssp2_2080 <- crop_to_east(map_climate_ssp2_2080)
map_climate_ssp5_2080 <- crop_to_east(map_climate_ssp5_2080)

#create accurate percentile bins for plots
raster_values_ssp2_2050 <- values(map_climate_ssp2_2050, na.rm = TRUE)
raster_values_ssp5_2050 <- values(map_climate_ssp5_2050, na.rm = TRUE)
raster_values_ssp2_2080 <- values(map_climate_ssp2_2080, na.rm = TRUE)
raster_values_ssp5_2080 <- values(map_climate_ssp5_2080, na.rm = TRUE)

# calculate means of rasters
mean(raster_values_ssp2_2050)
mean(raster_values_ssp5_2050)
mean(raster_values_ssp2_2080)
mean(raster_values_ssp5_2080)

# binning values ACROSS all time steps / scenarios to show differences in climate velocity between scenarios
all_values <- c(raster_values_ssp2_2050, raster_values_ssp5_2050, raster_values_ssp2_2080, raster_values_ssp5_2080)
brks <- unique(quantile(all_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(200))


## ---------------------- creating a legend

# Set up a blank plot and fill it with the gradient
png("gradient_palette.png", width = 800, height = 200)
par(mar = c(1, 1, 1, 1)) # Remove margins

# Create a gradient rectangle
plot(1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", xlim = c(0, 200), ylim = c(0, 1))
rect(xleft = 0:(length(pal) - 1), xright = 1:length(pal), ybottom = 0, ytop = 1, col = pal, border = NA)

# Close the graphics device
dev.off()


## ---------------------- saving each plot as an image

png("plots/climate_ssp2_2050.png", width = 3500, height = 4000, res = 1000)
plot(map_climate_ssp2_2050, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()

png("plots/climate_ssp5_2050.png", width = 3500, height = 4000, res = 1000)
plot(map_climate_ssp5_2050, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()

png("plots/climate_ssp2_2080.png", width = 3500, height = 4000, res = 1000)
plot(map_climate_ssp2_2080, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()

png("plots/climate_ssp5_2080.png", width = 3500, height = 4000, res = 1000)
plot(map_climate_ssp5_2080, 
     axes = F,
     breaks = brks, 
     col = pal, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()


## ---------------------- pulling in images and collating into a single plot

ssp2_2050_image <- rasterGrob(readPNG("plots/climate_ssp2_2050.png"), interpolate = TRUE)
ssp5_2050_image <- rasterGrob(readPNG("plots/climate_ssp5_2050.png"), interpolate = TRUE)
ssp2_2080_image <- rasterGrob(readPNG("plots/climate_ssp2_2080.png"), interpolate = TRUE)
ssp5_2080_image <- rasterGrob(readPNG("plots/climate_ssp5_2080.png"), interpolate = TRUE)
climate_velocity_legend <- rasterGrob(readPNG("plots/climate_velocity_legend.png"), interpolate = TRUE)

ssp2_2050_image <- ggdraw() +
  draw_grob(ssp2_2050_image) +
  draw_label("A", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("SSP2 2050", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")


ssp5_2050_image <- ggdraw() +
  draw_grob(ssp5_2050_image) +
  draw_label("B", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("SSP5 2050", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")


ssp2_2080_image <- ggdraw() +
  draw_grob(ssp2_2080_image) +
  draw_label("C", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("SSP2 2080", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")

ssp5_2080_image <- ggdraw() +
  draw_grob(ssp5_2080_image) +
  draw_label("D", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("SSP5 2080", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")


# arrange in a grid
plots_grid <- plot_grid(
  ssp2_2050_image, ssp5_2050_image,
  ssp2_2080_image, ssp5_2080_image,
  ncol = 2
)

# combine legend and grid
final_plot <- plot_grid(
  plots_grid, 
  ggdraw() + draw_grob(climate_velocity_legend), 
  ncol = 1, 
  rel_heights = c(1, 0.09) # Adjust relative heights of the grid and legend
)

ggsave("plots/Figure_S1.png", final_plot, width = 10, height = 11)








# ------------------------------------------------------------------------------

#        Figure S2 -- Maps of LULC  Projections in Different Scenarios

## -----------------------------------------------------------------------------

#load in file for ssp2 2050
time_step <- "2055"
LULC_scenario <- "SSP2_RCP45"
LULC_file_path <- file.path(paste0("LULC_change_data/global_", LULC_scenario, "_", time_step, ".tif"))
map_LULC_ssp2_2050 <- rast(LULC_file_path)

#load in file for ssp5 2050
time_step <- "2055"
LULC_scenario <- "SSP5_RCP85"
LULC_file_path <- file.path(paste0("LULC_change_data/global_", LULC_scenario, "_", time_step, ".tif"))
map_LULC_ssp5_2050 <- rast(LULC_file_path)

#load in file for ssp2 2080
time_step <- "2085"
LULC_scenario <- "SSP2_RCP45"
LULC_file_path <- file.path(paste0("LULC_change_data/global_", LULC_scenario, "_", time_step, ".tif"))
map_LULC_ssp2_2080 <- rast(LULC_file_path)

#load in file for ssp5 2080
time_step <- "2085"
LULC_scenario <- "SSP5_RCP85"
LULC_file_path <- file.path(paste0("LULC_change_data/global_", LULC_scenario, "_", time_step, ".tif"))
map_LULC_ssp5_2080 <- rast(LULC_file_path)

#crop files to Eastern US
map_LULC_ssp2_2050 <- crop_to_east(map_LULC_ssp2_2050)
map_LULC_ssp5_2050 <- crop_to_east(map_LULC_ssp5_2050)
map_LULC_ssp2_2080 <- crop_to_east(map_LULC_ssp2_2080)
map_LULC_ssp5_2080 <- crop_to_east(map_LULC_ssp5_2080)


# ------------- set colors for different LULC categories

category_colors <- c(
  "darkblue",  # Blue for water
  #"#66D1B0",  # Forest green for forest
  #"#70db65", 
  "#50cf61", 
  "#CC79A7",  # Lime-yellow for grassland
  "#E69F00",  # Orange for barren
  "#F0E442",  # Pinkish-purple for cropland F0E442
  "#D41159"   # Bright red for urban (stands out)
)

# ------------- Save each subplot as an image


png("plots/LULC_ssp2_2050.png", width = 3500, height = 4000, res = 1000)
plot(map_LULC_ssp2_2050, 
     axes = F,
     col = category_colors, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()

png("plots/LULC_ssp5_2050.png", width = 3500, height = 4000, res = 1000)
plot(map_LULC_ssp5_2050, 
     axes = F,
     col = category_colors, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()

png("plots/LULC_ssp2_2080.png", width = 3500, height = 4000, res = 1000)
plot(map_LULC_ssp2_2080, 
     axes = F,
     col = category_colors, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()

png("plots/LULC_ssp5_2080.png", width = 3500, height = 4000, res = 1000)
plot(map_LULC_ssp5_2080, 
     axes = F,
     col = category_colors, 
     legend = F, 
     box = F)
plot(east, border = "black", col = NA, add = TRUE)
dev.off()



## ---------------------- pull in subplots and collate into a single plot

ssp2_2050_image <- rasterGrob(readPNG("plots/LULC_ssp2_2050.png"), interpolate = TRUE)
ssp5_2050_image <- rasterGrob(readPNG("plots/LULC_ssp5_2050.png"), interpolate = TRUE)
ssp2_2080_image <- rasterGrob(readPNG("plots/LULC_ssp2_2080.png"), interpolate = TRUE)
ssp5_2080_image <- rasterGrob(readPNG("plots/LULC_ssp5_2080.png"), interpolate = TRUE)
LULC_velocity_legend <- rasterGrob(readPNG("plots/LULC_legend.png"), interpolate = TRUE)

ssp2_2050_image <- ggdraw() +
  draw_grob(ssp2_2050_image) +
  draw_label("A", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("SSP2 2050", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")

ssp5_2050_image <- ggdraw() +
  draw_grob(ssp5_2050_image) +
  draw_label("B", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") +
  draw_label("SSP5 2050", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")

ssp2_2080_image <- ggdraw() +
  draw_grob(ssp2_2080_image) +
  draw_label("C", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("SSP2 2080", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")

ssp5_2080_image <- ggdraw() +
  draw_grob(ssp5_2080_image) +
  draw_label("D", x = 0.12, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("SSP5 2080", x = 0.38, y = 0.93, hjust = 0, vjust = 1, size = 20, fontface = "bold")

# arrange in a grid
plots_grid <- plot_grid(
  ssp2_2050_image, ssp5_2050_image,
  ssp2_2080_image, ssp5_2080_image,
  ncol = 2
)

# combine legend and grid
final_plot <- plot_grid(
  plots_grid, 
  ggdraw() + draw_grob(LULC_velocity_legend), 
  ncol = 1, 
  rel_heights = c(1, 0.07)) # Adjust relative heights of the grid and legend)

ggsave("plots/Figure_S2.png", final_plot, width = 10, height = 11)





# ------------------------------------------------------------------------------

#    Figure S11 -- Spatial Correlation between Community-level Exposure Maps

## -----------------------------------------------------------------------------


# load in maps for selected time step/scenario
time_step <- "2080's"
scenario <- "ssp5"
climate_file_path <- file.path(paste0("community_maps/community_map_climate_", time_step, "_", scenario, ".tif"))
LULC_file_path <- file.path(paste0("community_maps/community_map_LULC_", time_step, "_", scenario, ".tif"))
community_map_climate <- rast(climate_file_path)
community_map_LULC <- rast(LULC_file_path)

# Ensure both rasters have the same extent, resolution, and projection
# If they are not aligned, resample or crop as needed
if (!compareGeom(community_map_climate, community_map_LULC, stopOnError = FALSE)) {
  community_map_LULC <- resample(community_map_LULC, community_map_climate)
}

# Extract cell values as vectors
climate_values <- values(community_map_climate)
LULC_values <- values(community_map_LULC)

# Remove NA values from both datasets to ensure they align properly
data <- data.frame(
  climate = climate_values,
  LULC = LULC_values
)
data <- na.omit(data)
colnames(data) <- c("climate", "LULC")

cor.test(data$climate, data$LULC, method = "pearson")



# Convert climate and LULC values to percentiles
n <- nrow(data)
data$climate <- rank(data$climate, ties.method = "average") / n
data$LULC    <- rank(data$LULC,    ties.method = "average") / n

cor.test(data$climate, data$LULC)

# Calculate the hexbin density
hex_data <- ggplot(data, aes(x = climate, y = LULC)) +
  stat_bin_hex(bins = 40) +
  scale_fill_viridis_c()

# Extract hexbin counts and convert them to density percentiles
hex_data <- ggplot_build(hex_data)$data[[1]] %>%
  mutate(density_percentile = ecdf(count)(count))  # Calculate the percentile of each hex count

pal <- colorRampPalette(brewer.pal(7, "Blues"))(200)

# Plot the hexbin with density percentiles and add geom_smooth
plot <- ggplot() +
  geom_hex(data = hex_data, aes(x = x, y = y, fill = density_percentile), stat = "identity") +
  scale_fill_gradientn(colors = pal, name = "Density Percentile\n", labels = percent_format()) +
  labs(x = "Community Climate Exposure Map Percentiles",
       y = "Community LULC Exposure Map Percentiles") +
  theme_minimal(17) + 
  theme(axis.text = element_text(size= 15), 
        legend.text = element_text(size = 14))

plot
ggsave("plots/Figure_S11.png", plot, width = 8, height = 5, dpi = 500)



# ------------------------------------------------------------------------------

#     Figure S16-- Species-level Exposure Maps in Each Scenario for WOOTHR

## -----------------------------------------------------------------------------


# Define Scenarios
scenarios <- list(
  list(time_step = "2050's", climate_scenario = "ssp2"),
  list(time_step = "2050's", climate_scenario = "ssp5"),
  list(time_step = "2080's", climate_scenario = "ssp2")
)

# Define species
species_code <- "woothr" # Wood Thrush
species_name <- "Wood Thrush"

# ------------------------------------------------------------------------------
# Loop through scenarios and generate maps
for (scenario in scenarios) {
  
  time_step <- scenario$time_step
  climate_scenario <- scenario$climate_scenario
  
  # Create file paths
  species_climate_file_path <- file.path(
    "species_climate_impact_maps", paste0(species_code, ".", time_step, ".", climate_scenario, ".tif")
  )
  species_LULC_file_path <- file.path(
    "species_LULC_impact_maps", paste0(species_code, ".", time_step, ".", climate_scenario, ".tif")
  )
  
  species_map_climate <- rast(species_climate_file_path)
  species_map_LULC <- rast(species_LULC_file_path)
  
  # ---------------------------------------------------------------------------
  #                MAPPING SPECIES-SPECIFIC CLIMATE IMPACT  
  
  
  raster_values <- values(species_map_climate, na.rm = TRUE)
  brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
  pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
  
  png(paste0("plots/species_CEI_map_", climate_scenario, "_", time_step, ".png"), width = 4000, height = 4000, res = 1000)
  
  plot(species_map_climate, 
       axes = F,
       breaks = brks, 
       col = pal, 
       legend = F, 
       box = F)
  plot(east, border = "black", col = NA, add = TRUE)
  
  dev.off()  # Closes the device
  
  # ---------------------------------------------------------------------------
  #                MAPPING SPECIES-SPECIFIC LULC IMPACT  
  
  
  raster_values <- values(species_map_LULC, na.rm = TRUE)
  brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
  
  png(paste0("plots/species_LEI_map_", climate_scenario, "_", time_step, ".png"), width = 4000, height = 4000, res = 1000)
  
  plot(species_map_LULC, 
       axes = F,
       breaks = brks, 
       col = pal, 
       legend = F, 
       box = F)
  plot(east, border = "black", col = NA, add = TRUE)
  
  dev.off()  # Closes the device
}

# ------------------------------------------------------------------------------
# Load images for visualization
CEI_map_image_ssp2_2050 <- rasterGrob(readPNG("plots/species_CEI_map_ssp2_2050's.png"), interpolate = TRUE)
LEI_map_image_ssp2_2050 <- rasterGrob(readPNG("plots/species_LEI_map_ssp2_2050's.png"), interpolate = TRUE)

CEI_map_image_ssp5_2050 <- rasterGrob(readPNG("plots/species_CEI_map_ssp5_2050's.png"), interpolate = TRUE)
LEI_map_image_ssp5_2050 <- rasterGrob(readPNG("plots/species_LEI_map_ssp5_2050's.png"), interpolate = TRUE)

CEI_map_image_ssp2_2080 <- rasterGrob(readPNG("plots/species_CEI_map_ssp2_2080's.png"), interpolate = TRUE)
LEI_map_image_ssp2_2080 <- rasterGrob(readPNG("plots/species_LEI_map_ssp2_2080's.png"), interpolate = TRUE)

legend_image <- rasterGrob(readPNG("plots/vulnerability_map_legend.png"), interpolate = TRUE)

# ------------------------------------------------------------------------------
# Create labeled plots
label_positions <- list(x = 0.1, y = 0.95)
ssp_positions <- list(x = 0.3, y = 0.95)

CEI_ssp2_2050 <- ggdraw() + draw_grob(CEI_map_image_ssp2_2050) +
  draw_label("A", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate SSP2 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")

LEI_ssp2_2050 <- ggdraw() + draw_grob(LEI_map_image_ssp2_2050) +
  draw_label("B", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("LULC SSP2 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")

CEI_ssp5_2050 <- ggdraw() + draw_grob(CEI_map_image_ssp5_2050) +
  draw_label("C", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate SSP5 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")

LEI_ssp5_2050 <- ggdraw() + draw_grob(LEI_map_image_ssp5_2050) +
  draw_label("D", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("LULC SSP5 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")

CEI_ssp2_2080 <- ggdraw() + draw_grob(CEI_map_image_ssp2_2080) +
  draw_label("E", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate SSP2 2080", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")

LEI_ssp2_2080 <- ggdraw() + draw_grob(LEI_map_image_ssp2_2080) +
  draw_label("F", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("LULC SSP2 2080", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")

# ------------------------------------------------------------------------------
# Arrange plots in 3 rows
final_plot <- plot_grid(
  CEI_ssp2_2050, LEI_ssp2_2050,
  CEI_ssp5_2050, LEI_ssp5_2050,
  CEI_ssp2_2080, LEI_ssp2_2080,
  nrow = 3, align = "v"
)

# Overlay the legend
final_plot <- ggdraw() +
  draw_plot(final_plot, x = 0, y = 0, width = 1, height = 1) + # Base maps
  draw_grob(legend_image, x = 0.38, y = 0.00, width = 0.22, height = 0.25) # Adjusted legend placement

# ------------------------------------------------------------------------------
# Display and save the final plot
final_plot
ggsave("plots/Figure_S16.png", final_plot, width = 12, height = 12)




# ------------------------------------------------------------------------------

#        Figure S_17 -- Community-level Exposure Maps in Each Scenario

## -----------------------------------------------------------------------------


# Scenario
time_step <- "2050's"
climate_scenario <- "ssp2"

# pull files
climate_file_path <- file.path("community_maps", paste0("community_map_climate_", time_step, "_", climate_scenario, ".tif"))
LULC_file_path <- file.path("community_maps", paste0("community_map_LULC_", time_step, "_", climate_scenario, ".tif"))
community_map_climate <- rast(climate_file_path)
community_map_LULC <- rast(LULC_file_path)

raster_values <- values(community_map_climate, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community Climate Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_CEI_map_ssp2_2050.png", width = 4000, height = 4000, res = 1000)

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


raster_values <- values(community_map_LULC, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community LULC Change Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_LEI_map_ssp2_2050.png", width = 4000, height = 4000, res = 1000)

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


#####

# Scenario
time_step <- "2050's"
climate_scenario <- "ssp5"

# pull files
climate_file_path <- file.path("community_maps", paste0("community_map_climate_", time_step, "_", climate_scenario, ".tif"))
LULC_file_path <- file.path("community_maps", paste0("community_map_LULC_", time_step, "_", climate_scenario, ".tif"))
community_map_climate <- rast(climate_file_path)
community_map_LULC <- rast(LULC_file_path)

raster_values <- values(community_map_climate, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community Climate Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_CEI_map_ssp5_2050.png", width = 4000, height = 4000, res = 1000)

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


raster_values <- values(community_map_LULC, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community LULC Change Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_LEI_map_ssp5_2050.png", width = 4000, height = 4000, res = 1000)

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



#####

# Scenario
time_step <- "2080's"
climate_scenario <- "ssp2"

# pull files
climate_file_path <- file.path("community_maps", paste0("community_map_climate_", time_step, "_", climate_scenario, ".tif"))
LULC_file_path <- file.path("community_maps", paste0("community_map_LULC_", time_step, "_", climate_scenario, ".tif"))
community_map_climate <- rast(climate_file_path)
community_map_LULC <- rast(LULC_file_path)

raster_values <- values(community_map_climate, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community Climate Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_CEI_map_ssp2_2080.png", width = 4000, height = 4000, res = 1000)

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


raster_values <- values(community_map_LULC, na.rm = TRUE)
brks <- unique(quantile(raster_values, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
pal <- colorRampPalette(brewer.pal(7, "GnBu"))(200)
title <- paste0("Community LULC Change Impact Map (", climate_scenario, ", ", time_step, ")")

png("plots/community_LEI_map_ssp2_2080.png", width = 4000, height = 4000, res = 1000)

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
CEI_map_image_ssp2_2050 <- rasterGrob(readPNG("plots/community_CEI_map_ssp2_2050.png"), interpolate = TRUE)
LEI_map_image_ssp2_2050 <- rasterGrob(readPNG("plots/community_LEI_map_ssp2_2050.png"), interpolate = TRUE)

CEI_map_image_ssp5_2050 <- rasterGrob(readPNG("plots/community_CEI_map_ssp5_2050.png"), interpolate = TRUE)
LEI_map_image_ssp5_2050 <- rasterGrob(readPNG("plots/community_LEI_map_ssp5_2050.png"), interpolate = TRUE)

CEI_map_image_ssp2_2080 <- rasterGrob(readPNG("plots/community_CEI_map_ssp2_2080.png"), interpolate = TRUE)
LEI_map_image_ssp2_2080 <- rasterGrob(readPNG("plots/community_LEI_map_ssp2_2080.png"), interpolate = TRUE)

legend_image <- rasterGrob(readPNG("plots/vulnerability_map_legend.png"), interpolate = TRUE)

# ------------------------------------------------------------------------------

# Create labeled plots
label_positions <- list(x = 0.1, y = 0.95)
ssp_positions <- list(x = 0.3, y = 0.95)

CEI_ssp2_2050 <- ggdraw() + draw_grob(CEI_map_image_ssp2_2050) +
  draw_label("A", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate SSP2 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")


LEI_ssp2_2050 <- ggdraw() + draw_grob(LEI_map_image_ssp2_2050) +
  draw_label("B", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("LULC SSP2 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")


CEI_ssp5_2050 <- ggdraw() + draw_grob(CEI_map_image_ssp5_2050) +
  draw_label("C", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate SSP5 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")


LEI_ssp5_2050 <- ggdraw() + draw_grob(LEI_map_image_ssp5_2050) +
  draw_label("D", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("LULC SSP5 2050", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")


CEI_ssp2_2080 <- ggdraw() + draw_grob(CEI_map_image_ssp2_2080) +
  draw_label("E", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("Climate SSP2 2080", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")


LEI_ssp2_2080 <- ggdraw() + draw_grob(LEI_map_image_ssp2_2080) +
  draw_label("F", x = label_positions$x, y = label_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold") + 
  draw_label("LULC SSP2 2080", x = ssp_positions$x, y = ssp_positions$y, hjust = 0, vjust = 1, size = 20, fontface = "bold")


# ------------------------------------------------------------------------------

# Arrange plots in 3 rows
final_plot <- plot_grid(
  CEI_ssp2_2050, LEI_ssp2_2050,
  CEI_ssp5_2050, LEI_ssp5_2050,
  CEI_ssp2_2080, LEI_ssp2_2080,
  nrow = 3, align = "v"
)

# Overlay the legend
final_plot <- ggdraw() +
  draw_plot(final_plot, x = 0, y = 0, width = 1, height = 1) + # Base maps
  draw_grob(legend_image, x = 0.38, y = 0.00, width = 0.22, height = 0.25) # Adjusted legend placement

# ------------------------------------------------------------------------------

# Display the final plot
final_plot

# Save the final combined plot
ggsave("plots/Figure_S13.png", final_plot, width = 12, height = 12)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------