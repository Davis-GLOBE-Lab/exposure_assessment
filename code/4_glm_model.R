# ------------------------------------------------------------------------------

# -----------    GROUP-LEVEL ANALYSIS with MODELING    -------------------------

# ------------------------------------------------------------------------------

# Uses data frames calculated by the community_raster_builder.r script to analyze
# differences in climate/LULC change exposure between groups of species using
# a linear model

# Last modified by: Matthew Gilbert
# Last modified on: 5 Aug 2026

# set working directory
setwd("~/Desktop/vulnerability_analysis")

# load packages
library(ggplot2)
library(terra)
library(car)
library(broom)
library(betareg)
library(patchwork)
library(cowplot)
library(MASS)
library(grid)
library(gridExtra)
library(png)
library(dplyr)
library(DHARMa)
library(glmmTMB)
library(gap)
library(gap.datasets)


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

# Retrieve output files from folder
df_2050_ssp2 <- read.csv("output_data_frames/output_df_2050's_ssp2.csv")
df_2050_ssp2 <- replace_cei_or_cvi(df_2050_ssp2, "CEI")

df_2050_ssp5 <- read.csv("output_data_frames/output_df_2050's_ssp5.csv")
df_2050_ssp5 <- replace_cei_or_cvi(df_2050_ssp5, "CEI")

df_2080_ssp2 <- read.csv("output_data_frames/output_df_2080's_ssp2.csv")
df_2080_ssp2 <- replace_cei_or_cvi(df_2080_ssp2, "CEI")

df_2080_ssp5 <- read.csv("output_data_frames/output_df_2080's_ssp5.csv")
df_2080_ssp5 <- replace_cei_or_cvi(df_2080_ssp5, "CEI")




# -------------------------   SETUP DATA    ------------------------------------


# Define the function to process data for each data frame

process_data <- function(df) {
  
  # Define Passerine families
  passerine_families <- c(
    "Alaudidae", "Bombycillidae", "Calcariidae", "Cardinalidae", "Certhiidae", 
    "Corvidae", "Fringillidae", "Hirundinidae", "Icteridae", "Icteriidae", 
    "Laniidae", "Mimidae", "Motacillidae", "Paridae", "Parulidae", 
    "Passerellidae", "Passeridae", "Polioptilidae", "Regulidae", "Sittidae", 
    "Sturnidae", "Troglodytidae", "Turdidae", "Tyrannidae", "Vireonidae")
  
  # Create a new column "order" to classify as "Passerine" or "Non-Passerine"
  df$order <- ifelse(df$family %in% passerine_families, "Passerine", "Non-Passerine")
  
  # Group all vertebrate hunters together
  df <- df %>%
    mutate(guild = case_when(
      guild %in% c("Mammal Hunter", "Fish Hunter", "Bird Hunter") ~ "Vertebrate Hunter",
      TRUE ~ guild))
  
  # Count the number of factors in each habitat and filter for habitats with at least 5 occurrences
  habitat_counts <- table(df$habitat)
  valid_habitats <- names(habitat_counts[habitat_counts >= 5])
  
  # call Prop_LULC_change LEI (Land Exposure Index) for simpler coding
  df$LEI <- df$prop_LULC_change
  
  # Count the number of factors in each guild and filter for guilds with at least 5 occurrences
  guild_counts <- table(df$guild)
  valid_guilds <- names(guild_counts[guild_counts >= 5])
  
  # Subset the data to include only valid habitats and guilds
  df_filtered <- df %>%
    filter(habitat %in% valid_habitats, guild %in% valid_guilds)
  
  # Ensure categorical columns are factors
  df_filtered$order <- as.factor(df_filtered$order)
  df_filtered$habitat <- as.factor(df_filtered$habitat)
  df_filtered$guild <- as.factor(df_filtered$guild)
  
  # Return the processed data frame
  return(df_filtered)
}

# Apply the function to each of your data frames
df_2050_ssp2 <- process_data(df_2050_ssp2)
df_2050_ssp5 <- process_data(df_2050_ssp5)
df_2080_ssp2 <- process_data(df_2080_ssp2)
df_2080_ssp5 <- process_data(df_2080_ssp5)



# --------------------- EVALUATE FOR MULTICOLLINEARITY  ----------------------------
calculate_vif <- function(df, scenario, date) {
  
  lat_c  <- scale(df$mean_latitude, center = TRUE, scale = FALSE)
  long_c <- scale(df$mean_longitude, center = TRUE, scale = FALSE)
  
  # Since response variable is not included in VIF analysis, CEI or LEI here 
  # produce equivalent results
  global_model <- lm(
    CEI ~ guild + habitat + order + 
      lat_c + I(lat_c^2) +
      long_c + I(long_c^2),
    data = df
  ) 
  
  vif_table <- as.data.frame(car::vif(global_model))
  
  vif_table$Variable <- rownames(vif_table)
  rownames(vif_table) <- NULL
  
  # Clean variable names
  vif_table$Variable <- dplyr::recode(
    vif_table$Variable,
    "I(lat_c^2)" = "lat_c²",
    "I(long_c^2)" = "long_c²"
  )
  
  # For factors, keep adjusted GVIF.
  # For 1-df terms, report ordinary VIF.
  vif_table <- vif_table %>%
    mutate(
      VIF = GVIF, 
      Adjusted_GVIF = ifelse(Df > 1, GVIF^(1/(2 * Df)), NA),
      Scenario = scenario,
      Date = date
    ) %>%
    dplyr::select(
      Scenario,
      Date,
      Variable,
      Df,
      VIF,
      Adjusted_GVIF
    )
  
  return(vif_table)
}

# Picking any scenario since predictors are the same for all models
vif_results <- calculate_vif(df_2050_ssp2, "SSP2", "2050")

print(vif_results)

write.csv(
  vif_results,
  "plots/Table_S2.csv",
  row.names = FALSE
)




# -------------------------   CEI FOR LATITUDE MODEL  ----------------------------------


create_CEI_model_and_predict <- function(df) {
  
  library(MuMIn)
  library(dplyr)
  
  # Required by MuMIn
  old_na <- getOption("na.action")
  options(na.action = "na.fail")
  on.exit(options(na.action = old_na))
  
  # -----------------------------
  # Global model
  # -----------------------------
  global_model <- lm(
    CEI ~ guild + habitat + order +
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    data = df
  )
  
  # -----------------------------
  # Model selection
  # -----------------------------
  dredge_table <- dredge(
    global_model,
    rank = "AIC",
    subset =
      dc(mean_latitude, I(mean_latitude^2)) &
      dc(mean_longitude, I(mean_longitude^2))
  )
  
  # Convert to data frame for saving / inspection
  model_comparison <- as.data.frame(dredge_table)
  
  model_comparison$R2 <- sapply(
    get.models(dredge_table, subset = TRUE),
    function(m) summary(m)$r.squared
  )
  
  model_comparison <- model_comparison %>%
    mutate(
      Structure = paste0(
        ifelse(!is.na(guild),   "guild+",   ""),
        ifelse(!is.na(order),   "order+",   ""),
        ifelse(!is.na(habitat), "habitat+", ""),
        ifelse(
          !is.na(mean_latitude) | !is.na(`I(mean_latitude^2)`),
          "lat+", ""
        ),
        ifelse(
          !is.na(mean_longitude) | !is.na(`I(mean_longitude^2)`),
          "long+", ""
        )
      ),
      Structure = sub("\\+$", "", Structure)
    ) %>%
    transmute(
      Structure,
      k = df,
      logLik = logLik,
      delta_AIC = delta,
      weight = weight, 
      R2 = round(R2, 2) 
      )
  
  # -----------------------------
  # Best model
  # -----------------------------
  CEI_model <- get.models(dredge_table, 1)[[1]]

  # -----------------------------
  # DHARMa diagnostics
  # -----------------------------
  simulationOutput <- simulateResiduals(
    fittedModel = CEI_model,
    plot = FALSE
  )
  
  plot(simulationOutput)
  
  # -----------------------------
  # Prediction grid
  # -----------------------------
  latitude_range <- seq(
    min(df$mean_latitude, na.rm = TRUE),
    max(df$mean_latitude, na.rm = TRUE),
    length.out = 100
  )
  
  prediction_data <- data.frame(
    guild = factor("Foliage Gleaner", levels = levels(df$guild)),
    habitat = factor("Forests", levels = levels(df$habitat)),
    order = factor("Passerine", levels = levels(df$order)),
    mean_latitude = latitude_range,
    mean_longitude = mean(df$mean_longitude, na.rm = TRUE)
  )
  
  preds <- predict(
    CEI_model,
    newdata = prediction_data,
    se.fit = TRUE
  )
  
  prediction_data <- prediction_data %>%
    mutate(
      fit   = preds$fit,
      lower = preds$fit - 1.96 * preds$se.fit,
      upper = preds$fit + 1.96 * preds$se.fit
    )
  
  # -----------------------------
  # Return
  # -----------------------------
  return(list(
    model = CEI_model,
    predictions = prediction_data,
    model_selection = model_comparison
  ))
}


# Apply the function to each data frame and store results
pred_2050_ssp2 <- create_CEI_model_and_predict(df_2050_ssp2)$predictions
pred_2050_ssp5 <- create_CEI_model_and_predict(df_2050_ssp5)$predictions
pred_2080_ssp2 <- create_CEI_model_and_predict(df_2080_ssp2)$predictions
pred_2080_ssp5 <- create_CEI_model_and_predict(df_2080_ssp5)$predictions


# Combine predictions into one data frame with a model identifier
pred_2050_ssp2$model <- "2050 SSP2"
pred_2050_ssp5$model <- "2050 SSP5"
pred_2080_ssp2$model <- "2080 SSP2"
pred_2080_ssp5$model <- "2080 SSP5"

all_predictions_CEI_lat <- bind_rows(pred_2050_ssp2, pred_2050_ssp5, pred_2080_ssp2, pred_2080_ssp5)

# Define colors for each model
model_colors <- c("2050 SSP2" = "#A0D8BE", "2050 SSP5" = "#FDB774", "2080 SSP2" = "#5989BE", "2080 SSP5" = "#D73027")

png("plots/CEI_lat.png", width = 3000, height = 3000, res = 1000)

# Plotting with customized font sizes
CEI_latitude <- ggplot(all_predictions_CEI_lat, aes(x = mean_latitude, y = fit, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +  # Removes the edge around the ribbon
  scale_color_manual(values = model_colors) +  # Manually set line colors for each model
  scale_fill_manual(values = model_colors) +   # Manually set fill colors for each model
  labs(
    #title = "Predicted CEI by Mean Latitude",
    x = "Mean Latitude (°N)",
    y = "Predicted CEI",
    color = "Model",
    fill = "Model"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(2, 6)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),         # Set title font size
    #axis.title.x = element_text(size = 14),       # Set x-axis title font size
    axis.title.y = element_text(size = 14),       # Set y-axis title font size
    panel.grid.minor = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_text(size = 12),        # Set x-axis text font size
    axis.text.y = element_text(size = 12),        # Set y-axis text font size
    legend.title = element_text(size = 13),       # Set legend title font size
    legend.text = element_text(size = 11))         # Set legend text font size

# show graph
CEI_latitude

dev.off()  # Closes the device

# ------------------------- SAVE & JOIN DIAGNOSTICS TABLES -------------------------

library(dplyr)
library(tidyr)

# -------------------------
# 1. Extract model tables
# -------------------------

t_2050_ssp2 <- create_CEI_model_and_predict(df_2050_ssp2)$model_selection
t_2050_ssp5 <- create_CEI_model_and_predict(df_2050_ssp5)$model_selection
t_2080_ssp2 <- create_CEI_model_and_predict(df_2080_ssp2)$model_selection
t_2080_ssp5 <- create_CEI_model_and_predict(df_2080_ssp5)$model_selection

# -------------------------
# 
# -------------------------


t_2050_ssp2 <- create_CEI_model_and_predict(df_2050_ssp2)$model_selection %>%
  mutate(Scenario = "SSP2 2050")

t_2050_ssp5 <- create_CEI_model_and_predict(df_2050_ssp5)$model_selection %>%
  mutate(Scenario = "SSP5 2050")

t_2080_ssp2 <- create_CEI_model_and_predict(df_2080_ssp2)$model_selection %>%
  mutate(Scenario = "SSP2 2080")

t_2080_ssp5 <- create_CEI_model_and_predict(df_2080_ssp5)$model_selection %>%
  mutate(Scenario = "SSP5 2080")

cei_model_selection_table <- bind_rows(
  t_2050_ssp2,
  t_2050_ssp5,
  t_2080_ssp2,
  t_2080_ssp5
) %>%
  arrange(Scenario, delta_AIC) %>%
  dplyr::select(
    Scenario,
    Structure,
    k,
    logLik,
    delta_AIC,
    weight, 
    R2
  )

cei_model_selection_table <- cei_model_selection_table %>%
  mutate(
    logLik = round(logLik, 2),
    delta_AIC = round(delta_AIC, 2),
    weight = round(weight, 2), 
    R_square = R2
  )

write.csv(
  cei_model_selection_table,
  "plots/Table_S3.csv",
  row.names = FALSE
)



# -------------------------   CEI by LONGITUDE ----------------------------------

create_CEI_model_and_predict <- function(df) {
  
  library(MuMIn)
  library(dplyr)
  
  # Required by MuMIn
  old_na <- getOption("na.action")
  options(na.action = "na.fail")
  on.exit(options(na.action = old_na))
  
  # -----------------------------
  # Global model
  # -----------------------------
  global_model <- lm(
    CEI ~ guild + habitat + order +
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    data = df
  )
  
  # -----------------------------
  # Model selection
  # -----------------------------
  dredge_table <- dredge(
    global_model,
    rank = "AIC",
    subset =
      dc(mean_latitude, I(mean_latitude^2)) &
      dc(mean_longitude, I(mean_longitude^2))
  )
  
  # Convert to data frame for saving / inspection
  model_comparison <- as.data.frame(dredge_table) %>%
    transmute(
      # ---- inclusion indicators ----
      guild   = as.integer(!is.na(guild)),
      order   = as.integer(!is.na(order)),
      habitat = as.integer(!is.na(habitat)),
      
      lat = as.integer(
        !is.na(mean_latitude) | !is.na(`I(mean_latitude^2)`)
      ),
      
      long = as.integer(
        !is.na(mean_longitude) | !is.na(`I(mean_longitude^2)`)
      ),
      
      # ---- diagnostics ----
      weight     = weight,
      delta_AIC = delta,
      df        = df
    )  
  
  # -----------------------------
  # Best model
  # -----------------------------
  CEI_model <- get.models(dredge_table, 1)[[1]]
  
  # -----------------------------
  # Diagnostics
  # -----------------------------
  par(mfrow = c(2, 2))
  plot(CEI_model, which = c(1:3, 5))
  par(mfrow = c(1, 1))
  
  # -----------------------------
  # Prediction grid
  # -----------------------------
  longitude_range <- seq(
    min(df$mean_longitude, na.rm = TRUE),
    max(df$mean_longitude, na.rm = TRUE),
    length.out = 100
  )
  
  longitude_range <- seq(min(df$mean_longitude, na.rm = TRUE), max(df$mean_longitude, na.rm = TRUE), length.out = 100)
  
  # Create a data frame for prediction
  prediction_data <- data.frame(
    guild = "Foliage Gleaner",   # Change as needed for specific guild or other default values
    habitat = "Forests", # Setting to first level as a placeholder, adapt as required
    order = "Passerine",           # Change as needed for specific order or other default values
    mean_latitude = mean(df$mean_latitude, na.rm = TRUE),
    mean_longitude = longitude_range
  )
  
  preds <- predict(
    CEI_model,
    newdata = prediction_data,
    se.fit = TRUE
  )
  
  prediction_data <- prediction_data %>%
    mutate(
      fit   = preds$fit,
      lower = preds$fit - 1.96 * preds$se.fit,
      upper = preds$fit + 1.96 * preds$se.fit
    )
  
  # -----------------------------
  # Return
  # -----------------------------
  return(prediction_data)
}

# Apply the function to each data frame and store results
pred_2050_ssp2 <- create_CEI_model_and_predict(df_2050_ssp2)
pred_2050_ssp5 <- create_CEI_model_and_predict(df_2050_ssp5)
pred_2080_ssp2 <- create_CEI_model_and_predict(df_2080_ssp2)
pred_2080_ssp5 <- create_CEI_model_and_predict(df_2080_ssp5)

# Combine predictions into one data frame with a model identifier
pred_2050_ssp2$model <- "2050 SSP2"
pred_2050_ssp5$model <- "2050 SSP5"
pred_2080_ssp2$model <- "2080 SSP2"
pred_2080_ssp5$model <- "2080 SSP5"

all_predictions_CEI_long <- bind_rows(pred_2050_ssp2, pred_2050_ssp5, pred_2080_ssp2, pred_2080_ssp5)

png("plots/CEI_long.png", width = 3000, height = 3000, res = 1000)


# Plotting with customized font sizes for Longitude
CEI_longitude <- ggplot(all_predictions_CEI_long, aes(x = mean_longitude, y = fit, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +  # Removes the edge around the ribbon
  scale_color_manual(values = model_colors) +  # Manually set line colors for each model
  scale_fill_manual(values = model_colors) +   # Manually set fill colors for each model
  labs(
    #title = "Predicted CEI by Mean Longitude",
    x = "Mean Longitude (°W)",
    y = "Predicted CEI",
    color = "Model",
    fill = "Model"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(2, 6)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),         # Set title font size
    #axis.title.x = element_text(size = 14),       # Set x-axis title font size
    #axis.title.y = element_text(size = 14),       # Set y-axis title font size
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 12),        # Set x-axis text font size
    axis.text.y = element_text(size = 12),        # Set y-axis text font size
    legend.title = element_text(size = 13),       # Set legend title font size
    legend.text = element_text(size = 11)         # Set legend text font size
  )

# Show graph
CEI_longitude

dev.off()



# -------------------------   LEI by LATITUDE ----------------------------------


create_LEI_model_and_predict <- function(df) {
  
  library(MuMIn)
  library(dplyr)
  
  # Required by MuMIn
  old_na <- getOption("na.action")
  options(na.action = "na.fail")
  on.exit(options(na.action = old_na))
  
  # -----------------------------
  # Global model
  # -----------------------------
  global_model <- lm(
    LEI ~ guild + habitat + order +
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    data = df
  )
  
  # -----------------------------
  # Model selection
  # -----------------------------
  dredge_table <- dredge(
    global_model,
    rank = "AIC",
    subset =
      dc(mean_latitude, I(mean_latitude^2)) &
      dc(mean_longitude, I(mean_longitude^2))
  )
  
  model_comparison <- as.data.frame(dredge_table)
  
  model_comparison$R2 <- sapply(
    get.models(dredge_table, subset = TRUE),
    function(x) summary(x)$r.squared
  )
  
  # Convert to data frame for saving / inspection
  model_comparison <- model_comparison %>%
    mutate(
      Structure = paste0(
        ifelse(!is.na(guild),   "guild+",   ""),
        ifelse(!is.na(order),   "order+",   ""),
        ifelse(!is.na(habitat), "habitat+", ""),
        ifelse(
          !is.na(mean_latitude) | !is.na(`I(mean_latitude^2)`),
          "lat+", ""
        ),
        ifelse(
          !is.na(mean_longitude) | !is.na(`I(mean_longitude^2)`),
          "long+", ""
        )
      ),
      Structure = sub("\\+$", "", Structure)
    ) %>%
    transmute(
      Structure,
      k = df,
      logLik = logLik,
      delta_AIC = delta,
      weight = weight, 
      R2 = round(R2, 2)
    ) 
  
  # -----------------------------
  # Best model
  # -----------------------------
  LEI_model <- get.models(dredge_table, 1)[[1]]
  
  # -----------------------------
  # Diagnostics
  # -----------------------------
  par(mfrow = c(2, 2))
  plot(LEI_model, which = c(1:3, 5))
  par(mfrow = c(1, 1))
  
  # -----------------------------
  # Prediction grid
  # -----------------------------
  latitude_range <- seq(
    min(df$mean_latitude, na.rm = TRUE),
    max(df$mean_latitude, na.rm = TRUE),
    length.out = 100
  )
  
  prediction_data <- data.frame(
    guild = factor("Foliage Gleaner", levels = levels(df$guild)),
    habitat = factor("Forests", levels = levels(df$habitat)),
    order = factor("Passerine", levels = levels(df$order)),
    mean_latitude = latitude_range,
    mean_longitude = mean(df$mean_longitude, na.rm = TRUE)
  )
  
  preds <- predict(
    LEI_model,
    newdata = prediction_data,
    se.fit = TRUE
  )
  
  prediction_data <- prediction_data %>%
    mutate(
      fit   = preds$fit,
      lower = preds$fit - 1.96 * preds$se.fit,
      upper = preds$fit + 1.96 * preds$se.fit
    )
  
  # -----------------------------
  # Return
  # -----------------------------
  return(list(
    model = LEI_model,
    predictions = prediction_data,
    model_selection = model_comparison
  ))
}


# Apply the function to each data frame and store results
pred_2050_ssp2 <- create_LEI_model_and_predict(df_2050_ssp2)$predictions
pred_2050_ssp5 <- create_LEI_model_and_predict(df_2050_ssp5)$predictions
pred_2080_ssp2 <- create_LEI_model_and_predict(df_2080_ssp2)$predictions
pred_2080_ssp5 <- create_LEI_model_and_predict(df_2080_ssp5)$predictions

# Combine predictions into one data frame with a model identifier
pred_2050_ssp2$model <- "2050 SSP2"
pred_2050_ssp5$model <- "2050 SSP5"
pred_2080_ssp2$model <- "2080 SSP2"
pred_2080_ssp5$model <- "2080 SSP5"

all_predictions_LEI_lat <- bind_rows(pred_2050_ssp2, pred_2050_ssp5, pred_2080_ssp2, pred_2080_ssp5)


png("plots/LEI_lat.png", width = 3000, height = 3000, res = 1000)

# Plotting with customized font sizes
LEI_latitude <- ggplot(all_predictions_LEI_lat, aes(x = mean_latitude, y = fit, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +  # Removes the edge around the ribbon
  scale_color_manual(values = model_colors) +  # Manually set line colors for each model
  scale_fill_manual(values = model_colors) +   # Manually set fill colors for each model
  labs(
    #title = "Predicted LEI by Mean Latitude",
    x = "Mean Latitude (°N)",
    y = "Predicted LEI",
    color = "Model",
    fill = "Model"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.4)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),         # Set title font size
    axis.title.x = element_text(size = 14),       # Set x-axis title font size
    axis.title.y = element_text(size = 14),       # Set y-axis title font size
    #axis.title.x = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 12),        # Set x-axis text font size
    axis.text.y = element_text(size = 12),        # Set y-axis text font size
    legend.title = element_text(size = 13),       # Set legend title font size
    legend.text = element_text(size = 11)         # Set legend text font size
  )

# Show graph
LEI_latitude

dev.off()


# ------------------------- SAVE & JOIN DIAGNOSTICS TABLES (LEI) -------------------------

library(dplyr)
library(tidyr)

# -------------------------
# 1. Extract model tables
# -------------------------

t_2050_ssp2 <- create_LEI_model_and_predict(df_2050_ssp2)$model_selection %>%
  mutate(Scenario = "SSP2 2050")

t_2050_ssp5 <- create_LEI_model_and_predict(df_2050_ssp5)$model_selection %>%
  mutate(Scenario = "SSP5 2050")

t_2080_ssp2 <- create_LEI_model_and_predict(df_2080_ssp2)$model_selection %>%
  mutate(Scenario = "SSP2 2080")

t_2080_ssp5 <- create_LEI_model_and_predict(df_2080_ssp5)$model_selection %>%
  mutate(Scenario = "SSP5 2080")

lei_model_selection_table <- bind_rows(
  t_2050_ssp2,
  t_2050_ssp5,
  t_2080_ssp2,
  t_2080_ssp5
) %>%
  arrange(Scenario, delta_AIC) %>%
  dplyr::select(
    Scenario,
    Structure,
    k,
    logLik,
    delta_AIC,
    weight, 
    R2
  )

lei_model_selection_table <- lei_model_selection_table %>%
  mutate(
    logLik = round(logLik, 2),
    delta_AIC = round(delta_AIC, 2),
    weight = round(weight, 2), 
    R_squared = R2
  )

write.csv(
  lei_model_selection_table,
  "plots/Table_S4.csv",
  row.names = FALSE
)


# -------------------------   LEI by LONGITUDE ----------------------------------


create_LEI_model_and_predict <- function(df) {
  
  library(MuMIn)
  library(dplyr)
  
  # Required by MuMIn
  old_na <- getOption("na.action")
  options(na.action = "na.fail")
  on.exit(options(na.action = old_na))
  
  # -----------------------------
  # Global model
  # -----------------------------
  global_model <- lm(
    LEI ~ guild + habitat + order +
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    data = df
  )
  
  # -----------------------------
  # Model selection
  # -----------------------------
  dredge_table <- dredge(
    global_model,
    rank = "AIC",
    subset =
      dc(mean_latitude, I(mean_latitude^2)) &
      dc(mean_longitude, I(mean_longitude^2))
  )
  
  # Convert to data frame for saving / inspection
  model_comparison <- as.data.frame(dredge_table) %>%
    transmute(
      # ---- inclusion indicators ----
      guild   = as.integer(!is.na(guild)),
      order   = as.integer(!is.na(order)),
      habitat = as.integer(!is.na(habitat)),
      
      lat = as.integer(
        !is.na(mean_latitude) | !is.na(`I(mean_latitude^2)`)
      ),
      
      long = as.integer(
        !is.na(mean_longitude) | !is.na(`I(mean_longitude^2)`)
      ),
      
      # ---- diagnostics ----
      weight     = weight,
      delta_AIC = delta,
      df        = df
    )  
  
  # -----------------------------
  # Best model
  # -----------------------------
  LEI_model <- get.models(dredge_table, 1)[[1]]
  
  # -----------------------------
  # Diagnostics
  # -----------------------------
  par(mfrow = c(2, 2))
  plot(LEI_model, which = c(1:3, 5))
  par(mfrow = c(1, 1))
  
  # -----------------------------
  # Prediction grid
  # -----------------------------
  longitude_range <- seq(
    min(df$mean_longitude, na.rm = TRUE),
    max(df$mean_longitude, na.rm = TRUE),
    length.out = 100
  )
  
  longitude_range <- seq(min(df$mean_longitude, na.rm = TRUE), max(df$mean_longitude, na.rm = TRUE), length.out = 100)
  
  # Create a data frame for prediction
  prediction_data <- data.frame(
    guild = "Foliage Gleaner",   # Change as needed for specific guild or other default values
    habitat = "Forests", # Setting to first level as a placeholder, adapt as required
    order = "Passerine",           # Change as needed for specific order or other default values
    mean_latitude = mean(df$mean_latitude, na.rm = TRUE),
    mean_longitude = longitude_range
  )
  
  preds <- predict(
    LEI_model,
    newdata = prediction_data,
    se.fit = TRUE
  )
  
  prediction_data <- prediction_data %>%
    mutate(
      fit   = preds$fit,
      lower = preds$fit - 1.96 * preds$se.fit,
      upper = preds$fit + 1.96 * preds$se.fit
    )
  
  # -----------------------------
  # Return
  # -----------------------------
  return(prediction_data)
}

# Apply the function to each data frame and store results
pred_2050_ssp2 <- create_LEI_model_and_predict(df_2050_ssp2)
pred_2050_ssp5 <- create_LEI_model_and_predict(df_2050_ssp5)
pred_2080_ssp2 <- create_LEI_model_and_predict(df_2080_ssp2)
pred_2080_ssp5 <- create_LEI_model_and_predict(df_2080_ssp5)

# Combine predictions into one data frame with a model identifier
pred_2050_ssp2$model <- "2050 SSP2"
pred_2050_ssp5$model <- "2050 SSP5"
pred_2080_ssp2$model <- "2080 SSP2"
pred_2080_ssp5$model <- "2080 SSP5"

all_predictions_LEI_long <- bind_rows(pred_2050_ssp2, pred_2050_ssp5, pred_2080_ssp2, pred_2080_ssp5)

png("plots/LEI_long.png", width = 3000, height = 3000, res = 1000)

# Plotting with customized font sizes for Longitude
LEI_longitude <- ggplot(all_predictions_LEI_long, aes(x = mean_longitude, y = fit, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +  # Removes the edge around the ribbon
  scale_color_manual(values = model_colors) +  # Manually set line colors for each model
  scale_fill_manual(values = model_colors) +   # Manually set fill colors for each model
  labs(
    #title = "Predicted LEI by Mean Longitude",
    x = "Mean Longitude (°W)",
    y = "Predicted LEI",
    color = "Model",
    fill = "Model"
  ) +
  coord_cartesian(ylim = c(0, 0.4)) + 
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),         # Set title font size
    axis.title.x = element_text(size = 14),       # Set x-axis title font size
    #axis.title.y = element_text(size = 14),       # Set y-axis title font size
    axis.title.y = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 12),        # Set x-axis text font size
    axis.text.y = element_text(size = 12),        # Set y-axis text font size
    legend.title = element_text(size = 13))       # Set legend title font size
    
LEI_longitude
dev.off()


# -------------------------   LAT/LONG FULL PLOT  ---------------------------


# Load images
CEI_lat <- rasterGrob(readPNG("plots/CEI_lat.png"), interpolate = TRUE)
CEI_long <- rasterGrob(readPNG("plots/CEI_long.png"), interpolate = TRUE)
LEI_lat <- rasterGrob(readPNG("plots/LEI_lat.png"), interpolate = TRUE)
LEI_long <- rasterGrob(readPNG("plots/LEI_long.png"), interpolate = TRUE)
glm_legend <- rasterGrob(readPNG("plots/glm_legend.png"), interpolate = TRUE)


# Create labeled plots using ggdraw
CEI_lat_labeled <- ggdraw() +
  draw_grob(CEI_lat) +
  draw_label("a", x = 0.05, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold")

CEI_long_labeled <- ggdraw() +
  draw_grob(CEI_long) +
  draw_label("b", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold")

LEI_lat_labeled <- ggdraw() +
  draw_grob(LEI_lat) +
  draw_label("c", x = 0.05, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold")

LEI_long_labeled <- ggdraw() +
  draw_grob(LEI_long) +
  draw_label("d", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold")

# arrange in a grid
plots_grid <- plot_grid(
  CEI_lat_labeled, CEI_long_labeled,
  LEI_lat_labeled, LEI_long_labeled,
  ncol = 2
)

# combine legend and grid
final_plot <- plot_grid(
  plots_grid, 
  ggdraw() + draw_grob(glm_legend), 
  ncol = 1, 
  rel_heights = c(1, 0.07) # Adjust relative heights of the grid and legend
)

# save final plot
ggsave("plots/Figure4.png", final_plot, width = 12, height = 10)






# -------------------------   CEI Predictions by Habitat for Each Scenario ---------------------------

# Define a function to generate predictions for a specific model and scenario
create_predictions_by_habitat <- function(df, model_name) {

  final_model <- lm(
    CEI ~ habitat + mean_latitude + I(mean_latitude^2) + 
      mean_longitude + I(mean_longitude^2), 
    data = df,)
  
  # Create a data frame with distinct habitat values and representative values for prediction
  habitat_df <- df %>%
    distinct(habitat) %>%
    mutate(
      mean_latitude = mean(df$mean_latitude, na.rm = TRUE),
      mean_longitude = mean(df$mean_longitude, na.rm = TRUE),
      guild = factor("Foliage Gleaner", levels = levels(df$guild)),
      order = factor("Passerine", levels = levels(df$order))
    )
  
  # Generate predictions with confidence intervals for each habitat
  habitat_predictions <- predict(final_model, newdata = habitat_df, se.fit = TRUE)
  habitat_df <- habitat_df %>%
    mutate(
      fit = (habitat_predictions$fit),
      lower = (habitat_predictions$fit - 1.96 * habitat_predictions$se.fit),
      upper = (habitat_predictions$fit + 1.96 * habitat_predictions$se.fit),
      model = model_name
    )
  
  return(habitat_df)
}

# Generate predictions for each scenario
pred_2050_ssp2 <- create_predictions_by_habitat(df_2050_ssp2, "2050 SSP2")
pred_2050_ssp5 <- create_predictions_by_habitat(df_2050_ssp5, "2050 SSP5")
pred_2080_ssp2 <- create_predictions_by_habitat(df_2080_ssp2, "2080 SSP2")
pred_2080_ssp5 <- create_predictions_by_habitat(df_2080_ssp5, "2080 SSP5")

all_predictions_df <- bind_rows(pred_2050_ssp2, pred_2050_ssp5, pred_2080_ssp2, pred_2080_ssp5)
predictions_2050_df <- bind_rows(pred_2050_ssp2, pred_2050_ssp5)
predictions_2080_df <- bind_rows(pred_2080_ssp2, pred_2080_ssp5)

# Define colors for each model
model_colors <- c("2050 SSP2" = "#5989BE", "2050 SSP5" = "#D73027", "2080 SSP2" = "#5989BE", "2080 SSP5" = "#D73027")


# Plotting side-by-side points and whiskers by habitat with gray shading
CEI_habitat_plot_2050 <- ggplot(predictions_2050_df, aes(x = habitat, y = fit, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = model_colors) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black", size = 0.5) +  # Add y=1 line
  geom_segment(data = predictions_2050_df, aes(x = habitat, xend = habitat, y = 0.95, yend = 1.05), 
               inherit.aes = FALSE, color = "black", size = 0.3) + 
  geom_segment(aes(x = -Inf, xend = -Inf, y = 1, yend = 5), linetype = "solid", color = "black", size = 0.75) +  # Line at y=1
  labs(
    #title = "Predicted CEI by Habitat across Scenarios",
    x = "Habitat",
    y = "Predicted CEI",
    color = "Model"
  ) +
  theme_minimal() +
  ylim(0.95, 5) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),
    #axis.title.x = element_text(size = 14),
    axis.title.x = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.title.y = element_text(size = 14),
    #axis.text.x = element_blank(), 
    axis.text.x = element_text(size = 12, angle = 70, hjust = 1, vjust = 1.1),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13), 
    panel.grid.major.x = element_blank()
  )

CEI_habitat_plot_2050

ggsave("plots/CEI_habitat_2050.png", CEI_habitat_plot_2050, width = 3.5, height = 3.5)




CEI_habitat_plot_2080 <- ggplot(predictions_2080_df, aes(x = habitat, y = fit, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = model_colors) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black", size = 0.5) +  # Add y=1 line
  geom_segment(data = predictions_2050_df, aes(x = habitat, xend = habitat, y = 0.95, yend = 1.05), 
               inherit.aes = FALSE, color = "black", size = 0.3) + 
  geom_segment(aes(x = -Inf, xend = -Inf, y = 1, yend = 5), linetype = "solid", color = "black", size = 0.75) +  # Line at y=1
  labs(
    #title = "Predicted CEI by Habitat across Scenarios",
    x = "Habitat",
    y = "Predicted CEI",
    color = "Model"
  ) +
  theme_minimal() +
  ylim(0.95, 5) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),
    #axis.title.x = element_text(size = 14),
    axis.title.x = element_blank(), 
    #axis.title.y = element_text(size = 14),
    panel.grid.minor = element_blank(), 
    axis.title.y = element_blank(), 
    #axis.text.x = element_blank(), 
    axis.text.x = element_text(size = 12, angle = 70, hjust = 1, vjust = 1.1),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13),
    panel.grid.major.x = element_blank()
  )

CEI_habitat_plot_2080

ggsave("plots/CEI_habitat_2080.png", CEI_habitat_plot_2080, width = 3.5, height = 3.5)







# -------------------------   LEI Predictions by Habitat for Each Scenario ---------------------------

# Define a function to generate predictions for a specific model and scenario
create_predictions_by_habitat <- function(df, model_name) {
  
  final_model <- lm(
    LEI ~ habitat + mean_latitude + I(mean_latitude^2) + 
      mean_longitude + I(mean_longitude^2), 
    data = df,)
  
  # Create a data frame with distinct habitat values and representative values for prediction
  habitat_df <- df %>%
    distinct(habitat) %>%
    mutate(
      mean_latitude = mean(df$mean_latitude, na.rm = TRUE),
      mean_longitude = mean(df$mean_longitude, na.rm = TRUE),
      guild = factor("Foliage Gleaner", levels = levels(df$guild)),
      order = factor("Passerine", levels = levels(df$order))
    )
  
  # Generate predictions with confidence intervals for each habitat
  habitat_predictions <- predict(final_model, newdata = habitat_df, se.fit = TRUE)
  habitat_df <- habitat_df %>%
    mutate(
      fit = (habitat_predictions$fit),
      lower = (habitat_predictions$fit - 1.96 * habitat_predictions$se.fit),
      upper = (habitat_predictions$fit + 1.96 * habitat_predictions$se.fit),
      model = model_name
    )
  
  return(habitat_df)
}

# Generate predictions for each scenario
pred_2050_ssp2 <- create_predictions_by_habitat(df_2050_ssp2, "2050 SSP2")
pred_2050_ssp5 <- create_predictions_by_habitat(df_2050_ssp5, "2050 SSP5")
pred_2080_ssp2 <- create_predictions_by_habitat(df_2080_ssp2, "2080 SSP2")
pred_2080_ssp5 <- create_predictions_by_habitat(df_2080_ssp5, "2080 SSP5")

# Combine all predictions
all_predictions_df <- bind_rows(pred_2050_ssp2, pred_2050_ssp5, pred_2080_ssp2, pred_2080_ssp5)
predictions_2050_df <- bind_rows(pred_2050_ssp2, pred_2050_ssp5)
predictions_2080_df <- bind_rows(pred_2080_ssp2, pred_2080_ssp5)

# Plotting side-by-side points and whiskers by habitat
LEI_habitat_plot_2050 <- ggplot(predictions_2050_df, aes(x = habitat, y = fit, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = model_colors) +
  geom_hline(yintercept = 0.2, linetype = "solid", color = "black", size = 0.5) +  # Add y=1 line
  geom_segment(data = predictions_2050_df, aes(x = habitat, xend = habitat, y = 0.195, yend = 0.205), 
               inherit.aes = FALSE, color = "black", size = 0.3) + 
  geom_segment(aes(x = -Inf, xend = -Inf, y = 0.2, yend = 0.6), linetype = "solid", color = "black", size = 0.75) +  # Line at y=1
  labs(
    #title = "Predicted LEI by Habitat across Scenarios",
    x = "Habitat (2050's)",
    y = "Predicted LEI",
    color = "Model"
  ) +
  theme_minimal() +
  ylim(0.195, 0.6) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),
    #axis.title.x = element_text(size = 14),
    axis.title.x = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.grid.major.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 70, hjust = 1, vjust = 1.1),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

LEI_habitat_plot_2050

ggsave("plots/LEI_habitat_2050.png", LEI_habitat_plot_2050, width = 3.5, height = 3.5)


# Plotting side-by-side points and whiskers by habitat
LEI_habitat_plot_2080 <- ggplot(predictions_2080_df, aes(x = habitat, y = fit, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = model_colors) +
  labs(
    #title = "Predicted LEI by Habitat across Scenarios",
    x = "Habitat (2080's)",
    y = "Predicted LEI",
    color = "Model"
  ) +
  theme_minimal() +
  ylim(0.195, 0.6) + 
  geom_hline(yintercept = 0.2, linetype = "solid", color = "black", size = 0.5) +  # Add y=1 line
  geom_segment(data = predictions_2050_df, aes(x = habitat, xend = habitat, y = 0.195, yend = 0.205), 
               inherit.aes = FALSE, color = "black", size = 0.3) + 
  geom_segment(aes(x = -Inf, xend = -Inf, y = 0.2, yend = 0.6), linetype = "solid", color = "black", size = 0.75) +  # Line at y=1
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18),
    #axis.title.x = element_text(size = 14),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 12, angle = 70, hjust = 1, vjust = 1.1),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

LEI_habitat_plot_2080

ggsave("plots/LEI_habitat_2080.png", LEI_habitat_plot_2080, width = 3.5, height = 3.5)






# -------------------------   HABITAT FULL PLOT  -------------------------------



# Load images
CEI_habitat_2050 <- rasterGrob(readPNG("plots/CEI_habitat_2050.png"), interpolate = TRUE)
CEI_habitat_2080 <- rasterGrob(readPNG("plots/CEI_habitat_2080.png"), interpolate = TRUE)
LEI_habitat_2050 <- rasterGrob(readPNG("plots/LEI_habitat_2050.png"), interpolate = TRUE)
LEI_habitat_2080 <- rasterGrob(readPNG("plots/LEI_habitat_2080.png"), interpolate = TRUE)
glm_legend <- rasterGrob(readPNG("plots/glm_legend_reduced.png"), interpolate = TRUE)


# Create labeled plots using ggdraw
CEI_lat_labeled <- ggdraw() +
  draw_grob(CEI_habitat_2050) +
  draw_label("a", x = 0.05, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold") + 
  draw_label("2050", x = 0.45, y = 1, hjust = 0, vjust = 1, size = 25, fontface = "bold")

CEI_long_labeled <- ggdraw() +
  draw_grob(CEI_habitat_2080) +
  draw_label("b", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold") + 
  draw_label("2080", x = 0.45, y = 1, hjust = 0, vjust = 1, size = 25, fontface = "bold")

LEI_lat_labeled <- ggdraw() +
  draw_grob(LEI_habitat_2050) +
  draw_label("c", x = 0.05, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold") +
  draw_label("2050", x = 0.45, y = 1, hjust = 0, vjust = 1, size = 25, fontface = "bold")

LEI_long_labeled <- ggdraw() +
  draw_grob(LEI_habitat_2080) +
  draw_label("d", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 25, fontface = "bold") + 
  draw_label("2080", x = 0.45, y = 1, hjust = 0, vjust = 1, size = 25, fontface = "bold")


# arrange in a grid
plots_grid <- plot_grid(
  CEI_lat_labeled, CEI_long_labeled,
  LEI_lat_labeled, LEI_long_labeled,
  ncol = 2
)

# combine legend and grid
final_plot <- plot_grid(
  plots_grid, 
  ggdraw() + draw_grob(glm_legend), 
  ncol = 1, 
  rel_heights = c(1, 0.05) # Adjust relative heights of the grid and legend
)

# save final plot
ggsave("plots/Figure5.png", final_plot, width = 12, height = 12)


