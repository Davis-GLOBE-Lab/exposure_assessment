# ------------------------------------------------------------------------------

# -----------    MODEL FAMILY COMPARISON USING MODEL DIAGNOSTICS    ------------

# ------------------------------------------------------------------------------

# Compares 3 model families (gaussian linear, log-transformed linear, and gamma) 
# for CEI and 2 model families (gaussian linear, beta) for LEI using DHARMa 
# model diagnostics

# Last modified by: Matthew Gilbert
# Last modified on: 28 Jan 2026

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

# ------------------------------------------------------------------------------
# ------------------   DATA PREPROCESSING  -----------------------
# ------------------------------------------------------------------------------


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

# ------------------------------------------------------------------------------
# ------------------   DEFINE AND RUN MODELS for CEI  -----------------------
# ------------------------------------------------------------------------------

#### LINEAR (CEI) --------------------------------------------------

linear_model <- function(df, scenario) {
  
  # Set model
  CEI_model <- lm(
    CEI ~ habitat + 
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    data = df
  )
  
  # DHARMa diagnostics
  simulationOutput <- simulateResiduals(
    fittedModel = CEI_model,
    plot = FALSE
  )
  plot(simulationOutput, title = NA)
}

#### GAMMA (CEI) --------------------------------------------------

gamma_model <- function(df, scenario) {
  
  # Set model
  CEI_model <- glmmTMB(
    CEI ~ habitat + 
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    family = Gamma(link = "log"),
    data = df
  )
  
  # DHARMa diagnostics
  simulationOutput <- simulateResiduals(
    fittedModel = CEI_model,
    plot = FALSE
  )
  plot(simulationOutput, title = NA)
}

#### LOG-TRANSFORM (CEI) --------------------------------------------------

log_transform_model <- function(df, scenario) {
  
  # Set model
  CEI_model <- glmmTMB(
    log(CEI) ~ habitat +
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    family = gaussian(),
    data = df
  )
  
  # DHARMa diagnostics
  simulationOutput <- simulateResiduals(
    fittedModel = CEI_model,
    plot = T
  )
  plot(simulationOutput, title = NA)
}


#### Run different models ------------------------------------------------------

# 2050 SSP2
png("plots/Figure_S3a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2050_ssp2, "2050 SSP2")
dev.off()

png("plots/Figure_S3b.png", width = 5400, height = 3000, res = 600)
gamma_model(df_2050_ssp2, "2050 SSP2")
dev.off()

png("plots/Figure_S3c.png", width = 5400, height = 3000, res = 600)
log_transform_model(df_2050_ssp2, "2050 SSP2")
dev.off()

# 2050 SSP5
png("plots/Figure_S4a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2050_ssp5, "2050 SSP5")
dev.off()

png("plots/Figure_S4b.png", width = 5400, height = 3000, res = 600)
gamma_model(df_2050_ssp5, "2050 SSP5")
dev.off()

png("plots/Figure_S4c.png", width = 5400, height = 3000, res = 600)
log_transform_model(df_2050_ssp5, "2050 SSP5")
dev.off()

# 2080 SSP2
png("plots/Figure_S5a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2080_ssp2, "2080 SSP2")
dev.off()

png("plots/Figure_S5b.png", width = 5400, height = 3000, res = 600)
gamma_model(df_2080_ssp2, "2080 SSP2")
dev.off()

png("plots/Figure_S5c.png", width = 5400, height = 3000, res = 600)
log_transform_model(df_2080_ssp2, "2080 SSP2")
dev.off()

# 2080 SSP5
png("plots/Figure_S6a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2080_ssp5, "2080 SSP5")
dev.off()

png("plots/Figure_S6b.png", width = 5400, height = 3000, res = 600)
gamma_model(df_2080_ssp5, "2080 SSP5")
dev.off()

png("plots/Figure_S6c.png", width = 5400, height = 3000, res = 600)
log_transform_model(df_2080_ssp5, "2080 SSP5")
dev.off()


# ------------------------------------------------------------------------------
# ------------------   DEFINE AND RUN MODELS for LEI  -----------------------
# ------------------------------------------------------------------------------

#### LINEAR (LEI) --------------------------------------------------

linear_model <- function(df, scenario) {
  
  # Set model
  LEI_model <- lm(
    LEI ~ habitat + 
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    data = df
  )
  
  # DHARMa diagnostics
  simulationOutput <- simulateResiduals(
    fittedModel = LEI_model,
    plot = FALSE
  )
  plot(simulationOutput, title = paste0("LEI DHARMa Linear ", scenario))
}

#### BETA (LEI) --------------------------------------------------

beta_model <- function(df, scenario) {
  
  # Set model
  LEI_model <- glmmTMB(
    LEI ~ habitat + 
      mean_latitude + I(mean_latitude^2) +
      mean_longitude + I(mean_longitude^2),
    family = beta_family(link = "logit"),
    data = df
  )
  
  # DHARMa diagnostics
  simulationOutput <- simulateResiduals(
    fittedModel = LEI_model,
    plot = FALSE
  )
  plot(simulationOutput, title = paste0("LEI DHARMa Beta ", scenario))
}


#### Run different models ------------------------------------------------------

# 2050 SSP2
png("plots/Figure_S7a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2050_ssp2, "2050 SSP2")
dev.off()

png("plots/Figure_S7b.png", width = 5400, height = 3000, res = 600)
beta_model(df_2050_ssp2, "2050 SSP2")
dev.off()

# 2050 SSP5
png("plots/Figure_S8a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2050_ssp5, "2050 SSP5")
dev.off()

png("plots/Figure_S8b.png", width = 5400, height = 3000, res = 600)
beta_model(df_2050_ssp5, "2050 SSP5")
dev.off()

# 2080 SSP2
png("plots/Figure_S9a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2080_ssp2, "2080 SSP2")
dev.off()

png("plots/Figure_S9b.png", width = 5400, height = 3000, res = 600)
beta_model(df_2080_ssp2, "2080 SSP2")
dev.off()

# 2080 SSP5
png("plots/Figure_S10a.png", width = 5400, height = 3000, res = 600)
linear_model(df_2080_ssp5, "2080 SSP5")
dev.off()

png("plots/Figure_S10b.png", width = 5400, height = 3000, res = 600)
beta_model(df_2080_ssp5, "2080 SSP5")
dev.off()
























