# ------------------------------------------------------------------------------

#                        Coefficient Estimate tables for Models

# ------------------------------------------------------------------------------

# Script to calculate model coefficients and statistics for top linear models 
# in predicting CEI and LEI by species ecological traits. 

# Last modified by: Matthew Gilbert
# Last modified on: 23 Aug 2026

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


# --------------------- LIST OF SUPPORTED MODELS ----------------------------


supported_models = c("LEI SSP2 2050 habitat+lat+long", 
"LEI SSP2 2050 guild+lat+long", 
"LEI SSP2 2050 order+habitat+lat+long", 
"LEI SSP2 2080 habitat+lat+long", 
"LEI SSP2 2080 guild+lat+long", 
"LEI SSP2 2080 order+habitat+lat+long", 
"LEI SSP5 2050 habitat+lat+long", 
"LEI SSP5 2050 order+habitat+lat+long", 
"LEI SSP5 2080 habitat+lat+long", 
"LEI SSP5 2080 guild+lat+long", 
"LEI SSP5 2080 order+habitat+lat+long", 
"CEI SSP2 2050 habitat+lat+long", 
"CEI SSP2 2050 order+habitat+lat+long", 
"CEI SSP2 2080 habitat+lat+long", 
"CEI SSP2 2080 order+habitat+lat+long", 
"CEI SSP5 2050 habitat+lat+long", 
"CEI SSP5 2050 order+habitat+lat+long", 
"CEI SSP5 2080 habitat+lat+long", 
"CEI SSP5 2080 order+habitat+lat+long")

# --------------------- FIT+ESTIMATE MODEL COEFFICIENTS ----------------------------


library(dplyr)
library(broom)
library(stringr)
library(purrr)

#----------------------------------------
# Match scenario names to data frames
#----------------------------------------

datasets <- list(
  "SSP2 2050" = df_2050_ssp2,
  "SSP5 2050" = df_2050_ssp5,
  "SSP2 2080" = df_2080_ssp2,
  "SSP5 2080" = df_2080_ssp5
)

#----------------------------------------
# Function to fit one supported model
#----------------------------------------

fit_supported_model <- function(model_string){
  
  # Parse model description
  parts <- str_split(model_string, " ", simplify = TRUE)
  
  response <- parts[1]
  scenario <- paste(parts[2], parts[3])
  structure <- parts[4]
  
  df <- datasets[[scenario]]
  
  # Build formula
  predictors <- c()
  
  if(str_detect(structure, "guild"))
    predictors <- c(predictors, "guild")
  
  if(str_detect(structure, "order"))
    predictors <- c(predictors, "order")
  
  if(str_detect(structure, "habitat"))
    predictors <- c(predictors, "habitat")
  
  if(str_detect(structure, "lat"))
    predictors <- c(
      predictors,
      "mean_latitude",
      "I(mean_latitude^2)"
    )
  
  if(str_detect(structure, "long"))
    predictors <- c(
      predictors,
      "mean_longitude",
      "I(mean_longitude^2)"
    )
  
  formula <- as.formula(
    paste(response,
          "~",
          paste(predictors, collapse = " + "))
  )
  
  model <- lm(formula, data = df)
  
  delta_AIC <- AIC(model)
  
  broom::tidy(model) %>%
    transmute(
      Response = response,
      Scenario = scenario,
      Model_Structure = structure,
      delta_AIC = delta_AIC,
      Term = term,
      Estimate = estimate,
      SE = std.error,
      t = statistic,
      p_value = p.value
    )
}

#----------------------------------------
# Create complete coefficient table
#----------------------------------------

coefficient_table <-
  map_dfr(supported_models, fit_supported_model)

write.csv(
  coefficient_table,
  "plots/Table_S5.csv",
  row.names = FALSE
)
