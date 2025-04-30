# Test script for spatial transformation functions using Chetty data

library(dplyr)
library(here)
library(janitor)
library(fixest)
library(geodist)

# If using development version, source the required functions
if (!requireNamespace("spuR", quietly = TRUE)) {
  source_dir <- "c:/Users/dgoett/Documents/spuR/R"
  source(file.path(source_dir, "transform.R"))
  source(file.path(source_dir, "make_transform.R"))
  source(file.path(source_dir, "getdistmat_euclidian.R"))
  source(file.path(source_dir, "getdistmat_lat_lon.R"))
  source(file.path(source_dir, "lvech.R"))
  source(file.path(source_dir, "get_sigma_lbm.R"))
  source(file.path(source_dir, "get_sigma_lbm_dm.R"))
  source(file.path(source_dir, "demean_sigma.R"))
  source(file.path(source_dir, "get_s_matrix.R"))
}

# Load Chetty data from .Rda file
load(here("data", "chetty.Rda"))  
chetty <- clean_names(chetty)
chetty = chetty %>% 
    filter(!(state %in% c("HI", "AK"))) %>%
    mutate(
        am = scale(am),
        frac_black = scale(frac_black),
    ) %>%
    filter(
        !is.na(am) & !is.na(frac_black) 
    ) 
print(head(chetty))  # Check the structure of the loaded data

# Apply spatial transformation (using LBM-GLS as default)
transformed_data <- spur_transform(
  data = chetty,
  vars = c("am", "frac_black"),
  lat = "lat",
  lon = "lon",
  transformation = "lbmgls",
)

# Run regression on original variables
print(summary(feols(tram ~ trfrac_black - 1, data = transformed_data, se = "HC1")))
