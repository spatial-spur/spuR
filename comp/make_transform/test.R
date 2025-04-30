# Source required functions
source("R/get_s_matrix.R")
source("R/make_transform.R")
source("R/getdistmat_euclidian.R")
source("R/getdistmat_lat_lon.R")
source("R/get_sigma_lbm.R")
source("R/get_sigma_lbm_dm.R")
source("R/demean_sigma.R")

# Create sample data - US cities with lat/long coordinates
cities <- data.frame(
  city = c("New York", "Los Angeles", "Chicago", "Houston", "Phoenix", "Philadelphia"),
  s_1 = c(40.7128, 34.0522, 41.8781, 29.7604, 33.4484, 39.9526),  # latitude
  s_2 = c(-74.0060, -118.2437, -87.6298, -95.3698, -112.0740, -75.1652)  # longitude
)

# Print the input data
cat("==== Input Data ====\n")
print(cities)

# Extract coordinates
coords <- cities[, c("s_1", "s_2")]

# Set output precision for better comparison
options(digits = 8)

# Test 1: LBM-GLS transformation (default)
cat("\n==== Testing LBM-GLS transformation ====\n")
H_lbmgls <- make_transform(coords, transformation = "lbmgls", latlong = TRUE)
cat("Dimension:", dim(H_lbmgls), "\n")
print("First 3x3 submatrix:")
print(H_lbmgls[1:3, 1:3])
cat("Matrix trace:", sum(diag(H_lbmgls)), "\n")
cat("Matrix determinant:", det(H_lbmgls), "\n")
