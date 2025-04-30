# Source required functions
source("R/make_transform.R")
source("R/getdistmat_euclidian.R")
source("R/getdistmat_lat_lon.R")
source("R/get_sigma_lbm.R")
source("R/get_sigma_lbm_dm.R")
source("R/demean_sigma.R")

# Set higher precision for output
options(digits = 15)

# Create test dataset with spatial coordinates
locations <- rbind(
  c(40.7128, -74.0060),   # New York
  c(34.0522, -118.2437),  # Los Angeles
  c(41.8781, -87.6298),   # Chicago
  c(29.7604, -95.3698),   # Houston
  c(25.7617, -80.1918)    # Miami
)

# Add row and column names for clarity
rownames(locations) <- c("New York", "Los Angeles", "Chicago", "Houston", "Miami")
colnames(locations) <- c("lat", "lon")

# Print input data
cat("Input coordinates:\n")
print(locations)

# Test 1: Using Euclidean distances (latlong = FALSE)
cat("\nTest 1: lbm_gls_matrix with Euclidean distances\n")
H_euclidean <- lbm_gls_matrix(locations, latlong = FALSE)

# Print results and matrix properties
cat("Matrix dimension:", dim(H_euclidean), "\n")
cat("Matrix trace:", sum(diag(H_euclidean)), "\n")
cat("First 3x3 submatrix:\n")
print(H_euclidean)

# Test 2: Using lat/long distances (latlong = TRUE)
cat("\nTest 2: lbm_gls_matrix with lat/long distances\n")
H_latlong <- lbm_gls_matrix(locations, latlong = TRUE)

# Print results and matrix properties
cat("Matrix dimension:", dim(H_latlong), "\n")
cat("Matrix trace:", sum(diag(H_latlong)), "\n")
print(H_latlong)


cat("\nResults saved to r_lbm_gls_euclidean.csv and r_lbm_gls_latlong.csv\n")