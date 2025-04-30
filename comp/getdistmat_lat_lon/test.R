# Simple test for getdistmat_lat_lon function to compare with Stata output

# Source the function to test
source("R/getdistmat_lat_lon.R")

# Create a data frame with city coordinates (lat, lon)
# Using small set for easier comparison
cities <- data.frame(
  city = c("New York", "Los Angeles", "Chicago", "Houston", "Miami"),
  lat = c(40.7128, 34.0522, 41.8781, 29.7604, 25.7617),
  lon = c(-74.0060, -118.2437, -87.6298, -95.3698, -80.1918)
)

# Print the input data
cat("Input coordinates:\n")
print(cities)

# Convert to matrix format expected by the function
s <- as.matrix(cities[, c("lat", "lon")])

# Calculate distance matrix
dist_mat <- getdistmat_lat_lon(s)

# Convert to km for easier comparison with Stata
dist_mat_km <- dist_mat * (2 * pi * 6371)

# Print the full distance matrix with 6 decimal precision
cat("\nDistance matrix in kilometers (R implementation):\n")
print(round(dist_mat_km, 6))
