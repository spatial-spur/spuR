

source("R/get_s_matrix.R")
# Create sample data - US cities with lat/long coordinates
cities <- data.frame(
  city = c("New York", "Los Angeles", "Chicago", "Houston", "Phoenix", "Philadelphia"),
  s_1 = c(40.7128, 34.0522, 41.8781, 29.7604, 33.4484, 39.9526),  # latitude
  s_2 = c(-74.0060, -118.2437, -87.6298, -95.3698, -112.0740, -75.1652)  # longitude
)

# Print the input data
cat("==== Input Data ====\n")
print(cities)

# Extract just the coordinates
coords <- cities[, c("s_1", "s_2")]

# Test with latlong = FALSE (Euclidean distance)
cat("\n==== Testing with latlong = FALSE ====\n")
result_euclidean <- get_s_matrix(coords, latlong = FALSE)
print("Return value:")
print(head(result_euclidean))
cat("Global flag value:", latlongflag, "\n")

# Test with latlong = TRUE (Spherical distance)
cat("\n==== Testing with latlong = TRUE ====\n")
result_spherical <- get_s_matrix(coords, latlong = TRUE)
print("Return value:")
print(head(result_spherical))
cat("Global flag value:", latlongflag, "\n")
