# Source the function to test
source("R/get_sigma_lbm.R")

# Set higher precision for output
options(digits = 15)

# Create a small test distance matrix
# Using a 5x5 matrix with some arbitrary distances
distmat <- matrix(c(
  0.0000, 0.6124, 0.2171, 0.4829, 0.8331,
  0.6124, 0.0000, 0.5251, 0.3798, 0.7481,
  0.2171, 0.5251, 0.0000, 0.3091, 0.6509,
  0.4829, 0.3798, 0.3091, 0.0000, 0.3553,
  0.8331, 0.7481, 0.6509, 0.3553, 0.0000
), nrow = 5, byrow = TRUE)

# Print the input distance matrix
cat("Input distance matrix:\n")
print(distmat)

# Compute the LBM covariance matrix
sigma_lbm <- get_sigma_lbm(distmat)

# Print the result
cat("\nComputed LBM covariance matrix (R implementation):\n")
print(sigma_lbm)

# Print specific values for easier comparison
cat("\nKey values for comparison:\n")
cat("First element [1,1]:", sigma_lbm[1,1], "\n")
cat("Element [2,3]:", sigma_lbm[2,3], "\n")
cat("Element [4,5]:", sigma_lbm[4,5], "\n")
cat("Matrix trace:", sum(diag(sigma_lbm)), "\n")