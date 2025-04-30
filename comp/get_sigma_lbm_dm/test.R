# Source the required functions
source("R/get_sigma_lbm.R")
source("R/demean_sigma.R")
source("R/get_sigma_lbm_dm.R")

# Create a test distance matrix
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

# Get the demeaned sigma matrix
sigma_lbm_dm <- get_sigma_lbm_dm(distmat)

# Print the result
cat("\nDemeaned sigma_lbm matrix (R implementation):\n")
print(sigma_lbm_dm)

# Check properties of demeaned matrix
cat("\nMatrix properties check:\n")
cat("Matrix dimension:", dim(sigma_lbm_dm), "\n")
cat("Row sums (should be ~0):", round(rowSums(sigma_lbm_dm), 10), "\n")
cat("Column sums (should be ~0):", round(colSums(sigma_lbm_dm), 10), "\n")
cat("Trace:", sum(diag(sigma_lbm_dm)), "\n")

# Print specific elements for easier comparison
cat("\nSpecific elements for comparison:\n")
cat("[1,1]:", sigma_lbm_dm[1,1], "\n")
cat("[2,3]:", sigma_lbm_dm[2,3], "\n")
cat("[4,5]:", sigma_lbm_dm[4,5], "\n")