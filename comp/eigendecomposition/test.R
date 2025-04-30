# Source the required functions
source("R/get_sigma_lbm.R")
source("R/demean_sigma.R")
source("R/get_sigma_lbm_dm.R")

# Set higher precision for output
options(digits = 15)

# Create a small test distance matrix
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

# Test 1: Get demeaned sigma_lbm matrix (direct function)
cat("\nTest 1: get_sigma_lbm_dm function\n")
sigma_lbm_dm <- get_sigma_lbm_dm(distmat)
cat("Result dimension:", dim(sigma_lbm_dm), "\n")
cat("Matrix trace:", sum(diag(sigma_lbm_dm)), "\n")
cat("Row sums (should be ~0):", rowSums(sigma_lbm_dm), "\n")
cat("Full matrix:\n")
print(sigma_lbm_dm)

# Test 2: Compare with two-step approach
cat("\nTest 2: Compare with two-step approach\n")
sigma_lbm <- get_sigma_lbm(distmat)
sigma_lbm_dm2 <- demean_sigma(sigma_lbm)
cat("Max difference between methods:", max(abs(sigma_lbm_dm - sigma_lbm_dm2)), "\n")

# Test 3: Check eigendecomposition (for comparison with Stata)
cat("\nTest 3: Check eigendecomposition\n")
eig <- eigen(sigma_lbm_dm, symmetric = TRUE)
values <- eig$values
vectors <- eig$vectors

# Sort eigenvalues decreasing (to match Stata's order() function)
idx <- order(values, decreasing = TRUE)
values <- values[idx]
vectors <- vectors[, idx]

# Apply same threshold as Stata
small <- 1e-10
idx <- values > small
values <- values[idx]
vectors <- vectors[, idx]

# Print eigenvalues for comparison
cat("Number of eigenvalues above threshold:", length(values), "\n")
cat("Eigenvalues (descending):", values, "\n")
cat("First eigenvector:", vectors[,1], "\n")
cat("Last eigenvector:", vectors[,length(values)], "\n")

# Test 4: Create and verify full GLS transformation matrix
dsi <- 1 / sqrt(values)
H <- vectors %*% diag(dsi) %*% t(vectors)

cat("\nTest 4: Final transformation matrix properties\n")
cat("Matrix shape:", dim(H), "\n")
cat("Row sums:", rowSums(H)[1:5], "\n")  # First 5 row sums
cat("First element [1,1]:", H[1,1], "\n")
cat("Element [2,3]:", H[2,3], "\n")
cat("Element [4,5]:", H[4,5], "\n")