# Create test distance matrix
n <- 5
coords <- expand.grid(1:n, 1:n)
dist_mat <- as.matrix(dist(coords))

# Test different rho values
rho_vals <- c(0.001, 0.01, 0.05, 0.1)

# Run getcbar for each rho value
results <- sapply(rho_vals, function(rho) {
  getcbar(rho, dist_mat)
})

# Print results
cat("R getcbar results:\n")
print(data.frame(
  rho = rho_vals,
  cbar = results,
  row.names = NULL
))