# Create a simple test matrix in R
test_matrix <- matrix(c(
  4, 2, 1,
  2, 5, 3,
  1, 3, 6
), nrow=3, byrow=TRUE)

# Test the get_sigma_lbm and demean_sigma functions
dist_mat <- matrix(c(
  0, 1, 2,
  1, 0, 3,
  2, 3, 0
), nrow=3, byrow=TRUE)

# Run your functions
sigma_lbm <- get_sigma_lbm(dist_mat)
sigma_lbm_dm <- demean_sigma(sigma_lbm)

# Then test Get_R
r_result <- Get_R(sigma_lbm_dm, 2)
print("R values:")
print(r_result$DS)
print(r_result$R)

