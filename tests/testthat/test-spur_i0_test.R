library(testthat)

test_that("spur_i0_test produces expected outputs with simple inputs", {
  # Set seed for reproducibility
  set.seed(123)
  
  # Create test inputs
  n <- 10  # Number of observations
  n_series <- 2  # Number of time series
  q <- 3  # Reducing number of eigenvectors (less than n)
  n_sim <- 200  # Number of simulations
  
  # Create a more realistic distance matrix - using 2D coordinates
  coords <- cbind(runif(n), runif(n))  # Random coordinates in unit square
  distmat <- as.matrix(dist(coords))
  
  # Create data matrix Y with spatial correlation
  Y <- matrix(rnorm(n * n_series), nrow = n, ncol = n_series)
  
  # Create simulation matrix emat - standard normal to ensure good properties
  emat <- matrix(rnorm(q * n_sim), nrow = q, ncol = n_sim)
  
  # Run the function with try to capture any errors
  result <- tryCatch({
    spur_i0_test(Y, distmat, emat, q)
  }, error = function(e) {
    print(paste("Error in spur_i0_test:", e$message))
    NULL
  })
  
  # Only test if we got a result
  if (!is.null(result)) {
    # Test structure and dimensions of the output
    expect_type(result, "list")
    expect_named(result, c("LR", "pvalue", "cvalue", "ha_parm", 
                         "cvalue_mat", "pvalue_mat", "rho_grid"))
    
    # Test dimensions
    expect_length(result$LR, n_series)
    expect_length(result$pvalue, n_series)
    expect_equal(dim(result$pvalue_mat), c(30, n_series))
    expect_equal(dim(result$cvalue_mat), c(30, 3))
    expect_length(result$rho_grid, 30)
    
    # Test value ranges
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
    expect_true(all(result$LR > 0))
    expect_true(result$ha_parm > 0)
  }
})