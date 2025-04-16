# Test script to verify lvech implementation matches MATLAB

# R implementation of lvech
lvech <- function(S) {
  return(S[lower.tri(S, diag = FALSE)])
}

# Create a simple test matrix (same structure in both languages)
test_matrix <- matrix(1:16, nrow=4)
print("Test matrix:")
print(test_matrix)

# Apply lvech function in R
r_result <- lvech(test_matrix)
print("R lvech result:")
print(r_result)
