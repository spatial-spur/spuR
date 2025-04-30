# Source the transform function
source("R/transform.R")  # Note: filename typo "tramsform.R" vs "transform.R"

# Create a small test dataset
test_data <- data.frame(
  x = c(1.5, 2.3, 3.7, 4.2, 5.8),
  y = c(0.5, 1.2, 1.8, 2.5, 3.1)
)

# Create three different test matrices
# 1. Identity matrix (no transformation)
H_identity <- diag(5)

# 2. Simple averaging matrix (each row is averaged with neighbors)
H_avg <- matrix(0, 5, 5)
for(i in 1:5) {
  if(i == 1) {
    H_avg[i, 1:2] <- c(0.7, 0.3)
  } else if(i == 5) {
    H_avg[i, 4:5] <- c(0.3, 0.7)
  } else {
    H_avg[i, (i-1):(i+1)] <- c(0.2, 0.6, 0.2)
  }
}

# 3. Symmetric matrix (like what might come from LBM-GLS)
H_sym <- matrix(c(
  2.1, -0.3, -0.5, -0.2, -0.1,
  -0.3, 1.9, -0.4, -0.1, -0.1,
  -0.5, -0.4, 2.5, -0.6, -0.3,
  -0.2, -0.1, -0.6, 1.8, -0.5,
  -0.1, -0.1, -0.3, -0.5, 2.0
), nrow=5, byrow=TRUE)

# Set precision for output
options(digits = 10)

# Test 1: Identity matrix, no demeaning
cat("Test 1: Identity matrix, 'x' variable, no demeaning\n")
result1 <- transform(test_data, "x", H_identity, "iso")
cat("Result:", result1, "\n\n")

# Test 2: Averaging matrix, no demeaning
cat("Test 2: Averaging matrix, 'x' variable, no demeaning\n")
result2 <- transform(test_data, "x", H_avg, "iso")
cat("Result:", result2, "\n\n")

# Test 3: Symmetric matrix with demeaning (lbmgls)
cat("Test 3: Symmetric matrix, 'y' variable, with demeaning (lbmgls)\n")
result3 <- transform(test_data, "y", H_sym, "lbmgls")
cat("Result:", result3, "\n")
cat("Mean of result3:", mean(result3), "\n\n")

# Test 4: With subsetting
use_rows <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
cat("Test 4: Identity matrix, 'x' variable, with subsetting\n")
result4 <- transform(test_data, "x", H_identity[c(1,3,4), c(1,3,4)], "iso", use_rows)
cat("Result:", result4, "\n\n")
