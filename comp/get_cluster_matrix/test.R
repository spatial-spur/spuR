# Source the function
source("R/get_cluster_matrix.R")

# Create test data with clustering variable
cities <- data.frame(
  city = c("New York", "Los Angeles", "Chicago", "Houston", "Phoenix", "Philadelphia"),
  cluster = c(1, 2, 1, 2, 2, 1)  # East=1, West=2 clustering
)

# Test case 1: All rows
cat("==== Test 1: All rows ====\n")
result1 <- get_cluster_matrix(cities$cluster)
print(result1)
print(class(result1))
print(dim(result1))

