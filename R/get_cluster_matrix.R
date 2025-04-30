#' Get Cluster Matrix
#'
#' Extracts the cluster variable values for specified observations.
#' This is a direct R translation of the Stata function.
#'
#' @param cluster_values A vector containing cluster identifiers
#' @param use_rows Logical vector indicating which rows to include
#' @return A matrix (column vector) of cluster values for selected observations
#' @export
get_cluster_matrix <- function(cluster_values, use_rows = NULL) {
  # If use_rows is provided, subset the cluster values
  if (!is.null(use_rows)) {
    cluster_values <- cluster_values[use_rows]
  }
  
  # Convert to a column matrix to match Stata's behavior
  # (Stata st_data returns a matrix even for a single column)
  cluster_matrix <- as.matrix(cluster_values)
  
  return(cluster_matrix)
}