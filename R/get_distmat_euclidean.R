#' Calculate Euclidean Distance Matrix
#'
#' Computes the matrix of Euclidean distances between all pairs of locations
#' in the input matrix.
#'
#' @param s A numeric matrix where rows represent locations and columns represent dimensions.
#' @return A numeric matrix of pairwise Euclidean distances between locations.
get_distmat_euclidean <- function(s) {
  # Use R's built-in dist function and convert to matrix
  distmat <- as.matrix(stats::dist(s, method = "euclidean"))
  
  return(distmat)
}

#' @rdname get_distmat_euclidean
getdistmat_euclidian <- function(s) {
  get_distmat_euclidean(s = s)
}
