#' Create Transformation Matrix
#'
#' Creates a transformation matrix based on the specified transformation method.
#'
#' @param s A numeric matrix of spatial coordinates.
#' @param transformation A string specifying the transformation type: "nn", "iso", "cluster", or "lbmgls".
#' @param radius A numeric scalar, required for "iso" transformation.
#' @param cluster A vector of cluster IDs, required for "cluster" transformation.
#' @param latlong Logical; if TRUE, coordinates are latitude/longitude.
#' @return A numeric matrix representing the spatial transformation.
#' @export
make_transform <- function(s, transformation, radius = NULL, cluster = NULL, latlong = FALSE) {
  H <- NULL
  
  if (transformation == "nn") {
    H <- nn_matrix(s, latlong)
  } else if (transformation == "iso") {
    if (is.null(radius)) {
      stop("Radius required for iso transformation")
    }
    H <- iso_matrix(s, radius, latlong)
  } else if (transformation == "cluster") {
    if (is.null(cluster)) {
      stop("Cluster vector required for cluster transformation")
    }
    H <- cluster_matrix(cluster)
  } else if (transformation == "lbmgls") {
    H <- lbm_gls_matrix(s, latlong)
  } else {
    stop("Invalid transformation.")
  }
  
  return(H)
}

# Helper transformation functions that need to be implemented:
# nn_matrix, iso_matrix, cluster_matrix, lbm_gls_matrix

#' Create Nearest Neighbor Transformation Matrix
#'
#' @param s A numeric matrix of spatial coordinates.
#' @return A numeric matrix representing nearest neighbor transformation.
#' @keywords internal
nn_matrix <- function(s) {
  # Implement nearest neighbor transformation
  stop("nn_matrix not yet implemented")
}

#' Create Isotropic Transformation Matrix
#'
#' @param s A numeric matrix of spatial coordinates.
#' @param radius A numeric scalar defining the neighborhood radius.
#' @return A numeric matrix representing isotropic transformation.
#' @keywords internal
iso_matrix <- function(s, radius, latlong = FALSE) {
  # Calculate distance matrix
  if (ncol(s) == 2 && latlong) {
    distmat <- getdistmat_lat_lon(s)
  } else {
    distmat <- getdistmat_euclidian(s)
  }
  
  # Create isotropic weights
  n <- nrow(s)
  H <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    # Find neighbors within radius
    neighbors <- which(distmat[i, ] <= radius)
    if (length(neighbors) > 0) {
      # Assign equal weights to all neighbors
      H[i, neighbors] <- 1 / length(neighbors)
    }
  }
  
  return(H)
}

#' Create Cluster Transformation Matrix
#'
#' @param cluster A vector of cluster IDs.
#' @return A numeric matrix representing cluster transformation.
#' @keywords internal
cluster_matrix <- function(cluster) {
  n <- length(cluster)
  H <- matrix(0, nrow = n, ncol = n)
  
  # Get unique cluster IDs
  unique_clusters <- unique(cluster)
  
  # Create cluster weights
  for (c in unique_clusters) {
    members <- which(cluster == c)
    if (length(members) > 0) {
      for (i in members) {
        H[i, members] <- 1 / length(members)
      }
    }
  }
  
  return(H)
}

#' Create LBM-GLS Transformation Matrix
#'
#' @param s A numeric matrix of spatial coordinates.
#' @return A numeric matrix representing LBM-GLS transformation.
#' @keywords internal
lbm_gls_matrix <- function(s, latlong = FALSE) {
  # Calculate appropriate distance matrix
  if (ncol(s) == 2 && latlong) {
    distmat <- getdistmat_lat_lon(s) # slightly different from stata 
  } else {
    distmat <- getdistmat_euclidian(s)
  }
  
  # Normalize distance matrix
  distmat <- distmat / max(distmat)
  
  # Get LBM covariance matrix
  sigma_lbm_dm = get_sigma_lbm_dm(distmat) # exactly identical
  
  # Compute eigendecomposition
  eig <- eigen(sigma_lbm_dm, symmetric = TRUE) 
  values <- eig$values
  vectors <- eig$vectors
  
  # Sort eigenvalues decreasing (to match Stata's order() function)
  idx <- order(values, decreasing = TRUE)
  values <- values[idx]
  vectors <- vectors[, idx]
  
  # Filter out small or negative eigenvalues
  idx <- values > 1e-10
  values <- values[idx]
  vectors <- vectors[, idx]
  
  # Create GLS transformation matrix
  dsi <- 1 / sqrt(values)
  H <- vectors %*% diag(dsi) %*% t(vectors)
  
  return(H)
}