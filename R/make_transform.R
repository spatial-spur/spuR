#' Create Transformation Matrix
#'
#' Internal constructor for SPUR spatial transformation matrices.
#'
#' @param s A numeric matrix of spatial coordinates.
#' @param transformation One of `"nn"`, `"iso"`, `"cluster"`, or `"lbmgls"`.
#' @param radius Numeric radius for the isotropic transformation.
#' @param cluster Numeric cluster ids for the cluster transformation.
#' @param latlong Logical; if `TRUE`, treat coordinates as latitude/longitude.
#' @return A numeric transformation matrix.
#' @keywords internal
make_transform <- function(s,
                           transformation,
                           radius = NULL,
                           cluster = NULL,
                           latlong = FALSE) {
  if (transformation == "nn") {
    return(nn_matrix(s = s, latlong = latlong))
  }
  if (transformation == "iso") {
    if (is.null(radius)) {
      stop("Radius required for iso transformation.")
    }
    return(iso_matrix(s = s, radius = radius, latlong = latlong))
  }
  if (transformation == "cluster") {
    if (is.null(cluster)) {
      stop("Cluster vector required for cluster transformation.")
    }
    return(cluster_matrix(cluster = cluster))
  }
  if (transformation == "lbmgls") {
    return(lbm_gls_matrix(s = s, latlong = latlong))
  }
  stop("Invalid transformation.")
}

#' Nearest-Neighbor Transformation Matrix
#'
#' Direct translation of `nn_matrix.mata`.
#'
#' @param s A numeric matrix of spatial coordinates.
#' @param latlong Logical; if `TRUE`, coordinates are latitude/longitude.
#' @return A numeric matrix.
#' @keywords internal
nn_matrix <- function(s, latlong = FALSE) {
  if (isTRUE(latlong)) {
    distmat <- get_distmat_lat_lon(s)
  } else {
    distmat <- get_distmat_euclidean(s)
  }
  distmat <- distmat / max(distmat)

  distmat_nozeros <- distmat + (distmat == 0) * 1e10
  rowmins <- apply(distmat_nozeros, 1L, min)
  rowmins_mat <- distmat_nozeros == matrix(rowmins, nrow = nrow(distmat_nozeros), ncol = ncol(distmat_nozeros))
  rowmins_mat <- rowmins_mat / rowSums(rowmins_mat)

  diag(nrow(s)) - rowmins_mat
}

#' Isotropic Transformation Matrix
#'
#' Direct translation of `iso_matrix.mata`.
#'
#' @param s A numeric matrix of spatial coordinates.
#' @param radius Numeric scalar radius (meters when `latlong = TRUE`).
#' @param latlong Logical; if `TRUE`, coordinates are latitude/longitude.
#' @return A numeric matrix.
#' @keywords internal
iso_matrix <- function(s, radius, latlong = FALSE) {
  distmat <- if (isTRUE(latlong)) {
    get_distmat_lat_lon(s)
  } else {
    get_distmat_euclidean(s)
  }

  if (isTRUE(latlong)) {
    distmat <- distmat * 3.14159265359 * 6371000.009 * 2
  }

  distmat_nozeros <- distmat + (distmat == 0) * 1e10
  dist_below_b <- distmat_nozeros <= radius
  dist_below_b <- dist_below_b / rowSums(dist_below_b)

  iso_mat <- diag(nrow(distmat)) - dist_below_b
  iso_mat[is.na(iso_mat)] <- 0
  iso_mat
}

#' Cluster Transformation Matrix
#'
#' Direct translation of `cluster_matrix.mata`.
#'
#' @param cluster A vector of cluster ids.
#' @return A numeric matrix.
#' @keywords internal
cluster_matrix <- function(cluster) {
  cluster <- as.vector(cluster)
  n <- length(cluster)
  clust_mat <- outer(cluster, cluster, "==")
  clust_mat <- clust_mat / rowSums(clust_mat)
  diag(n) - clust_mat
}

#' LBM-GLS Transformation Matrix
#'
#' Direct translation of `lbm_gls_matrix.mata`.
#'
#' @param s A numeric matrix of spatial coordinates.
#' @param latlong Logical; if `TRUE`, coordinates are latitude/longitude.
#' @return A numeric matrix.
#' @keywords internal
lbm_gls_matrix <- function(s, latlong = FALSE) {
  distmat <- if (isTRUE(latlong)) {
    get_distmat_lat_lon(s)
  } else {
    get_distmat_euclidean(s)
  }
  distmat <- distmat / max(distmat)

  sigma_lbm_dm <- get_sigma_lbm_dm(distmat)

  eig <- eigen(sigma_lbm_dm, symmetric = TRUE)
  values <- eig$values
  vectors <- eig$vectors

  idx <- order(values, decreasing = TRUE)
  values <- values[idx]
  vectors <- vectors[, idx, drop = FALSE]

  keep <- values > 1e-10
  values <- values[keep]
  vectors <- vectors[, keep, drop = FALSE]

  vectors %*% diag(1 / sqrt(values)) %*% t(vectors)
}
