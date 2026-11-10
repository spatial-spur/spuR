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
  storage.mode(distmat) <- "double"
  .Call(cpp_nn_matrix, distmat)
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
  storage.mode(distmat) <- "double"
  .Call(cpp_iso_matrix, distmat, as.numeric(radius))
}

#' Cluster Transformation Matrix
#'
#' Direct translation of `cluster_matrix.mata`.
#'
#' @param cluster A vector of cluster ids.
#' @return A numeric matrix.
#' @keywords internal
cluster_matrix <- function(cluster) {
  .Call(cpp_cluster_matrix, as.numeric(cluster))
}

#' LBM-GLS Transformation Matrix
#'
#' Direct translation of `lbm_gls_matrix.mata`.
#'
#' @param s A numeric matrix of spatial coordinates.
#' @param latlong Logical; if `TRUE`, coordinates are latitude/longitude.
#' @return A numeric matrix.
#' @keywords internal
lbm_gls_matrix <- function(s, latlong = FALSE, distmat = NULL, sigma_lbm_dm = NULL) {
  if (is.null(sigma_lbm_dm)) {
    if (is.null(distmat)) {
      distmat <- if (isTRUE(latlong)) {
        get_distmat_lat_lon(s)
      } else {
        get_distmat_euclidean(s)
      }
      distmat <- distmat / max(distmat)
    }
    sigma_lbm_dm <- get_sigma_lbm_dm(distmat)
  }
  storage.mode(sigma_lbm_dm) <- "double"
  .Call(cpp_lbm_gls_from_sigma, sigma_lbm_dm)
}
