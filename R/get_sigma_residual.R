#' Residual Covariance Matrix
#'
#' Translation of SPUR Mata `get_sigma_residual.mata`.
#'
#' @param distmat A distance matrix.
#' @param c Spatial persistence parameter.
#' @param M Residual-maker matrix.
#' @return Residual covariance matrix.
#' @keywords internal
get_sigma_residual <- function(distmat, c, M) {
  sigma <- exp(-c * distmat)
  M %*% sigma %*% t(M)
}
