#' Compute Demeaned Sigma LBm
#'
#' Computes the sigma_lbm matrix using the provided distance matrix and then demeans it.
#'
#' @param distmat A numeric matrix of distances.
#' @return A numeric matrix representing the demeaned sigma_lbm.
get_sigma_lbm_dm <- function(distmat) {
  sigma_lbm <- get_sigma_lbm(distmat)
  sigma_lbm_dm <- demean_sigma(sigma_lbm)
  return(sigma_lbm_dm)
}
