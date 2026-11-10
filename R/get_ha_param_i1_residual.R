#' Compute Ha Parameter for I(1) Residual Test
#'
#' Translation of SPUR Mata `get_ha_parm_I1_residual.mata`.
#'
#' @param om_ho Null covariance matrix in low-frequency space.
#' @param distmat Distance matrix.
#' @param R Eigenvector matrix.
#' @param e Monte Carlo draw matrix.
#' @param M Residual-maker matrix.
#' @param R_M Optional basis after projecting out the regression design.
#' @return Scalar persistence parameter.
#' @keywords internal
get_ha_param_i1_residual <- function(om_ho, distmat, R, e, M = NULL, R_M = NULL) {
  if (nrow(e) < 2L) stop("I(1) tests require at least two low-frequency directions.")
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e)
  if (is.null(R_M)) R_M <- crossprod(M, R)
  basis <- R_M
  storage.mode(distmat) <- "double"
  storage.mode(basis) <- "double"
  cbars <- get_cbar(c(0.95, 0.01), distmat)
  evaluate <- function(c) pow_eval(.Call(cpp_om_grid, distmat, basis, as.numeric(c)))
  .calibrate_i1_power(evaluate, cbars[1L], cbars[2L])
}
