#' Compute Ha Parameter for I(1) Residual Test
#'
#' Translation of SPUR Mata `get_ha_parm_I1_residual.mata`.
#'
#' @param om_ho Null covariance matrix in low-frequency space.
#' @param distmat Distance matrix.
#' @param R Eigenvector matrix.
#' @param e Monte Carlo draw matrix.
#' @param M Residual-maker matrix.
#' @return Scalar persistence parameter.
#' @keywords internal
get_ha_param_i1_residual <- function(om_ho, distmat, R, e, M = NULL, R_M = NULL) {
  pow50 <- 0.5
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e)
  if (is.null(R_M)) {
    R_M <- crossprod(M, R)
  }
  storage.mode(distmat) <- "double"
  storage.mode(R_M) <- "double"
  cbars <- get_cbar(c(0.95, 0.01), distmat)

  pow <- 1
  ctry <- cbars[1L]
  while (pow > pow50) {
    c <- ctry
    om_c <- .Call(cpp_om_grid, distmat, R_M, as.numeric(c))
    pow <- pow_eval(om_c)
    ctry <- ctry / 2
  }
  c1 <- c

  pow <- 0
  ctry <- cbars[2L]
  while (pow < pow50) {
    c <- ctry
    om_c <- .Call(cpp_om_grid, distmat, R_M, as.numeric(c))
    pow <- pow_eval(om_c)
    ctry <- 2 * ctry
  }
  c2 <- c

  ii <- 1
  while (abs(pow - pow50) > 0.01) {
    c <- (c1 + c2) / 2
    om_c <- .Call(cpp_om_grid, distmat, R_M, as.numeric(c))
    pow <- pow_eval(om_c)

    if (pow > pow50) {
      c2 <- c
    } else if (pow < pow50) {
      c1 <- c
    }

    ii <- ii + 1
    if (ii > 20) {
      break
    }
  }

  c
}
