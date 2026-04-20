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
get_ha_param_i1_residual <- function(om_ho, distmat, R, e, M) {
  pow50 <- 0.5
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e)

  pow <- 1
  ctry <- get_cbar(0.95, distmat)
  while (pow > pow50) {
    c <- ctry
    sigdm_c <- get_sigma_residual(distmat, c, M)
    om_c <- t(R) %*% sigdm_c %*% R
    pow <- pow_eval(om_c)
    ctry <- ctry / 2
  }
  c1 <- c

  pow <- 0
  ctry <- get_cbar(0.01, distmat)
  while (pow < pow50) {
    c <- ctry
    sigdm_c <- get_sigma_residual(distmat, c, M)
    om_c <- t(R) %*% sigdm_c %*% R
    pow <- pow_eval(om_c)
    ctry <- 2 * ctry
  }
  c2 <- c

  ii <- 1
  while (abs(pow - pow50) > 0.01) {
    c <- (c1 + c2) / 2
    sigdm_c <- get_sigma_residual(distmat, c, M)
    om_c <- t(R) %*% sigdm_c %*% R
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
