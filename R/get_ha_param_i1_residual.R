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

  pow <- 1
  ctry <- get_cbar(0.95, distmat)
  while (pow > pow50) {
    c <- ctry
    sigdm_c <- get_sigma_residual(distmat, c, M)
    om_c <- t(R) %*% sigdm_c %*% R
    pow <- get_pow_qf(om_ho, om_c, e)
    ctry <- ctry / 2
  }
  c1 <- c

  pow <- 0
  ctry <- get_cbar(0.01, distmat)
  while (pow < pow50) {
    c <- ctry
    sigdm_c <- get_sigma_residual(distmat, c, M)
    om_c <- t(R) %*% sigdm_c %*% R
    pow <- get_pow_qf(om_ho, om_c, e)
    ctry <- 2 * ctry
  }
  c2 <- c

  ii <- 1
  while (abs(pow - pow50) > 0.01) {
    c <- (c1 + c2) / 2
    sigdm_c <- get_sigma_residual(distmat, c, M)
    om_c <- t(R) %*% sigdm_c %*% R
    pow <- get_pow_qf(om_ho, om_c, e)

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

#' @rdname get_ha_param_i1_residual
#' @keywords internal
get_ha_parm_I1_residual <- function(om_ho, distmat, R, e, M) {
  get_ha_param_i1_residual(
    om_ho = om_ho,
    distmat = distmat,
    R = R,
    e = e,
    M = M
  )
}
