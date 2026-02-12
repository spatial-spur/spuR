#' Compute Ha Parameter for I1 Test
#'
#' Determines the parameter \code{g} such that the power of the test is approximately 50%.
#' The power is computed by the function \code{get_pow_qf()}, and \code{g} is determined
#' through an iterative bisection procedure.
#'
#' @param om_ho A numeric matrix.
#' @param distmat A numeric distance matrix.
#' @param R A numeric matrix of eigenvectors.
#' @param e A numeric matrix (or vector) used in the computation of power.
#' @return A numeric scalar representing the parameter \code{g}.
get_ha_param_i1 <- function(om_ho, distmat, R, e) {
  ## Goal: choose c so that the quadratic‑form test has ≈ 50 % power

  ## ---- lower bracket: power just under 0.5 ------------------
  pow  <- 1
  ctry <- get_cbar(0.95, distmat)
  while (pow > 0.5) {
    c        <- ctry
    sigdm_c  <- get_sigma_dm(distmat, c)
    om_c     <- t(R) %*% sigdm_c %*% R
    pow      <- get_pow_qf(om_ho, om_c, e)
    ctry     <- ctry / 2
  }
  c1 <- c        # store lower bound

  ## ---- upper bracket: power just over 0.5 -------------------
  pow  <- 0
  ctry <- get_cbar(0.01, distmat)
  while (pow < 0.5) {
    c        <- ctry
    sigdm_c  <- get_sigma_dm(distmat, c)
    om_c     <- t(R) %*% sigdm_c %*% R
    pow      <- get_pow_qf(om_ho, om_c, e)
    ctry     <- 2 * ctry
  }
  c2 <- c        # store upper bound

  ## ---- binary search for 50 % power -------------------------
  ii  <- 1
  while (abs(pow - 0.5) > 0.01) {
    c        <- (c1 + c2) / 2
    sigdm_c  <- get_sigma_dm(distmat, c)
    om_c     <- t(R) %*% sigdm_c %*% R
    pow      <- get_pow_qf(om_ho, om_c, e)

    if (pow > 0.5) {
      c2 <- c
    } else {          # pow < 0.5
      c1 <- c
    }
    ii <- ii + 1
    if (ii > 20) break
  }

  ## return parameter that achieves ≈50 % power
  c
}
#-------------------------------------------------------------

#' @rdname get_ha_param_i1
get_ha_parm_I1 <- function(om_ho, distmat, R, e) {
  get_ha_param_i1(om_ho = om_ho, distmat = distmat, R = R, e = e)
}
