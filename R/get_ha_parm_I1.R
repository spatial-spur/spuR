#' Compute Ha Parameter for I1 Test
#'
#' Determines the parameter \code{g} such that the power of the test is approximately 50%.
#' The power is computed by the function \code{getpow_qf}, and \code{g} is determined
#' through an iterative bisection procedure.
#'
#' @param om_ho A numeric matrix.
#' @param om_i0 A numeric matrix.
#' @param om_bm A numeric matrix.
#' @param e A numeric matrix (or vector) used in the computation of power.
#' @return A numeric scalar representing the parameter \code{g}.
#' @export
get_ha_parm_I1 <- function(om_ho, distmat, R, e) {
  ## Goal: choose c so that the quadratic‑form test has ≈ 50 % power

  ## ---- lower bracket: power just under 0.5 ------------------
  pow  <- 1
  ctry <- getcbar(0.95, distmat)
  while (pow > 0.5) {
    c        <- ctry
    sigdm_c  <- get_sigma_dm(distmat, c)
    om_c     <- t(R) %*% sigdm_c %*% R
    pow      <- getpow_qf(om_ho, om_c, e)
    ctry     <- ctry / 2
  }
  c1 <- c        # store lower bound

  ## ---- upper bracket: power just over 0.5 -------------------
  pow  <- 0
  ctry <- getcbar(0.01, distmat)
  while (pow < 0.5) {
    c        <- ctry
    sigdm_c  <- get_sigma_dm(distmat, c)
    om_c     <- t(R) %*% sigdm_c %*% R
    pow      <- getpow_qf(om_ho, om_c, e)
    ctry     <- 2 * ctry
  }
  c2 <- c        # store upper bound

  ## ---- binary search for 50 % power -------------------------
  ii  <- 1
  repeat {
    if (abs(pow - 0.5) <= 0.01 || ii > 20) break
    c        <- (c1 + c2) / 2
    sigdm_c  <- get_sigma_dm(distmat, c)
    om_c     <- t(R) %*% sigdm_c %*% R
    pow      <- getpow_qf(om_ho, om_c, e)

    if (pow > 0.5) {
      c2 <- c
    } else {          # pow < 0.5
      c1 <- c
    }
    ii <- ii + 1
  }

  ## return parameter that achieves ≈50 % power
  c
}
#-------------------------------------------------------------
