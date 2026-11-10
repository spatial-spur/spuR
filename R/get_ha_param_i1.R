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
#' @keywords internal
get_ha_param_i1 <- function(om_ho, distmat, R, e, nthreads = NULL) {
  nthreads <- .get_spur_threads(nthreads)
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e, nthreads = nthreads)
  R_dm <- sweep(R, 2L, colMeans(R), "-")
  storage.mode(distmat) <- "double"
  storage.mode(R_dm) <- "double"
  cbars <- get_cbar(c(0.95, 0.01), distmat)

  ## Goal: choose c so that the quadratic‑form test has ≈ 50 % power

  ## ---- lower bracket: power just under 0.5 ------------------
  pow  <- 1
  ctry <- cbars[1L]
  while (pow > 0.5) {
    c        <- ctry
    om_c     <- rcpp_om_grid(distmat, R_dm, as.numeric(c), as.integer(nthreads))
    pow      <- pow_eval(om_c)
    ctry     <- ctry / 2
  }
  c1 <- c        # store lower bound

  ## ---- upper bracket: power just over 0.5 -------------------
  pow  <- 0
  ctry <- cbars[2L]
  while (pow < 0.5) {
    c        <- ctry
    om_c     <- rcpp_om_grid(distmat, R_dm, as.numeric(c), as.integer(nthreads))
    pow      <- pow_eval(om_c)
    ctry     <- 2 * ctry
  }
  c2 <- c        # store upper bound

  ## ---- binary search for 50 % power -------------------------
  ii  <- 1
  while (abs(pow - 0.5) > 0.01) {
    c        <- (c1 + c2) / 2
    om_c     <- rcpp_om_grid(distmat, R_dm, as.numeric(c), as.integer(nthreads))
    pow      <- pow_eval(om_c)

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
