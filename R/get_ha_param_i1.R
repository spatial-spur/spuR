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
get_ha_param_i1 <- function(om_ho, distmat, R, e) {
  if (nrow(e) < 2L) stop("I(1) tests require at least two low-frequency directions.")
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e)
  basis <- sweep(R, 2L, colMeans(R), "-")
  storage.mode(distmat) <- "double"
  storage.mode(basis) <- "double"
  cbars <- get_cbar(c(0.95, 0.01), distmat)
  evaluate <- function(c) pow_eval(.Call(cpp_om_grid, distmat, basis, as.numeric(c)))
  .calibrate_i1_power(evaluate, cbars[1L], cbars[2L])
}
# Calibrate against simulated power, with finite brackets and an explicit failure
# when the selected low-frequency experiment cannot achieve the target power.
.calibrate_i1_power <- function(evaluate, lower, upper) {
  bracket <- function(c, direction) {
    for (attempt in seq_len(40L)) {
      power <- evaluate(c)
      if (!is.finite(power)) stop("I(1) power calibration produced a non-finite probability.")
      if ((direction < 1 && power <= 0.5) || (direction > 1 && power >= 0.5)) {
        return(list(c = c, power = power))
      }
      c <- c * direction
    }
    stop("Cannot bracket 50% power for the I(1) alternative; increase q or nrep.")
  }
  lo <- bracket(lower, 0.5)
  hi <- bracket(upper, 2)
  best <- if (abs(lo$power - 0.5) < abs(hi$power - 0.5)) lo else hi
  for (iteration in seq_len(60L)) {
    c <- (lo$c + hi$c) / 2
    power <- evaluate(c)
    if (!is.finite(power)) stop("I(1) power calibration produced a non-finite probability.")
    if (abs(power - 0.5) < abs(best$power - 0.5)) best <- list(c = c, power = power)
    if (abs(power - 0.5) <= 0.01) return(c)
    if (power > 0.5) hi <- list(c = c, power = power) else lo <- list(c = c, power = power)
    if (hi$c - lo$c <= 1e-10 * hi$c) break
  }
  best$c
}
