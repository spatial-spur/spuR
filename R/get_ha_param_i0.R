#' Compute Ha Parameter for I0 Test
#'
#' Determines the parameter \code{g} such that the power of the test is approximately 50%.
#' The power is computed by the function \code{get_pow_qf()}, and \code{g} is determined
#' through an iterative bisection procedure.
#'
#' @param om_ho A numeric matrix.
#' @param om_i0 A numeric matrix.
#' @param om_bm A numeric matrix.
#' @param e A numeric matrix (or vector) used in the computation of power.
#' @return A numeric scalar representing the parameter \code{g}.
#' @keywords internal
get_ha_param_i0 <- function(om_ho, om_i0, om_bm, e) {
  lower <- t(chol(om_ho))
  no_shape_information <- function(om) {
    left <- forwardsolve(lower, om)
    whitened <- t(forwardsolve(lower, t(left)))
    scale <- mean(diag(whitened))
    is.finite(scale) && scale > 0 &&
      max(abs(whitened / scale - diag(nrow(whitened)))) <= 1e-10
  }
  if (no_shape_information(om_i0) && no_shape_information(om_bm)) {
    stop("Cannot bracket 50% power: the I(0) alternative has no covariance-shape information.")
  }
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e)
  evaluate <- function(g) {
    power <- pow_eval(om_i0 + g * om_bm)
    if (!is.finite(power)) stop("I(0) power calibration produced a non-finite probability.")
    power
  }
  unbracketed <- function() {
    stop("Cannot bracket 50% power for the I(0) alternative; increase q or nrep.")
  }

  # Finite brackets also cover covariance experiments whose limiting power
  # cannot reach 50%, even when q is greater than one.
  gtry <- 1
  for (attempt in seq_len(64L)) {
    g <- gtry
    pow <- evaluate(g)
    if (pow <= 0.5) break
    if (attempt == 64L) unbracketed()
    gtry <- g / 2
  }
  g1 <- g

  gtry <- 30
  for (attempt in seq_len(64L)) {
    g <- gtry
    pow <- evaluate(g)
    if (pow >= 0.5) break
    if (attempt == 64L) unbracketed()
    gtry <- g * 2
  }
  g2 <- g

  # Preserve the calibrated alternative on ordinary inputs, and never return
  # an alternative that failed to attain the requested empirical power.
  for (iteration in seq_len(60L)) {
    if (abs(pow - 0.5) <= 0.01) return(g)
    g <- (g1 + g2) / 2
    pow <- evaluate(g)
    if (pow > 0.5) {
      g2 <- g
    } else if (pow < 0.5) {
      g1 <- g
    }
  }
  if (abs(pow - 0.5) <= 0.01) return(g)
  stop("Cannot calibrate the I(0) alternative to 50% power; increase q or nrep.")
}
