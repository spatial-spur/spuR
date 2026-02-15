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
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e)

  # First phase: find g1 such that power is below 50%
  pow <- 1
  gtry <- 1
  while (pow > 0.5) {
    g <- gtry
    pow <- pow_eval(om_i0 + g * om_bm)
    gtry <- g / 2
  }
  g1 <- g
  
  # Second phase: find g2 such that power is above 50%
  pow <- 0
  gtry <- 30
  while (pow < 0.5) {
    g <- gtry
    pow <- pow_eval(om_i0 + g * om_bm)
    gtry <- g * 2
  }
  g2 <- g
  
  # Bisection: refine g so that power is approximately 50%
  ii <- 1
  while (abs(pow - 0.5) > 0.01) {
    g <- (g1 + g2) / 2
    pow <- pow_eval(om_i0 + g * om_bm)
    if (pow > 0.5) {
      g2 <- g
    } else if (pow < 0.5) {
      g1 <- g
    }
    ii <- ii + 1
    if (ii > 20) break
  }
  
  return(g)
}
