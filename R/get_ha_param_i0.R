#' Compute Ha Parameter for I0 Test
#'
#' Determines the parameter \code{g} such that the power of the test is approximately 50%.
#' An iterative bisection procedure computes \code{g} using the power of a quadratic-form test.
#'
#' @param om_ho A numeric matrix.
#' @param om_i0 A numeric matrix.
#' @param om_bm A numeric matrix.
#' @param e A numeric matrix (or vector) used in the computation of power.
#' @return A numeric scalar representing the parameter \code{g}.
#' @keywords internal
get_ha_param_i0 <- function(om_ho, om_i0, om_bm, e) {
  .fit_spur_i0(om_ho, om_i0, om_bm, e)$rate
}

.fit_spur_i0 <- function(om_ho, om_i0, om_bm, e) {
  pow_eval <- .make_pow_qf_evaluator(om0 = om_ho, e = e)

  # First phase: find g1 such that power is below 50%
  pow <- 1
  gtry <- 1
  bracket_steps <- 0L
  while (pow > 0.5) {
    bracket_steps <- bracket_steps + 1L
    if (bracket_steps > 200L) {
      stop("get_ha_param_i0 step-1 bracketing did not converge.")
    }
    g <- gtry
    pow <- pow_eval(om_i0 + g * om_bm)
    gtry <- g / 2
  }
  if (pow == 0.5) return(list(rate = g, power = NA_real_))
  g1 <- g
  
  # Second phase: find g2 such that power is above 50%
  best_g <- g
  best_power <- pow
  limit_shape <- om_bm / sum(diag(om_bm))
  pow <- 0
  gtry <- if (bracket_steps > 1L) 2 * g1 else 30
  bracket_steps <- 0L
  while (pow < 0.5) {
    bracket_steps <- bracket_steps + 1L
    if (bracket_steps > 200L) {
      return(list(rate = best_g, power = best_power))
    }
    g <- gtry
    mixture <- om_i0 + g * om_bm
    pow <- pow_eval(mixture)
    if (pow > best_power) {
      best_power <- pow
      best_g <- g
    }
    if (pow < 0.5 && max(abs(mixture / sum(diag(mixture)) - limit_shape)) <= 8 * .Machine$double.eps) {
      return(list(rate = best_g, power = best_power))
    }
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
  
  list(rate = g, power = NA_real_)
}
