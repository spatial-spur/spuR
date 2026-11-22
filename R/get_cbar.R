#' Compute cbar Value
#'
#' Determines the value of \code{cbar} such that \eqn{mean(exp(-cbar \times vd))} is approximately equal to
#' a given \code{rhobar}, where \code{vd} is the vectorized lower-triangular portion (excluding the diagonal)
#' of the distance matrix \code{distmat}. This function uses iterative adjustments and bisection to find
#' an appropriate value.
#'
#' @param rhobar A numeric scalar representing the target value.
#' @param distmat A numeric matrix of distances.
#' @return A numeric scalar \code{cbar} satisfying the desired property.
#' @keywords internal
get_cbar <- function(rhobar, distmat) {
  .make_spur_decay_solver(distmat)(rhobar)
}

.make_spur_decay_solver <- function(distmat) {
  distances <- lvech(distmat)
  if (!length(distances) || any(!is.finite(distances) | distances < 0)) {
    stop("SPUR decay calibration requires finite nonnegative pairwise distances.")
  }
  correlation_floor <- mean(distances == 0)
  averages <- new.env(parent = emptyenv())
  roots <- new.env(parent = emptyenv())
  average <- function(c) {
    key <- sprintf("%a", c)
    out <- averages[[key]]
    if (is.null(out)) {
      out <- mean(exp(-c * distances))
      averages[[key]] <- out
    }
    out
  }

  function(rho) {
    if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||
        rho <= correlation_floor || rho >= 1) {
      stop("SPUR target correlation must lie above the co-located-pair fraction and below one.")
    }
    key <- sprintf("%a", rho)
    out <- roots[[key]]
    if (is.null(out)) {
      out <- .solve_spur_decay(rho, average)
      roots[[key]] <- out
    }
    out
  }
}

.solve_spur_decay <- function(rhobar, average) {
  c0 <- 10
  c1 <- 10

  while (average(c0) <= rhobar) {
    c1 <- c0
    c0 <- c0 / 2
    if (c0 == 0) {
      stop("SPUR decay calibration could not represent a lower bracket.")
    }
  }

  while (average(c1) > rhobar) {
    c0 <- c1
    c1 <- min(2 * c1, .Machine$double.xmax)
    if (c1 == c0) {
      stop("SPUR decay calibration could not represent an upper bracket.")
    }
  }

  # Keep the replication tolerance for an ordinary bracket. An extended
  # bracket uses relative width and a midpoint that cannot overflow.
  extended <- c1 > 10240
  tolerance <- if (extended) 1e-12 * c1 else 0.001
  midpoint <- function() if (extended) c0 + 0.5 * (c1 - c0) else sqrt(c0 * c1)
  while (c1 - c0 > tolerance) {
    cm <- midpoint()
    if (cm <= c0 || cm >= c1) break
    if (average(cm) < rhobar) {
      c1 <- cm
    } else {
      c0 <- cm
    }
  }

  midpoint()
}
