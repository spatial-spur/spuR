#' Set and get thread count for spuR
#'
#' Configure the number of threads used for parallel computation in \pkg{spuR}.
#'
#' @param num_threads Integer >= 1 specifying the number of threads, or \code{"auto"}
#'   / \code{NULL} to restore automatic thread detection.
#'
#' @return \code{spur_set_num_threads} invisibly returns the previous thread count setting.
#'   \code{spur_get_num_threads} returns the currently configured or active thread count.
#'
#' @details
#' \pkg{spuR} supports multithreading via \pkg{RcppParallel} for grid evaluations,
#' Monte Carlo persistence intervals, and distance calculations.
#'
#' Thread resolution follows this priority:
#' \enumerate{
#'   \item Explicit argument \code{nthreads} in \code{\link{spur}} and test functions
#'   \item Package option \code{options(spuR.threads = N)} (set directly or via \code{spur_set_num_threads(N)})
#'   \item Environment variable \code{SPUR_NUM_THREADS}
#'   \item CRAN test limit \code{_R_CHECK_LIMIT_CORES_} (capped at 2 if set)
#'   \item General R parallel option \code{options(mc.cores = N)}
#'   \item \pkg{RcppParallel} environment variable \code{RCPP_PARALLEL_NUM_THREADS}
#'   \item Automatic hardware detection (all available logical cores)
#' }
#'
#' Setting \code{num_threads = 1} (or \code{options(spuR.threads = 1)}) completely bypasses
#' background task submission, ensuring pure single-threaded execution. This is recommended
#' when parallelising outer loops with \code{parallel::mclapply}, \pkg{future}, or when running
#' on shared cluster nodes.
#'
#' @examples
#' \dontrun{
#' # Check currently active threads
#' spur_get_num_threads()
#'
#' # Force single-threaded execution (e.g. inside mclapply or cluster jobs)
#' old_th <- spur_set_num_threads(1)
#'
#' # Restore previous setting
#' spur_set_num_threads(old_th)
#'
#' # Or restore default automatic threading
#' spur_set_num_threads("auto")
#' }
#'
#' @export
spur_set_num_threads <- function(num_threads = "auto") {
  old <- spur_get_num_threads()
  if (is.null(num_threads) || identical(num_threads, "auto")) {
    options(spuR.threads = NULL)
  } else {
    if (!is.numeric(num_threads) || length(num_threads) != 1L ||
        !is.finite(num_threads) || num_threads < 1 ||
        as.integer(num_threads) != num_threads) {
      stop("`num_threads` must be a single positive integer, NULL, or \"auto\".")
    }
    options(spuR.threads = as.integer(num_threads))
  }
  invisible(old)
}

#' @rdname spur_set_num_threads
#' @export
spur_get_num_threads <- function() {
  .get_spur_threads()
}

.get_spur_threads <- function(nthreads = NULL) {
  if (!is.null(nthreads)) {
    if (!is.numeric(nthreads) || length(nthreads) != 1L || !is.finite(nthreads) ||
        nthreads < 1 || as.integer(nthreads) != nthreads) {
      stop("`nthreads` must be a single positive integer or NULL.")
    }
    return(as.integer(nthreads))
  }
  opt <- getOption("spuR.threads")
  if (!is.null(opt) && is.numeric(opt) && length(opt) == 1L && is.finite(opt) && opt >= 1) {
    return(as.integer(opt))
  }
  env <- Sys.getenv("SPUR_NUM_THREADS", "")
  if (nzchar(env)) {
    val <- suppressWarnings(as.integer(env))
    if (!is.na(val) && val >= 1L) return(val)
  }
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    return(2L)
  }
  mc <- getOption("mc.cores")
  if (!is.null(mc) && is.numeric(mc) && length(mc) == 1L && is.finite(mc) && mc >= 1) {
    return(as.integer(mc))
  }
  env_rcpp <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")
  if (nzchar(env_rcpp)) {
    val <- suppressWarnings(as.integer(env_rcpp))
    if (!is.na(val) && val >= 1L) return(val)
  }
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    return(as.integer(RcppParallel::defaultNumThreads()))
  }
  1L
}
