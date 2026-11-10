test_that("thread control helpers set and retrieve thread configuration correctly in spuR", {
  old_th <- spur_get_num_threads()
  on.exit({
    spur_set_num_threads("auto")
    options(spuR.threads = NULL)
    options(mc.cores = NULL)
    Sys.unsetenv("SPUR_NUM_THREADS")
    Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  }, add = TRUE)

  # 1. Setting specific thread counts via helper
  spur_set_num_threads(2)
  expect_equal(spur_get_num_threads(), 2L)
  expect_equal(getOption("spuR.threads"), 2L)

  # 2. Restoring to auto
  spur_set_num_threads("auto")
  expect_null(getOption("spuR.threads"))

  # 3. Direct options() control
  options(spuR.threads = 1L)
  expect_equal(spur_get_num_threads(), 1L)

  # 4. mc.cores fallback
  options(spuR.threads = NULL)
  options(mc.cores = 3L)
  expect_equal(spur_get_num_threads(), 3L)

  # 5. CRAN _R_CHECK_LIMIT_CORES_ policy cap at 2
  options(spuR.threads = NULL)
  options(mc.cores = NULL)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_equal(spur_get_num_threads(), 2L)
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")

  # 6. Invalid inputs throw error
  expect_error(spur_set_num_threads(0))
  expect_error(spur_set_num_threads(-2))
  expect_error(spur_set_num_threads("many"))
})

test_that("spur(..., nthreads = 1) parameter and single-threaded execution match multi-threaded execution", {
  on.exit(spur_set_num_threads("auto"), add = TRUE)

  data(spur_example, package = "spuR")

  # Run multi-threaded via auto
  spur_set_num_threads("auto")
  res_multi <- spurtest_i0(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42, q = 6, nrep = 1000)

  # Run via nthreads = 1 argument
  res_arg1 <- spurtest_i0(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42, q = 6, nrep = 1000, nthreads = 1)

  # Run single-threaded via setter
  spur_set_num_threads(1)
  res_single <- spurtest_i0(am ~ 1, data = spur_example, lon = "lon", lat = "lat", seed = 42, q = 6, nrep = 1000)

  expect_equal(res_arg1$teststat, res_multi$teststat, tolerance = 1e-14)
  expect_equal(res_arg1$pvalue, res_multi$pvalue, tolerance = 1e-14)
  expect_equal(res_arg1$cv, res_multi$cv, tolerance = 1e-14)
  expect_equal(res_single$teststat, res_multi$teststat, tolerance = 1e-14)
  expect_equal(res_single$pvalue, res_multi$pvalue, tolerance = 1e-14)
  expect_equal(res_single$cv, res_multi$cv, tolerance = 1e-14)
})

test_that("full spur() wrapper respects nthreads argument", {
  on.exit(spur_set_num_threads("auto"), add = TRUE)

  data(spur_example, package = "spuR")

  fit_auto <- spur(am ~ fracblack, data = spur_example, lon = "lon", lat = "lat", q = 8L, nrep = 500L, seed = 42L)
  fit_t1 <- spur(am ~ fracblack, data = spur_example, lon = "lon", lat = "lat", q = 8L, nrep = 500L, seed = 42L, nthreads = 1)

  expect_equal(fit_t1$tests$i0$teststat, fit_auto$tests$i0$teststat, tolerance = 1e-14)
  expect_equal(fit_t1$tests$i0$pvalue, fit_auto$tests$i0$pvalue, tolerance = 1e-14)
  expect_equal(fit_t1$fits$levels$scpc$scpcstats, fit_auto$fits$levels$scpc$scpcstats, tolerance = 1e-14)
  expect_equal(fit_t1$fits$transformed$scpc$scpcstats, fit_auto$fits$transformed$scpc$scpcstats, tolerance = 1e-14)
})
