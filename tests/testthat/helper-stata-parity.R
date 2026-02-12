spur_stata_bin <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }

    env_bin <- Sys.getenv("SPUR_STATA_BIN", unset = "")
    candidates <- c(
      env_bin,
      Sys.which(c("stata-se", "stata-mp", "stata")),
      "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se",
      "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp",
      "/Applications/Stata/Stata.app/Contents/MacOS/stata"
    )
    candidates <- unique(candidates[nzchar(candidates)])
    hit <- candidates[file.exists(candidates)]
    cached <<- if (length(hit)) hit[[1]] else ""
    cached
  }
})

spur_stata_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

spur_stata_rc_from_log <- function(log_lines) {
  matches <- regmatches(log_lines, gregexpr("r\\(([0-9]+)\\);", log_lines, perl = TRUE))
  flat <- unlist(matches, use.names = FALSE)
  if (!length(flat)) {
    return(0L)
  }
  last <- flat[[length(flat)]]
  as.integer(sub("^r\\(([0-9]+)\\);$", "\\1", last))
}

spur_run_stata_do <- function(do_lines) {
  wd <- tempfile("stata_parity_")
  dir.create(wd, recursive = TRUE)
  do_file <- file.path(wd, "run.do")
  log_file <- file.path(wd, "run.log")
  writeLines(do_lines, do_file)

  old_wd <- setwd(wd)
  on.exit(setwd(old_wd), add = TRUE)
  suppressWarnings(system2(spur_stata_bin(), c("-b", "do", basename(do_file))))

  log_lines <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character(0)
  rc <- spur_stata_rc_from_log(log_lines)
  list(rc = rc, log_lines = log_lines, log_file = log_file, wd = wd)
}

spur_skip_if_no_stata_or_spur <- local({
  checked <- FALSE
  available <- FALSE

  function() {
    if (checked && !available) {
      testthat::skip("SPUR Stata package not available in current Stata setup.")
    }
    if (checked && available) {
      return(invisible(NULL))
    }

    if (!nzchar(spur_stata_bin())) {
      checked <<- TRUE
      available <<- FALSE
      testthat::skip("Stata executable not found.")
    }

    probe <- spur_run_stata_do(c(
      "clear all",
      "cap which spurtest",
      "if _rc exit 200",
      "cap which spurtransform",
      "if _rc exit 201",
      "cap which spurhalflife",
      "if _rc exit 202",
      "exit, clear"
    ))

    checked <<- TRUE
    available <<- probe$rc == 0L
    if (!available) {
      testthat::skip("SPUR Stata package not available in current Stata setup.")
    }
    invisible(NULL)
  }
})

spur_fixture_data <- function(n = 180L) {
  candidates <- c(
    testthat::test_path("fixtures", "chetty_data_1.csv"),
    file.path("tests", "testthat", "fixtures", "chetty_data_1.csv")
  )
  fixture <- candidates[file.exists(candidates)][1]
  if (is.na(fixture) || !nzchar(fixture)) {
    stop("Missing Chetty fixture: tests/fixtures/chetty_data_1.csv")
  }
  df <- utils::read.csv(fixture, check.names = FALSE, stringsAsFactors = FALSE)
  df <- df[seq_len(min(n, nrow(df))), , drop = FALSE]
  df$row_id <- seq_len(nrow(df))
  df$s_1 <- df$lat
  df$s_2 <- df$lon
  df
}

spur_run_stata_spurtest <- function(df, test_type, varlist, q, nrep, latlong) {
  in_csv <- tempfile(fileext = ".csv")
  out_csv <- tempfile(fileext = ".csv")
  utils::write.csv(df, in_csv, row.names = FALSE, na = "")

  latlong_opt <- if (isTRUE(latlong)) " latlong" else ""
  do_lines <- c(
    "clear all",
    sprintf("import delimited using \"%s\", varnames(1) clear", spur_stata_path(in_csv)),
    "set seed 123",
    sprintf(
      "spurtest %s %s, q(%d) nrep(%d)%s",
      test_type,
      varlist,
      as.integer(q),
      as.integer(nrep),
      latlong_opt
    ),
    "scalar p_res = r(p)",
    "scalar teststat_res = r(teststat)",
    "scalar ha_res = r(ha_param)",
    "matrix cv_res = r(cv)",
    "clear",
    "set obs 1",
    "gen pvalue = p_res",
    "gen teststat = teststat_res",
    "gen ha_param = ha_res",
    "gen cv1 = cv_res[1,1]",
    "gen cv2 = cv_res[2,1]",
    "gen cv3 = cv_res[3,1]",
    sprintf("export delimited using \"%s\", replace", spur_stata_path(out_csv)),
    "exit, clear"
  )

  run <- spur_run_stata_do(do_lines)
  if (run$rc != 0L) {
    stop(sprintf("Stata spurtest run failed (rc=%s). Log: %s", run$rc, run$log_file))
  }

  out <- utils::read.csv(out_csv, check.names = FALSE)
  out[1, , drop = FALSE]
}

spur_run_stata_spurhalflife <- function(df, var, q, nrep, level, latlong, normdist) {
  in_csv <- tempfile(fileext = ".csv")
  out_csv <- tempfile(fileext = ".csv")
  utils::write.csv(df, in_csv, row.names = FALSE, na = "")

  latlong_opt <- if (isTRUE(latlong)) " latlong" else ""
  normdist_opt <- if (isTRUE(normdist)) " normdist" else ""
  do_lines <- c(
    "clear all",
    sprintf("import delimited using \"%s\", varnames(1) clear", spur_stata_path(in_csv)),
    "set seed 123",
    sprintf(
      "spurhalflife %s, q(%d) nrep(%d) level(%s)%s%s",
      var,
      as.integer(q),
      as.integer(nrep),
      format(level),
      latlong_opt,
      normdist_opt
    ),
    "scalar ci_l_res = r(ci_l)",
    "scalar ci_u_res = r(ci_u)",
    "scalar max_dist_res = r(max_dist)",
    "clear",
    "set obs 1",
    "gen ci_l = ci_l_res",
    "gen ci_u = ci_u_res",
    "gen max_dist = max_dist_res",
    sprintf("export delimited using \"%s\", replace", spur_stata_path(out_csv)),
    "exit, clear"
  )

  run <- spur_run_stata_do(do_lines)
  if (run$rc != 0L) {
    stop(sprintf("Stata spurhalflife run failed (rc=%s). Log: %s", run$rc, run$log_file))
  }

  out <- utils::read.csv(out_csv, check.names = FALSE)
  out[1, , drop = FALSE]
}

spur_run_stata_spurtransform <- function(df,
                                         vars,
                                         prefix,
                                         transformation,
                                         latlong,
                                         radius = NULL,
                                         clustvar = NULL,
                                         separately = FALSE,
                                         replace = FALSE) {
  in_csv <- tempfile(fileext = ".csv")
  out_csv <- tempfile(fileext = ".csv")
  utils::write.csv(df, in_csv, row.names = FALSE, na = "")

  opts <- c(
    sprintf("prefix(\"%s\")", prefix),
    sprintf("transformation(%s)", transformation)
  )
  if (!is.null(radius)) {
    opts <- c(opts, sprintf("radius(%s)", format(radius)))
  }
  if (!is.null(clustvar)) {
    opts <- c(opts, sprintf("clustvar(%s)", clustvar))
  }
  if (isTRUE(latlong)) {
    opts <- c(opts, "latlong")
  }
  if (isTRUE(replace)) {
    opts <- c(opts, "replace")
  }
  if (isTRUE(separately)) {
    opts <- c(opts, "separately")
  }

  out_vars <- paste0(prefix, vars)
  do_lines <- c(
    "clear all",
    sprintf("import delimited using \"%s\", varnames(1) clear", spur_stata_path(in_csv)),
    sprintf("spurtransform %s, %s", paste(vars, collapse = " "), paste(opts, collapse = " ")),
    sprintf("keep row_id %s", paste(out_vars, collapse = " ")),
    "sort row_id",
    sprintf("export delimited using \"%s\", replace", spur_stata_path(out_csv)),
    "exit, clear"
  )

  run <- spur_run_stata_do(do_lines)
  if (run$rc != 0L) {
    stop(sprintf("Stata spurtransform run failed (rc=%s). Log: %s", run$rc, run$log_file))
  }

  utils::read.csv(out_csv, check.names = FALSE)
}

spur_run_stata_command_rc <- function(df, command_line) {
  in_csv <- tempfile(fileext = ".csv")
  utils::write.csv(df, in_csv, row.names = FALSE, na = "")

  do_lines <- c(
    "clear all",
    sprintf("import delimited using \"%s\", varnames(1) clear", spur_stata_path(in_csv)),
    command_line,
    "exit, clear"
  )
  spur_run_stata_do(do_lines)$rc
}
