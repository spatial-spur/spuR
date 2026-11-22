spur_mw_root <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }

    env_root <- Sys.getenv("SPUR_MW_ROOT", unset = "")
    test_root <- tryCatch(
      testthat::test_path("..", "..", "ReplicationPackage_MS21654_MuellerWatson_r2"),
      error = function(...) "ReplicationPackage_MS21654_MuellerWatson_r2"
    )
    candidates <- c(
      env_root,
      test_root,
      file.path(getwd(), "ReplicationPackage_MS21654_MuellerWatson_r2"),
      "ReplicationPackage_MS21654_MuellerWatson_r2"
    )
    candidates <- unique(candidates[nzchar(candidates)])
    hit <- candidates[file.exists(candidates)]
    cached <<- if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = FALSE) else ""
    cached
  }
})

spur_matlab_bin <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }

    env_bin <- Sys.getenv("SPUR_MATLAB_BIN", unset = "")
    app_bins <- character(0)
    if (dir.exists("/Applications")) {
      apps <- sort(
        dir("/Applications", pattern = "^MATLAB_R[0-9]{4}[ab]\\.app$", full.names = TRUE),
        decreasing = TRUE
      )
      app_bins <- file.path(apps, "bin", "matlab")
    }
    linux_bins <- Sys.glob("/usr/local/MATLAB/*/bin/matlab")

    candidates <- c(
      env_bin,
      Sys.which(c("matlab", "octave", "octave-cli")),
      app_bins,
      linux_bins
    )
    candidates <- unique(candidates[nzchar(candidates)])
    hit <- candidates[file.exists(candidates)]
    cached <<- if (length(hit)) hit[[1]] else ""
    cached
  }
})

spur_matlab_engine <- function() {
  if (grepl("octave", basename(spur_matlab_bin()), ignore.case = TRUE)) {
    return("octave")
  }
  "matlab"
}

spur_matlab_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

spur_matlab_quote <- function(x) {
  gsub("'", "''", x, fixed = TRUE)
}

spur_matlab_addpath_line <- function(path) {
  sprintf("addpath('%s');", spur_matlab_quote(spur_matlab_path(path)))
}

spur_matlab_distmat_lines <- function(latlong_flag,
                                      s_var = "s",
                                      n_var = "n_obs",
                                      dist_raw_var = "distmat_raw",
                                      dist_var = "distmat") {
  if (isTRUE(as.logical(latlong_flag))) {
    return(c(
      sprintf("%s = size(%s,1);", n_var, s_var),
      sprintf("%s = NaN(%s,%s);", dist_raw_var, n_var, n_var),
      "c_const = 3.14159265359/180;",
      sprintf("for jj = 1:%s", n_var),
      sprintf("  tmp = repmat(%s(jj,:), %s, 1);", s_var, n_var),
      sprintf("  d = 0.5 * c_const * (tmp - %s);", s_var),
      sprintf("  %s(:,jj) = asin(sqrt((sin(d(:,1))).^2 + cos(c_const * %s(jj,1)) .* (cos(c_const * %s(:,1)) .* (sin(d(:,2))).^2))) / 3.14159265359;", dist_raw_var, s_var, s_var),
      "end",
      sprintf("%s = %s / max(max(%s));", dist_var, dist_raw_var, dist_raw_var)
    ))
  }

  c(
    sprintf("%s = size(%s,1);", n_var, s_var),
    sprintf("d_dim = size(%s,2);", s_var),
    sprintf("dist_acc = zeros(%s,%s);", n_var, n_var),
    "for jj = 1:d_dim",
    sprintf("  dist_acc = dist_acc + (repmat(%s(:,jj),1,%s) - repmat(%s(:,jj)',%s,1)).^2;", s_var, n_var, s_var, n_var),
    "end",
    sprintf("%s = sqrt(dist_acc);", dist_raw_var),
    sprintf("%s = %s / max(max(%s));", dist_var, dist_raw_var, dist_raw_var)
  )
}

spur_run_matlab_process <- function(args, wd) {
  old_wd <- setwd(wd)
  on.exit(setwd(old_wd), add = TRUE)

  output <- suppressWarnings(
    system2(spur_matlab_bin(), args = args, stdout = TRUE, stderr = TRUE)
  )
  rc <- attr(output, "status")
  if (is.null(rc)) {
    rc <- 0L
  }

  list(
    rc = as.integer(rc),
    log_lines = output,
    wd = wd
  )
}

spur_run_matlab_script <- function(script_lines) {
  wd <- tempfile("matlab_parity_")
  dir.create(wd, recursive = TRUE)
  script_file <- file.path(wd, "run_parity.m")
  writeLines(script_lines, script_file)
  script_stem <- tools::file_path_sans_ext(basename(script_file))

  engine <- spur_matlab_engine()
  run <- if (identical(engine, "octave")) {
    spur_run_matlab_process(
      args = c("--quiet", "--no-gui", basename(script_file)),
      wd = wd
    )
  } else {
    run_batch <- spur_run_matlab_process(
      args = c("-batch", script_stem),
      wd = wd
    )
    if (run_batch$rc != 0L && any(grepl("-batch", run_batch$log_lines, ignore.case = TRUE))) {
      legacy_cmd <- sprintf(
        "try, %s; catch ME, disp(ME.message); exit(1); end; exit(0);",
        script_stem
      )
      spur_run_matlab_process(
        args = c("-nodisplay", "-nosplash", "-nodesktop", "-r", legacy_cmd),
        wd = wd
      )
    } else {
      run_batch
    }
  }

  run$script_file <- script_file
  run
}

spur_skip_if_no_matlab_or_mw <- local({
  checked <- FALSE
  available <- FALSE

  function() {
    if (checked && !available) {
      testthat::skip("Matlab + MW replication package not available in current setup.")
    }
    if (checked && available) {
      return(invisible(NULL))
    }

    if (!nzchar(spur_matlab_bin())) {
      checked <<- TRUE
      available <<- FALSE
      testthat::skip("Matlab executable not found.")
    }

    mw_root <- spur_mw_root()
    needed_file <- file.path(mw_root, "matlab_functions", "lbm_gls_matrix.m")
    if (!nzchar(mw_root) || !file.exists(needed_file)) {
      checked <<- TRUE
      available <<- FALSE
      testthat::skip("MW Matlab replication package not found.")
    }

    probe <- spur_run_matlab_script(c(
      "try",
      spur_matlab_addpath_line(file.path(mw_root, "Matlab", "functions")),
      spur_matlab_addpath_line(file.path(mw_root, "matlab_functions")),
      spur_matlab_addpath_line(file.path(mw_root, "cscpc")),
      "ok = exist('lbm_gls_matrix', 'file') == 2;",
      "if ~ok, error('MW Matlab SPUR functions are not available on path.'); end",
      "disp('MW Matlab parity probe OK');",
      "exit(0);",
      "catch ME",
      "disp(ME.message);",
      "exit(1);",
      "end"
    ))

    checked <<- TRUE
    available <<- probe$rc == 0L
    if (!available) {
      testthat::skip("Matlab + MW replication package not available in current setup.")
    }
    invisible(NULL)
  }
})

spur_run_matlab_spurtransform <- function(df,
                                          vars,
                                          prefix,
                                          transformation,
                                          latlong,
                                          radius = NULL,
                                          clustvar = NULL,
                                          separately = FALSE,
                                          replace = FALSE) {
  if (!identical(transformation, "lbmgls")) {
    stop("Matlab 3-way transform parity currently supports transformation = 'lbmgls' only.")
  }
  if (!is.null(radius)) {
    stop("radius is not used for transformation = 'lbmgls'.")
  }
  if (!is.null(clustvar)) {
    stop("clustvar is not used for transformation = 'lbmgls'.")
  }
  if (!all(c("row_id", "s_1", "s_2", vars) %in% names(df))) {
    stop("Missing required columns for Matlab spurtransform parity run.")
  }

  in_csv <- tempfile(fileext = ".csv")
  out_csv <- tempfile(fileext = ".csv")
  utils::write.csv(df, in_csv, row.names = FALSE, na = "")

  latlong_flag <- if (isTRUE(latlong)) 1L else 0L
  vars_q <- vapply(vars, spur_matlab_quote, character(1))
  vars_cell <- sprintf("{%s}", paste(sprintf("'%s'", vars_q), collapse = ", "))

  var_blocks <- unlist(lapply(vars, function(v) {
    out_var <- paste0(prefix, v)
    c(
      sprintf("varname = '%s';", spur_matlab_quote(v)),
      sprintf("outname = '%s';", spur_matlab_quote(out_var)),
      if (isTRUE(separately)) {
        "touse = ~isnan(in_tbl.(varname));"
      } else {
        "touse = touse_base;"
      },
      "if sum(touse) < 2, error('Not enough observations after filtering for transform.'); end",
      "y = in_tbl.(varname);",
      "y_use = y(touse);",
      "s_use = s_all(touse,:);",
      spur_matlab_distmat_lines(latlong_flag, s_var = "s_use", n_var = "n_use", dist_raw_var = "dist_use_raw", dist_var = "dist_use"),
      "sigma_lbm = 0.5 * (repmat(dist_use(:,1),1,n_use) + repmat(dist_use(1,:),n_use,1) - dist_use);",
      "sigma_lbm_dm = sigma_lbm - repmat(mean(sigma_lbm,1),n_use,1);",
      "sigma_lbm_dm = sigma_lbm_dm - repmat(mean(sigma_lbm_dm,2),1,n_use);",
      "[V,d] = eig(sigma_lbm_dm,'vector');",
      "[eval_sorted, idx] = sort(d,'descend');",
      "evec = V(:,idx);",
      "keep = eval_sorted > 1.0e-10;",
      "eval_sorted = eval_sorted(keep == 1);",
      "evec = evec(:,keep == 1);",
      "H = evec * diag(1 ./ sqrt(eval_sorted)) * evec';",
      "hy = H * y_use;",
      "hy = hy - mean(hy);",
      "out_col = NaN(height(in_tbl),1);",
      "out_col(touse) = hy;",
      "out_tbl.(outname) = out_col;"
    )
  }), use.names = FALSE)

  script <- c(
    "try",
    sprintf("in_tbl = readtable('%s');", spur_matlab_quote(spur_matlab_path(in_csv))),
    "if height(in_tbl) == 0, error('Input table is empty.'); end",
    "if ~ismember('row_id', in_tbl.Properties.VariableNames), error('row_id missing from input table.'); end",
    "s_all = [in_tbl.s_1, in_tbl.s_2];",
    sprintf("vars = %s;", vars_cell),
    "touse_base = true(height(in_tbl),1);",
    "for iv = 1:numel(vars)",
    "  touse_base = touse_base & ~isnan(in_tbl.(vars{iv}));",
    "end",
    "out_tbl = table(in_tbl.row_id, 'VariableNames', {'row_id'});",
    var_blocks,
    sprintf("writetable(out_tbl, '%s');", spur_matlab_quote(spur_matlab_path(out_csv))),
    "exit(0);",
    "catch ME",
    "disp(ME.message);",
    "exit(1);",
    "end"
  )

  run <- spur_run_matlab_script(script)
  if (run$rc != 0L) {
    log_excerpt <- if (length(run$log_lines)) {
      paste(run$log_lines, collapse = " | ")
    } else {
      "<no matlab output>"
    }
    stop(sprintf(
      "Matlab spurtransform run failed (rc=%s). Script: %s. Log: %s",
      run$rc,
      run$script_file,
      log_excerpt
    ))
  }

  utils::read.csv(out_csv, check.names = FALSE)
}
