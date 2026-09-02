## Benchmark AR(1) covariance performance with repeated full fits and retapes.
##
## Run from the repository root, for example:
##   Rscript performance/covariance/benchmark-ar1-extensive.R
##
## Optional controls:
##   AR1_EXTENSIVE_REPS=5
##   AR1_EXTENSIVE_BLOCKS=60,180
##   AR1_EXTENSIVE_GROUPS=8,18
##   AR1_EXTENSIVE_EVAL_MAX=100

script_path <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
  } else {
    normalizePath("performance/covariance/benchmark-ar1-extensive.R",
                  mustWork = TRUE)
  }
})

repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                           mustWork = TRUE)
pkg_dir <- file.path(repo_root, "glmmTMB")

parse_ints <- function(x, default) {
  x <- Sys.getenv(x, default)
  ans <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(ans)) || any(ans < 1L)) {
    stop("Expected a comma-separated list of positive integers.")
  }
  ans
}

reps <- as.integer(Sys.getenv("AR1_EXTENSIVE_REPS", "3"))
block_sizes <- parse_ints("AR1_EXTENSIVE_BLOCKS", "60,180")
group_sizes <- parse_ints("AR1_EXTENSIVE_GROUPS", "8,18")
eval_max <- as.integer(Sys.getenv("AR1_EXTENSIVE_EVAL_MAX", "100"))

if (is.na(reps) || reps < 1L) {
  stop("AR1_EXTENSIVE_REPS must be a positive integer.")
}
if (length(block_sizes) != length(group_sizes)) {
  stop("AR1_EXTENSIVE_BLOCKS and AR1_EXTENSIVE_GROUPS must have equal length.")
}
if (is.na(eval_max) || eval_max < 1L) {
  stop("AR1_EXTENSIVE_EVAL_MAX must be a positive integer.")
}

make_ar1_data <- function(block_n, groups, seed) {
  set.seed(seed)
  d <- expand.grid(
    time = seq_len(block_n),
    group = factor(seq_len(groups))
  )
  d$time_fac <- factor(d$time, levels = seq_len(block_n), ordered = TRUE)
  d$x <- stats::rnorm(nrow(d))
  group_eff <- stats::rnorm(groups, sd = 0.8)
  ar1_eff <- unlist(
    replicate(
      groups,
      as.numeric(stats::arima.sim(list(ar = 0.65), n = block_n, sd = 0.7)),
      simplify = FALSE
    )
  )
  d$y <- 2 + 0.4 * d$x + group_eff[d$group] + ar1_eff +
    stats::rnorm(nrow(d), sd = 0.3)
  d
}

use_rtmb <- function(value) {
  ns <- asNamespace("glmmTMB")
  if (exists("useRTMB", envir = ns, inherits = FALSE)) {
    get("useRTMB", envir = ns)(value)
  }
  invisible(value)
}

fit_one <- function(d, backend) {
  use_rtmb(backend == "local_rtmb")
  glmmTMB::glmmTMB(
    y ~ x + (1 | group) + ar1(time_fac + 0 | group),
    data = d,
    control = glmmTMB::glmmTMBControl(
      optCtrl = list(iter.max = eval_max, eval.max = eval_max),
      conv_check = "skip"
    )
  )
}

time_backend <- function(backend) {
  backend_label <- switch(
    backend,
    installed_tmb = "installed TMB",
    local_rtmb = "local RTMB",
    stop("Unknown backend: ", backend)
  )

  if (backend == "local_rtmb") {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Install the 'pkgload' package before running local RTMB benchmarks.")
    }
    pkgload::load_all(pkg_dir, quiet = TRUE)
  } else {
    library(glmmTMB)
  }

  out <- list()
  k <- 0L
  for (case_id in seq_along(block_sizes)) {
    block_n <- block_sizes[case_id]
    groups <- group_sizes[case_id]
    d <- make_ar1_data(block_n, groups, seed = 1000L + case_id)

    ## Warmup one fit per case/backend before timed repetitions.
    invisible(fit_one(d, backend))

    for (run in seq_len(reps)) {
      fit_time <- system.time(fit <- fit_one(d, backend))[["elapsed"]]
      retape_time <- system.time(fit$obj$env$retape())[["elapsed"]]
      evaluations <- fit$fit$evaluations
      k <- k + 1L
      out[[k]] <- data.frame(
        backend = backend_label,
        case = paste0("block=", block_n, ", groups=", groups),
        block_n = block_n,
        groups = groups,
        observations = nrow(d),
        run = run,
        fit_sec = fit_time,
        retape_sec = retape_time,
        fn_evals = if (length(evaluations) >= 1L) evaluations[[1L]] else NA,
        gr_evals = if (length(evaluations) >= 2L) evaluations[[2L]] else NA,
        convergence = fit$fit$convergence,
        stringsAsFactors = FALSE
      )
      cat(
        backend_label, " ", out[[k]]$case, " run ", run,
        ": fit=", sprintf("%.3f", fit_time),
        "s retape=", sprintf("%.3f", retape_time), "s\n",
        sep = ""
      )
    }
  }
  do.call(rbind, out)
}

run_child <- function() {
  backend <- Sys.getenv("AR1_EXTENSIVE_CHILD_BACKEND")
  output <- Sys.getenv("AR1_EXTENSIVE_CHILD_OUTPUT")
  if (!backend %in% c("installed_tmb", "local_rtmb")) {
    stop("Child backend must be 'installed_tmb' or 'local_rtmb'.")
  }
  if (!nzchar(output)) {
    stop("Child output path is missing.")
  }
  saveRDS(time_backend(backend), output)
}

se <- function(x) stats::sd(x) / sqrt(length(x))

summarize <- function(rows) {
  cases <- unique(rows$case)
  backends <- unique(rows$backend)
  out <- list()
  k <- 0L
  for (case in cases) {
    for (backend in backends) {
      z <- rows[rows$case == case & rows$backend == backend, ]
      k <- k + 1L
      out[[k]] <- data.frame(
        case = case,
        backend = backend,
        block_n = z$block_n[1L],
        groups = z$groups[1L],
        observations = z$observations[1L],
        fit_mean_sec = mean(z$fit_sec),
        fit_se_sec = se(z$fit_sec),
        retape_mean_sec = mean(z$retape_sec),
        retape_se_sec = se(z$retape_sec),
        mean_fn_evals = mean(z$fn_evals, na.rm = TRUE),
        mean_gr_evals = mean(z$gr_evals, na.rm = TRUE),
        status = if (all(z$convergence == 0L)) "OK" else "CHECK",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

wide_summary <- function(summary_rows) {
  tmb <- summary_rows[summary_rows$backend == "installed TMB", ]
  rtmb <- summary_rows[summary_rows$backend == "local RTMB", ]
  names(tmb)[names(tmb) %in% c("fit_mean_sec", "fit_se_sec",
                               "retape_mean_sec", "retape_se_sec")] <-
    paste0("tmb_", names(tmb)[names(tmb) %in% c("fit_mean_sec", "fit_se_sec",
                                                "retape_mean_sec",
                                                "retape_se_sec")])
  names(rtmb)[names(rtmb) %in% c("fit_mean_sec", "fit_se_sec",
                                 "retape_mean_sec", "retape_se_sec")] <-
    paste0("rtmb_", names(rtmb)[names(rtmb) %in% c("fit_mean_sec",
                                                   "fit_se_sec",
                                                   "retape_mean_sec",
                                                   "retape_se_sec")])
  ans <- merge(
    tmb[c("case", "block_n", "groups", "observations", "tmb_fit_mean_sec",
          "tmb_fit_se_sec", "tmb_retape_mean_sec", "tmb_retape_se_sec",
          "mean_fn_evals", "mean_gr_evals", "status")],
    rtmb[c("case", "rtmb_fit_mean_sec", "rtmb_fit_se_sec",
           "rtmb_retape_mean_sec", "rtmb_retape_se_sec", "mean_fn_evals",
           "mean_gr_evals", "status")],
    by = "case",
    suffixes = c("_tmb", "_rtmb")
  )
  ans$fit_ratio <- ans$rtmb_fit_mean_sec / ans$tmb_fit_mean_sec
  ans$retape_ratio <- ans$rtmb_retape_mean_sec / ans$tmb_retape_mean_sec
  ans
}

fmt <- function(x) formatC(x, format = "f", digits = 3)

print_summary <- function(x) {
  display <- x[order(-x$fit_ratio), ]
  display <- data.frame(
    case = display$case,
    obs = display$observations,
    tmb_fit = paste0(fmt(display$tmb_fit_mean_sec), " +/- ",
                     fmt(display$tmb_fit_se_sec)),
    rtmb_fit = paste0(fmt(display$rtmb_fit_mean_sec), " +/- ",
                      fmt(display$rtmb_fit_se_sec)),
    fit_ratio = paste0(fmt(display$fit_ratio), "x"),
    tmb_retape = paste0(fmt(display$tmb_retape_mean_sec), " +/- ",
                        fmt(display$tmb_retape_se_sec)),
    rtmb_retape = paste0(fmt(display$rtmb_retape_mean_sec), " +/- ",
                         fmt(display$rtmb_retape_se_sec)),
    retape_ratio = paste0(fmt(display$retape_ratio), "x"),
    fn_evals = paste0(round(display$mean_fn_evals_tmb), "/",
                      round(display$mean_fn_evals_rtmb)),
    gr_evals = paste0(round(display$mean_gr_evals_tmb), "/",
                      round(display$mean_gr_evals_rtmb)),
    status = paste(display$status_tmb, display$status_rtmb, sep = "/"),
    stringsAsFactors = FALSE
  )
  print(display, row.names = FALSE)
}

run_parent <- function() {
  if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    stop("Run this script from the repository root.")
  }

  tmb_file <- tempfile(fileext = ".rds")
  rtmb_file <- tempfile(fileext = ".rds")
  on.exit(unlink(c(tmb_file, rtmb_file)), add = TRUE)

  for (backend in c("installed_tmb", "local_rtmb")) {
    output <- if (backend == "installed_tmb") tmb_file else rtmb_file
    status <- system2(
      "Rscript",
      c("--vanilla", script_path),
      env = c(
        paste0("AR1_EXTENSIVE_CHILD_BACKEND=", backend),
        paste0("AR1_EXTENSIVE_CHILD_OUTPUT=", output),
        paste0("AR1_EXTENSIVE_REPS=", reps),
        paste0("AR1_EXTENSIVE_BLOCKS=", paste(block_sizes, collapse = ",")),
        paste0("AR1_EXTENSIVE_GROUPS=", paste(group_sizes, collapse = ",")),
        paste0("AR1_EXTENSIVE_EVAL_MAX=", eval_max)
      )
    )
    if (!identical(status, 0L)) {
      stop(backend, " benchmark failed with exit status ", status)
    }
  }

  rows <- rbind(readRDS(tmb_file), readRDS(rtmb_file))
  cat("\nRaw timings\n")
  print(rows, row.names = FALSE)

  cat("\nAR1 extensive benchmark summary\n")
  cat("Fit columns are mean +/- standard error over ", reps,
      " timed runs.\n", sep = "")
  cat("fn_evals and gr_evals are shown as TMB/RTMB means.\n\n")
  print_summary(wide_summary(summarize(rows)))

  invisible(rows)
}

if (nzchar(Sys.getenv("AR1_EXTENSIVE_CHILD_BACKEND"))) {
  run_child()
} else {
  run_parent()
}
