## Benchmark individual glmmTMB tests under TMB and RTMB.
##
## Run from the repository root:
##   Rscript performance/test-suite/benchmark-tests.R
##
## Optional controls:
##   RTMB_TEST_BENCHMARK_FILTER=rtmb-gaussian
##   RTMB_TEST_BENCHMARK_TOP=30
##   RTMB_TEST_BENCHMARK_RATIO_MIN_TMB=0.05

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package before running this benchmark.")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install the 'testthat' package before running this benchmark.")
}
if (!requireNamespace("R6", quietly = TRUE)) {
  stop("Install the 'R6' package before running this benchmark.")
}

script_path <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
  } else {
    normalizePath("performance/test-suite/benchmark-tests.R", mustWork = TRUE)
  }
})

repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                           mustWork = TRUE)
pkg_dir <- file.path(repo_root, "glmmTMB")
test_dir <- file.path(pkg_dir, "tests", "testthat")

filter <- Sys.getenv("RTMB_TEST_BENCHMARK_FILTER", "")
top_n <- as.integer(Sys.getenv("RTMB_TEST_BENCHMARK_TOP", "30"))
ratio_min_tmb <- as.numeric(Sys.getenv("RTMB_TEST_BENCHMARK_RATIO_MIN_TMB",
                                       "0.05"))

if (is.na(top_n) || top_n < 1L) {
  stop("RTMB_TEST_BENCHMARK_TOP must be a positive integer.")
}
if (is.na(ratio_min_tmb) || ratio_min_tmb < 0) {
  stop("RTMB_TEST_BENCHMARK_RATIO_MIN_TMB must be non-negative.")
}

TimingReporter <- R6::R6Class(
  "TimingReporter",
  inherit = testthat::Reporter,
  public = list(
    rows = NULL,
    current_file = NULL,
    current_context = NULL,
    current_test = NULL,
    current_start = NULL,
    current_counts = NULL,
    test_index = 0L,

    initialize = function(file = stdout()) {
      super$initialize(file = file)
      self$rows <- list()
      self$current_counts <- self$empty_counts()
    },

    empty_counts = function() {
      c(expectations = 0L, failures = 0L, errors = 0L, warnings = 0L,
        skips = 0L)
    },

    start_file = function(filename) {
      self$current_file <- basename(filename)
    },

    start_test = function(context, test) {
      self$test_index <- self$test_index + 1L
      self$current_context <- context
      self$current_test <- test
      self$current_start <- proc.time()[["elapsed"]]
      self$current_counts <- self$empty_counts()
    },

    add_result = function(context, test, result) {
      self$current_counts[["expectations"]] <-
        self$current_counts[["expectations"]] + 1L

      if (inherits(result, "expectation_skip")) {
        self$current_counts[["skips"]] <- self$current_counts[["skips"]] + 1L
      } else if (inherits(result, "expectation_warning")) {
        self$current_counts[["warnings"]] <-
          self$current_counts[["warnings"]] + 1L
      } else if (inherits(result, "expectation_error")) {
        self$current_counts[["errors"]] <- self$current_counts[["errors"]] + 1L
      } else if (inherits(result, "expectation_failure")) {
        self$current_counts[["failures"]] <-
          self$current_counts[["failures"]] + 1L
      }
    },

    end_test = function(context, test) {
      elapsed <- proc.time()[["elapsed"]] - self$current_start
      counts <- as.list(self$current_counts)
      self$rows[[length(self$rows) + 1L]] <- data.frame(
        test_index = self$test_index,
        file = self$current_file,
        context = self$current_context,
        test = self$current_test,
        seconds = elapsed,
        expectations = counts$expectations,
        failures = counts$failures,
        errors = counts$errors,
        warnings = counts$warnings,
        skips = counts$skips,
        stringsAsFactors = FALSE
      )
    }
  )
)

run_backend <- function(backend) {
  pkgload::load_all(pkg_dir, quiet = TRUE)

  old_use_rtmb <- glmmTMB::useRTMB()
  on.exit(glmmTMB::useRTMB(old_use_rtmb), add = TRUE)
  glmmTMB::useRTMB(backend == "RTMB")

  reporter <- TimingReporter$new()
  testthat::test_dir(
    test_dir,
    filter = if (nzchar(filter)) filter else NULL,
    reporter = reporter,
    load_helpers = TRUE,
    stop_on_failure = FALSE,
    stop_on_warning = FALSE,
    package = "glmmTMB",
    load_package = "none"
  )

  out <- if (length(reporter$rows) > 0L) {
    do.call(rbind, reporter$rows)
  } else {
    data.frame(
      test_index = integer(),
      file = character(),
      context = character(),
      test = character(),
      seconds = numeric(),
      expectations = integer(),
      failures = integer(),
      errors = integer(),
      warnings = integer(),
      skips = integer(),
      stringsAsFactors = FALSE
    )
  }
  out$backend <- backend
  out
}

run_child <- function() {
  backend <- Sys.getenv("RTMB_TEST_BENCHMARK_CHILD_BACKEND")
  output <- Sys.getenv("RTMB_TEST_BENCHMARK_CHILD_OUTPUT")
  if (!backend %in% c("TMB", "RTMB")) {
    stop("Child backend must be TMB or RTMB.")
  }
  if (!nzchar(output)) {
    stop("Child output path is missing.")
  }
  saveRDS(run_backend(backend), output)
}

clip <- function(x, width) {
  x <- as.character(x)
  too_wide <- nchar(x, type = "width") > width
  x[too_wide] <- paste0(substr(x[too_wide], 1L, width - 1L), "~")
  x
}

format_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

print_table <- function(x, columns, widths) {
  x <- x[columns]
  for (nm in names(widths)) {
    x[[nm]] <- clip(x[[nm]], widths[[nm]])
  }
  header <- mapply(format, names(x), width = widths[names(x)],
                   justify = "left", USE.NAMES = FALSE)
  cat(paste(header, collapse = "  "), "\n", sep = "")
  cat(paste(vapply(widths[names(x)], function(w) paste(rep("-", w),
                                                       collapse = ""),
                   character(1L)), collapse = "  "), "\n", sep = "")
  for (i in seq_len(nrow(x))) {
    row <- mapply(format, x[i, ], width = widths[names(x)],
                  justify = "left", USE.NAMES = FALSE)
    cat(paste(row, collapse = "  "), "\n", sep = "")
  }
}

compare_backends <- function(tmb, rtmb) {
  tmb$status <- ifelse(tmb$failures + tmb$errors > 0L, "FAIL", "OK")
  rtmb$status <- ifelse(rtmb$failures + rtmb$errors > 0L, "FAIL", "OK")
  names(tmb)[names(tmb) == "seconds"] <- "tmb_sec"
  names(rtmb)[names(rtmb) == "seconds"] <- "rtmb_sec"
  names(tmb)[names(tmb) == "status"] <- "tmb_status"
  names(rtmb)[names(rtmb) == "status"] <- "rtmb_status"

  ans <- merge(
    tmb[c("test_index", "file", "context", "test", "tmb_sec", "tmb_status")],
    rtmb[c("test_index", "file", "context", "test", "rtmb_sec",
           "rtmb_status")],
    by = c("test_index", "file", "context", "test"),
    all = TRUE
  )
  ans$extra_sec <- ans$rtmb_sec - ans$tmb_sec
  ans$ratio <- ans$rtmb_sec / ans$tmb_sec
  ans$status <- ifelse(ans$tmb_status == "OK" & ans$rtmb_status == "OK",
                       "OK", paste(ans$tmb_status, ans$rtmb_status, sep = "/"))
  ans
}

run_parent <- function() {
  if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    stop("Run this script from the repository root.")
  }

  tmb_file <- tempfile(fileext = ".rds")
  rtmb_file <- tempfile(fileext = ".rds")
  on.exit(unlink(c(tmb_file, rtmb_file)), add = TRUE)

  for (backend in c("TMB", "RTMB")) {
    cat("Timing ", backend, " tests", if (nzchar(filter)) {
      paste0(" matching filter '", filter, "'")
    } else {
      ""
    }, "...\n", sep = "")
    output <- if (backend == "TMB") tmb_file else rtmb_file
    status <- system2(
      "Rscript",
      c("--vanilla", script_path),
      env = c(
        paste0("RTMB_TEST_BENCHMARK_CHILD_BACKEND=", backend),
        paste0("RTMB_TEST_BENCHMARK_CHILD_OUTPUT=", output),
        paste0("RTMB_TEST_BENCHMARK_FILTER=", filter)
      ),
      stdout = FALSE,
      stderr = FALSE
    )
    if (!identical(status, 0L)) {
      stop(backend, " benchmark failed with exit status ", status)
    }
  }

  tmb <- readRDS(tmb_file)
  rtmb <- readRDS(rtmb_file)
  comparison <- compare_backends(tmb, rtmb)
  comparison <- comparison[order(-comparison$extra_sec), ]

  tmb_total <- sum(tmb$seconds)
  rtmb_total <- sum(rtmb$seconds)
  n_tests <- nrow(comparison)

  cat("\nRTMB vs TMB test-suite timing\n")
  cat("=============================\n")
  cat("Tests timed: ", n_tests, "\n", sep = "")
  cat("Filter: ", if (nzchar(filter)) filter else "<none>", "\n", sep = "")
  cat("TMB total:  ", format_num(tmb_total), " sec\n", sep = "")
  cat("RTMB total: ", format_num(rtmb_total), " sec\n", sep = "")
  cat("Overall RTMB/TMB ratio: ", format_num(rtmb_total / tmb_total, 2L),
      "x\n", sep = "")
  cat("TMB failures/errors:  ",
      sum(tmb$failures + tmb$errors), "\n", sep = "")
  cat("RTMB failures/errors: ",
      sum(rtmb$failures + rtmb$errors), "\n\n", sep = "")

  display <- head(comparison, top_n)
  display$rank <- seq_len(nrow(display))
  display$tmb_sec <- format_num(display$tmb_sec)
  display$rtmb_sec <- format_num(display$rtmb_sec)
  display$extra_sec <- format_num(display$extra_sec)
  display$ratio <- paste0(format_num(display$ratio, 2L), "x")

  cat("Highest-priority RTMB slowdowns, ranked by extra seconds\n")
  print_table(
    display,
    c("rank", "file", "test", "tmb_sec", "rtmb_sec", "extra_sec", "ratio",
      "status"),
    c(rank = 4L, file = 24L, test = 54L, tmb_sec = 8L, rtmb_sec = 8L,
      extra_sec = 9L, ratio = 7L, status = 9L)
  )

  ratio_candidates <- comparison[
    is.finite(comparison$ratio) & comparison$tmb_sec >= ratio_min_tmb,
  ]
  ratio_candidates <- ratio_candidates[order(-ratio_candidates$ratio), ]
  ratio_display <- head(ratio_candidates, min(top_n, nrow(ratio_candidates)))
  if (nrow(ratio_display) > 0L) {
    ratio_display$rank <- seq_len(nrow(ratio_display))
    ratio_display$tmb_sec <- format_num(ratio_display$tmb_sec)
    ratio_display$rtmb_sec <- format_num(ratio_display$rtmb_sec)
    ratio_display$extra_sec <- format_num(ratio_display$extra_sec)
    ratio_display$ratio <- paste0(format_num(ratio_display$ratio, 2L), "x")

    cat("\nLargest RTMB/TMB ratios, excluding tiny TMB tests below ",
        ratio_min_tmb, " sec\n", sep = "")
    print_table(
      ratio_display,
      c("rank", "file", "test", "tmb_sec", "rtmb_sec", "extra_sec",
        "ratio", "status"),
      c(rank = 4L, file = 24L, test = 54L, tmb_sec = 8L, rtmb_sec = 8L,
        extra_sec = 9L, ratio = 7L, status = 9L)
    )
  }

  invisible(comparison)
}

if (nzchar(Sys.getenv("RTMB_TEST_BENCHMARK_CHILD_BACKEND"))) {
  run_child()
} else {
  run_parent()
}
