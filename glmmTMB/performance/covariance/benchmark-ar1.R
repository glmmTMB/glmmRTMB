## Benchmark the AR(1) covariance structure at two block sizes.
## This follows issue #9 by timing the full glmmTMB() call for
## ar1(row + 0 | Subject). It compares installed-package TMB from a clean
## Rscript process against local-source RTMB after pkgload::load_all().

reps <- as.integer(Sys.getenv("RTMB_BENCHMARK_TIMES", "3"))
low_n <- as.integer(Sys.getenv("AR1_BENCHMARK_LOW_N", "10"))
high_n <- as.integer(Sys.getenv("AR1_BENCHMARK_HIGH_N", "180"))
if (is.na(reps) || reps < 1L) {
  stop("RTMB_BENCHMARK_TIMES must be a positive integer.")
}
if (is.na(low_n) || low_n < 2L) {
  stop("AR1_BENCHMARK_LOW_N must be at least 2.")
}
if (is.na(high_n) || high_n < 2L) {
  stop("AR1_BENCHMARK_HIGH_N must be at least 2.")
}

installed_tmb_times <- function(n, seed) {
  code <- sprintf(
    paste(
      "library(glmmTMB)",
      "reps <- %dL",
      "n <- %dL",
      "seed <- %dL",
      "make_ar1_data <- function(n, seed) {",
      "  set.seed(seed)",
      "  d <- data.frame(",
      "    Reaction = 300 + rnorm(n, sd = 30),",
      "    row = factor(seq_len(n), levels = seq_len(n), ordered = TRUE),",
      "    Subject = factor(rep(1, n))",
      "  )",
      "  d",
      "}",
      "d <- make_ar1_data(n, seed)",
      "fit <- function() glmmTMB(",
      "  Reaction ~ (1 | Subject) + ar1(row + 0 | Subject),",
      "  data = d",
      ")",
      "fit()",
      "times <- numeric(reps)",
      "for (i in seq_len(reps)) times[i] <- system.time(fit())[[\"elapsed\"]]",
      "cat(\"INSTALLED_TMB_PACKAGE=\", as.character(packageVersion(\"glmmTMB\")), \"\\n\", sep = \"\")",
      "cat(\"INSTALLED_TMB_TIMES=\", paste(times, collapse = \",\"), \"\\n\", sep = \"\")",
      sep = "\n"
    ),
    reps, n, seed
  )

  output <- system2(
    "Rscript",
    c("--vanilla", "-e", shQuote(code)),
    stdout = TRUE,
    stderr = TRUE
  )
  package_line <- grep("^INSTALLED_TMB_PACKAGE=", output, value = TRUE)
  times_line <- grep("^INSTALLED_TMB_TIMES=", output, value = TRUE)
  if (length(times_line) != 1L) {
    cat(paste(output, collapse = "\n"), "\n")
    stop("Could not parse installed-package TMB benchmark output.")
  }
  package_version <- sub("^INSTALLED_TMB_PACKAGE=", "", package_line[1L])
  times <- as.numeric(strsplit(sub("^INSTALLED_TMB_TIMES=", "", times_line),
                               ",", fixed = TRUE)[[1L]])
  attr(times, "package_version") <- package_version
  times
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package before running this benchmark.")
}

pkg_dir <- if (file.exists("DESCRIPTION")) "." else "glmmTMB"
pkgload::load_all(pkg_dir, quiet = TRUE)

old_use_rtmb <- glmmTMB::useRTMB()
on.exit(glmmTMB::useRTMB(old_use_rtmb), add = TRUE)

make_ar1_data <- function(n, seed) {
  set.seed(seed)
  data.frame(
    Reaction = 300 + rnorm(n, sd = 30),
    row = factor(seq_len(n), levels = seq_len(n), ordered = TRUE),
    Subject = factor(rep(1, n))
  )
}

benchmark_case <- function(d, label) {
  installed_tmb <- installed_tmb_times(nrow(d), attr(d, "seed"))

  fit_rtmb <- function() {
    glmmTMB::useRTMB(TRUE)
    stopifnot(identical(glmmTMB::useRTMB(), TRUE))
    glmmTMB(
      Reaction ~ (1 | Subject) + ar1(row + 0 | Subject),
      data = d
    )
  }

  ## Build the local RTMB objective before timing, as in issue #9.
  fit_rtmb()

  rtmb_times <- numeric(reps)
  for (i in seq_len(reps)) {
    rtmb_times[i] <- system.time(fit_rtmb())[["elapsed"]]
  }

  cat("\n", label, "\n", sep = "")
  print(data.frame(
    backend = rep(c("installed TMB", "local RTMB"), each = reps),
    run = rep(seq_len(reps), times = 2L),
    elapsed_seconds = c(installed_tmb, rtmb_times),
    row.names = NULL
  ), row.names = FALSE)

  result <- data.frame(
    case = label,
    n = nrow(d),
    installed_tmb_version = attr(installed_tmb, "package_version"),
    installed_tmb_mean_seconds = mean(installed_tmb),
    local_rtmb_mean_seconds = mean(rtmb_times),
    local_rtmb_vs_installed_tmb_ratio = mean(rtmb_times) / mean(installed_tmb),
    row.names = NULL
  )
  print(result, row.names = FALSE)
  result
}

low_rank_data <- make_ar1_data(low_n, seed = 1L)
attr(low_rank_data, "seed") <- 1L
high_rank_data <- make_ar1_data(high_n, seed = 2L)
attr(high_rank_data, "seed") <- 2L

results <- rbind(
  benchmark_case(low_rank_data, "low-dimensional AR1 block"),
  benchmark_case(high_rank_data, "high-dimensional AR1 block")
)

print(results, row.names = FALSE)
