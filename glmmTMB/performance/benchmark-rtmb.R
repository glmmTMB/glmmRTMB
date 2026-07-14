## Compare TMB and RTMB fitting performance for the same model.
## Run from the repository root or the glmmTMB package directory:
##   Rscript glmmTMB/benchmark-rtmb.R
##   Rscript benchmark-rtmb.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package before running this benchmark.")
}
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("Install the 'microbenchmark' package before running this benchmark.")
}

pkg_dir <- if (file.exists("DESCRIPTION")) {
  "."
} else if (file.exists(file.path("glmmTMB", "DESCRIPTION"))) {
  "glmmTMB"
} else {
  stop("Run this script from the repository root or glmmTMB package directory.")
}

pkgload::load_all(pkg_dir, quiet = TRUE)

old_use_rtmb <- glmmTMB:::useRTMB()
on.exit(glmmTMB:::useRTMB(old_use_rtmb), add = TRUE)

times <- as.integer(Sys.getenv("RTMB_BENCHMARK_TIMES", "100"))
if (is.na(times) || times < 1L) {
  stop("RTMB_BENCHMARK_TIMES must be a positive integer.")
}

benchmark_fit <- function(label, fit) {
  fit_tmb <- function() {
    glmmTMB:::useRTMB(FALSE)
    fit()
    invisible(NULL)
  }

  fit_rtmb <- function() {
    glmmTMB:::useRTMB(TRUE)
    fit()
    invisible(NULL)
  }

  ## Construct each backend's objective before collecting timings.
  fit_tmb()
  fit_rtmb()
  invisible(gc())

  set.seed(101)
  benchmark <- microbenchmark::microbenchmark(
    TMB = fit_tmb(),
    RTMB = fit_rtmb(),
    times = times,
    unit = "ms",
    control = list(order = "random")
  )

  timings <- summary(benchmark)
  tmb_mean <- timings$mean[timings$expr == "TMB"]
  rtmb_mean <- timings$mean[timings$expr == "RTMB"]
  tmb_median <- timings$median[timings$expr == "TMB"]
  rtmb_median <- timings$median[timings$expr == "RTMB"]

  cat("\n", label, "\n", sep = "")
  print(
    timings[, c("expr", "min", "lq", "mean", "median", "uq", "max")],
    row.names = FALSE
  )
  cat(
    sprintf(
      "RTMB/TMB ratio: %.3f (mean), %.3f (median)\n",
      rtmb_mean / tmb_mean,
      rtmb_median / tmb_median
    )
  )

  invisible(benchmark)
}

fit_poisson <- function() {
  glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
    data = Salamanders
  )
}

data("sleepstudy", package = "lme4")

fit_gaussian <- function() {
  glmmTMB(
    Reaction ~ Days + (Days | Subject),
    family = gaussian,
    data = sleepstudy
  )
}

benchmark_fit("Poisson: Salamanders random intercept", fit_poisson)
benchmark_fit("Gaussian: sleepstudy correlated random slope", fit_gaussian)
