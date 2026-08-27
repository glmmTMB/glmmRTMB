## Benchmark the Ornstein-Uhlenbeck covariance structure at two block sizes.
## This follows issue #9 by timing the full glmmTMB() call for
## ou(time + 0 | group). It reports three baselines:
##   1. installed-package TMB, from a clean Rscript process;
##   2. local-source TMB, after pkgload::load_all();
##   3. local-source RTMB, after pkgload::load_all().

reps <- as.integer(Sys.getenv("RTMB_BENCHMARK_TIMES", "3"))
high_n <- as.integer(Sys.getenv("OU_BENCHMARK_HIGH_N", "300"))
if (is.na(reps) || reps < 1L) {
  stop("RTMB_BENCHMARK_TIMES must be a positive integer.")
}
if (is.na(high_n) || high_n < 2L) {
  stop("OU_BENCHMARK_HIGH_N must be at least 2.")
}

installed_tmb_times <- function(n, seed) {
  code <- sprintf(
    paste(
      "library(glmmTMB)",
      "reps <- %dL",
      "n <- %dL",
      "seed <- %dL",
      "make_ou_data <- function(n, seed) {",
      "  set.seed(seed)",
      "  times <- sort(cumsum(runif(n, min = 0.5, max = 1.5)))",
      "  d <- data.frame(",
      "    z = 1 + 0.02 * times + rnorm(n, sd = 1),",
      "    time = glmmTMB::numFactor(times),",
      "    group = factor(rep(1, n))",
      "  )",
      "  d",
      "}",
      "d <- make_ou_data(n, seed)",
      "fit <- function() glmmTMB(z ~ 1 + ou(time + 0 | group), data = d)",
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

make_ou_data <- function(n, seed) {
  set.seed(seed)
  times <- sort(cumsum(runif(n, min = 0.5, max = 1.5)))
  data.frame(
    z = 1 + 0.02 * times + rnorm(n, sd = 1),
    time = glmmTMB::numFactor(times),
    group = factor(rep(1, n))
  )
}

benchmark_case <- function(d, label) {
  installed_tmb <- installed_tmb_times(nrow(d), attr(d, "seed"))

  fit_tmb <- function() {
    glmmTMB::useRTMB(FALSE)
    stopifnot(identical(glmmTMB::useRTMB(), FALSE))
    glmmTMB(z ~ 1 + ou(time + 0 | group), data = d)
  }

  fit_rtmb <- function() {
    glmmTMB::useRTMB(TRUE)
    stopifnot(identical(glmmTMB::useRTMB(), TRUE))
    glmmTMB(z ~ 1 + ou(time + 0 | group), data = d)
  }

  ## Compile/build both objective functions before timing, as in issue #9.
  fit_tmb()
  fit_rtmb()

  tmb_times <- numeric(reps)
  rtmb_times <- numeric(reps)
  for (i in seq_len(reps)) {
    tmb_times[i] <- system.time(fit_tmb())[["elapsed"]]
    rtmb_times[i] <- system.time(fit_rtmb())[["elapsed"]]
  }

  cat("\n", label, "\n", sep = "")
  print(data.frame(
    backend = rep(c("installed TMB", "local TMB", "local RTMB"), each = reps),
    run = rep(seq_len(reps), times = 3L),
    elapsed_seconds = c(installed_tmb, tmb_times, rtmb_times),
    row.names = NULL
  ), row.names = FALSE)

  result <- data.frame(
    case = label,
    n = nrow(d),
    installed_tmb_version = attr(installed_tmb, "package_version"),
    installed_tmb_mean_seconds = mean(installed_tmb),
    local_tmb_mean_seconds = mean(tmb_times),
    local_rtmb_mean_seconds = mean(rtmb_times),
    local_rtmb_vs_installed_tmb_ratio = mean(rtmb_times) / mean(installed_tmb),
    local_rtmb_vs_local_tmb_ratio = mean(rtmb_times) / mean(tmb_times),
    row.names = NULL
  )
  print(result, row.names = FALSE)
  result
}

low_rank_data <- make_ou_data(100L, seed = 1L)
attr(low_rank_data, "seed") <- 1L
high_rank_data <- make_ou_data(high_n, seed = 2L)
attr(high_rank_data, "seed") <- 2L

results <- rbind(
  benchmark_case(low_rank_data, "low-dimensional OU block"),
  benchmark_case(high_rank_data, "high-dimensional OU block")
)

print(results, row.names = FALSE)
