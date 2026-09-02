## Benchmark the Toeplitz covariance structure at two block sizes.
## This follows issue #9 by timing the full glmmTMB() call for
## toep(time + 0 | group). It compares local-source TMB against
## local-source RTMB, both after pkgload::load_all().
##
## pkgload::load_all() defaults to debug = TRUE, which compiles with
## "-UNDEBUG -Wall -pedantic -g -O0" appended after R's own "-O2" -- since
## gcc/g++ take the last -O flag on the command line, this silently builds
## an unoptimized binary and makes local-TMB timings meaningless. We pass
## debug = FALSE below to get an optimized build comparable to a normal
## R CMD INSTALL (and to CRAN).

reps <- as.integer(Sys.getenv("RTMB_BENCHMARK_TIMES", "3"))
low_n <- as.integer(Sys.getenv("TOEP_BENCHMARK_LOW_N", "100"))
high_n <- as.integer(Sys.getenv("TOEP_BENCHMARK_HIGH_N", "300"))
if (is.na(reps) || reps < 1L) {
  stop("RTMB_BENCHMARK_TIMES must be a positive integer.")
}
if (is.na(low_n) || low_n < 2L) {
  stop("TOEP_BENCHMARK_LOW_N must be at least 2.")
}
if (is.na(high_n) || high_n < 2L) {
  stop("TOEP_BENCHMARK_HIGH_N must be at least 2.")
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package before running this benchmark.")
}
if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  stop("Install the 'RhpcBLASctl' package before running this benchmark.")
}

## Pin BLAS/OpenMP threading to 1. Sys.setenv(OPENBLAS_NUM_THREADS = ...)
## does NOT work for this process's own threading -- OpenBLAS reads that
## env var once, at library-load time, before any of our own code runs.
## RhpcBLASctl calls the underlying C APIs directly, which can resize the
## thread pool at any point, including here.
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

pkgload::load_all("glmmTMB", quiet = TRUE, debug = FALSE)

## Re-pin after load_all(), in case compiling/loading the package reset
## the thread pool.
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

old_use_rtmb <- glmmTMB::useRTMB()
on.exit(glmmTMB::useRTMB(old_use_rtmb), add = TRUE)

make_toep_data <- function(n, seed) {
  set.seed(seed)
  time <- factor(seq_len(n), levels = seq_len(n), ordered = TRUE)
  data.frame(
    z = 1 + 0.02 * seq_len(n) + rnorm(n, sd = 1),
    time = time,
    group = factor(rep(1, n))
  )
}

benchmark_case <- function(d, label) {
  fit_tmb <- function() {
    glmmTMB::useRTMB(FALSE)
    stopifnot(identical(glmmTMB::useRTMB(), FALSE))
    glmmTMB(z ~ 1 + toep(time + 0 | group), data = d)
  }

  fit_rtmb <- function() {
    glmmTMB::useRTMB(TRUE)
    stopifnot(identical(glmmTMB::useRTMB(), TRUE))
    glmmTMB(z ~ 1 + toep(time + 0 | group), data = d)
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
    backend = rep(c("local TMB", "local RTMB"), each = reps),
    run = rep(seq_len(reps), times = 2L),
    elapsed_seconds = c(tmb_times, rtmb_times),
    row.names = NULL
  ), row.names = FALSE)

  result <- data.frame(
    case = label,
    n = nrow(d),
    local_tmb_mean_seconds = mean(tmb_times),
    local_rtmb_mean_seconds = mean(rtmb_times),
    local_rtmb_vs_local_tmb_ratio = mean(rtmb_times) / mean(tmb_times),
    row.names = NULL
  )
  print(result, row.names = FALSE)
  result
}

low_rank_data <- make_toep_data(low_n, seed = 1L)
high_rank_data <- make_toep_data(high_n, seed = 2L)

results <- rbind(
  benchmark_case(low_rank_data, "low-dimensional Toeplitz block"),
  benchmark_case(high_rank_data, "high-dimensional Toeplitz block")
)

print(results, row.names = FALSE)

results_file <- "performance/covariance/benchmark-toep-results.rds"
saveRDS(results, results_file)
cat("\nSaved results table to ", results_file, "\n", sep = "")
