## Benchmark covariance structures that use the unstructured-correlation path.
## This compares local-source TMB against local-source RTMB, both after
## pkgload::load_all().
##
## These structures scale quadratically in the block dimension, so the default
## high dimension is intentionally much smaller than the AR(1)/spatial examples.
##
## pkgload::load_all() defaults to debug = TRUE, which compiles with
## "-UNDEBUG -Wall -pedantic -g -O0" appended after R's own "-O2" -- since
## gcc/g++ take the last -O flag on the command line, this silently builds
## an unoptimized binary and makes local-TMB timings meaningless. We pass
## debug = FALSE below to get an optimized build comparable to a normal
## R CMD INSTALL (and to CRAN).

reps <- as.integer(Sys.getenv("RTMB_BENCHMARK_TIMES", "3"))
low_n <- as.integer(Sys.getenv("UNSTRUCT_BENCHMARK_LOW_N", "5"))
high_n <- as.integer(Sys.getenv("UNSTRUCT_BENCHMARK_HIGH_N", "20"))
if (is.na(reps) || reps < 1L) {
  stop("RTMB_BENCHMARK_TIMES must be a positive integer.")
}
if (is.na(low_n) || low_n < 2L) {
  stop("UNSTRUCT_BENCHMARK_LOW_N must be at least 2.")
}
if (is.na(high_n) || high_n < 2L) {
  stop("UNSTRUCT_BENCHMARK_HIGH_N must be at least 2.")
}

make_unstructured_data <- function(n, seed) {
  set.seed(seed)
  id <- factor(seq_len(n), levels = seq_len(n))
  data.frame(
    z = 1 + rnorm(n, sd = 1),
    id = id,
    group = factor(rep(1, n))
  )
}

make_covariance_matrix <- function(n) {
  nm <- paste0("id", seq_len(n))
  V <- diag(n)
  dimnames(V) <- list(nm, nm)
  V
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

benchmark_case <- function(structure, n, seed, label) {
  d <- make_unstructured_data(n, seed)
  V <- make_covariance_matrix(n)

  fit_case <- switch(
    structure,
    us = function() glmmTMB(z ~ 1 + us(0 + id | group), data = d),
    propto = function() glmmTMB(z ~ 1 + propto(0 + id | group, V), data = d),
    equalto = function() glmmTMB(z ~ 1 + equalto(0 + id | group, V), data = d),
    stop("Unknown covariance structure: ", structure)
  )

  fit_tmb <- function() {
    glmmTMB::useRTMB(FALSE)
    stopifnot(identical(glmmTMB::useRTMB(), FALSE))
    fit_case()
  }

  fit_rtmb <- function() {
    glmmTMB::useRTMB(TRUE)
    stopifnot(identical(glmmTMB::useRTMB(), TRUE))
    fit_case()
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
    structure = structure,
    backend = rep(c("local TMB", "local RTMB"), each = reps),
    run = rep(seq_len(reps), times = 2L),
    elapsed_seconds = c(tmb_times, rtmb_times),
    row.names = NULL
  ), row.names = FALSE)

  result <- data.frame(
    structure = structure,
    case = label,
    n = n,
    local_tmb_mean_seconds = mean(tmb_times),
    local_rtmb_mean_seconds = mean(rtmb_times),
    local_rtmb_vs_local_tmb_ratio = mean(rtmb_times) / mean(tmb_times),
    row.names = NULL
  )
  print(result, row.names = FALSE)
  result
}

structures <- c("us", "propto", "equalto")
results <- do.call(
  rbind,
  lapply(structures, function(structure) {
    rbind(
      benchmark_case(
        structure, low_n, seed = 1L,
        label = paste("low-dimensional", structure, "block")
      ),
      benchmark_case(
        structure, high_n, seed = 2L,
        label = paste("high-dimensional", structure, "block")
      )
    )
  })
)

print(results, row.names = FALSE)

results_file <- "performance/covariance/benchmark-unstructured-results.rds"
saveRDS(results, results_file)
cat("\nSaved results table to ", results_file, "\n", sep = "")
