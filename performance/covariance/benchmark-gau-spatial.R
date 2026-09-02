## Benchmark the spatial Gaussian covariance structure at two block sizes.
## This follows the exp benchmark and times the full glmmTMB() call for
## gau(pos + 0 | group), comparing local TMB against local RTMB (both
## after pkgload::load_all()).
##
## pkgload::load_all() defaults to debug = TRUE, which compiles with
## "-UNDEBUG -Wall -pedantic -g -O0" appended after R's own "-O2" -- since
## gcc/g++ take the last -O flag on the command line, this silently builds
## an unoptimized binary and makes local-TMB timings meaningless. We pass
## debug = FALSE below to get an optimized build comparable to a normal
## R CMD INSTALL (and to CRAN).

reps <- as.integer(Sys.getenv("RTMB_BENCHMARK_TIMES", "3"))
high_n <- as.integer(Sys.getenv("GAU_BENCHMARK_HIGH_N", "300"))
if (is.na(reps) || reps < 1L) stop("RTMB_BENCHMARK_TIMES must be positive.")
if (is.na(high_n) || high_n < 2L) stop("GAU_BENCHMARK_HIGH_N must be at least 2.")

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

make_data <- function(n, seed) {
  set.seed(seed)
  if (n == 100L) {
    d <- data.frame(z = as.vector(volcano), x = as.vector(row(volcano)),
                    y = as.vector(col(volcano)))
    d$z <- d$z + rnorm(length(d$z), sd = 15)
    d <- d[sample(nrow(d), 100), ]
  } else {
    side <- ceiling(sqrt(n)); xy <- expand.grid(x = seq_len(side), y = seq_len(side))
    d <- xy[seq_len(n), ]; d$z <- rnorm(n)
  }
  d$pos <- glmmTMB::numFactor(d$x, d$y)
  d$group <- factor(rep(1, nrow(d)))
  d
}

benchmark_case <- function(d, label) {
  fit_tmb <- function() {
    glmmTMB::useRTMB(FALSE); stopifnot(!glmmTMB::useRTMB())
    glmmTMB(z ~ 1 + gau(pos + 0 | group), data = d)
  }
  fit_rtmb <- function() {
    glmmTMB::useRTMB(TRUE); stopifnot(glmmTMB::useRTMB())
    glmmTMB(z ~ 1 + gau(pos + 0 | group), data = d)
  }
  fit_tmb(); fit_rtmb()
  tmb <- rtmb <- numeric(reps)
  for (i in seq_len(reps)) {
    tmb[i] <- system.time(fit_tmb())[["elapsed"]]
    rtmb[i] <- system.time(fit_rtmb())[["elapsed"]]
  }
  print(data.frame(case = label,
                   backend = rep(c("local TMB", "local RTMB"), each = reps),
                   run = rep(seq_len(reps), times = 2L),
                   elapsed_seconds = c(tmb, rtmb)), row.names = FALSE)
  data.frame(case = label, n = nrow(d),
             local_tmb_mean_seconds = mean(tmb),
             local_rtmb_mean_seconds = mean(rtmb),
             local_rtmb_vs_local_tmb_ratio = mean(rtmb) / mean(tmb))
}

low <- make_data(100L, 1L)
high <- make_data(high_n, 2L)

results <- rbind(
  benchmark_case(low, "low-dimensional volcano example"),
  benchmark_case(high, "high-dimensional spatial block")
)

print(results, row.names = FALSE)

results_file <- "performance/covariance/benchmark-gau-results.rds"
saveRDS(results, results_file)
cat("\nSaved results table to ", results_file, "\n", sep = "")
