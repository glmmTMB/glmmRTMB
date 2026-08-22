## Benchmark the spatial Gaussian covariance structure at two block sizes.
## This follows the exp benchmark and times the full glmmTMB() call for
## gau(pos + 0 | group) using installed TMB, local TMB, and local RTMB.

reps <- as.integer(Sys.getenv("RTMB_BENCHMARK_TIMES", "3"))
high_n <- as.integer(Sys.getenv("GAU_BENCHMARK_HIGH_N", "300"))
if (is.na(reps) || reps < 1L) stop("RTMB_BENCHMARK_TIMES must be positive.")
if (is.na(high_n) || high_n < 2L) stop("GAU_BENCHMARK_HIGH_N must be at least 2.")

installed_tmb_times <- function(n, seed) {
  code <- sprintf(paste(
    "library(glmmTMB)", "reps <- %dL", "n <- %dL", "seed <- %dL",
    "make_data <- function(n, seed) {",
    "  set.seed(seed)",
    "  if (n == 100L) {",
    "    d <- data.frame(z = as.vector(volcano), x = as.vector(row(volcano)),",
    "                    y = as.vector(col(volcano)))",
    "    d$z <- d$z + rnorm(length(d$z), sd = 15)",
    "    d <- d[sample(nrow(d), 100), ]",
    "  } else {",
    "    side <- ceiling(sqrt(n)); xy <- expand.grid(x = seq_len(side), y = seq_len(side))",
    "    d <- xy[seq_len(n), ]; d$z <- rnorm(n)",
    "  }",
    "  d$pos <- glmmTMB::numFactor(d$x, d$y); d$group <- factor(rep(1, nrow(d))); d",
    "}",
    "d <- make_data(n, seed)",
    "fit <- function() glmmTMB(z ~ 1 + gau(pos + 0 | group), data = d)",
    "fit(); times <- numeric(reps)",
    "for (i in seq_len(reps)) times[i] <- system.time(fit())[[\"elapsed\"]]",
    "cat(\"INSTALLED_GAU_PACKAGE=\", as.character(packageVersion(\"glmmTMB\")), \"\\n\", sep=\"\")",
    "cat(\"INSTALLED_GAU_TIMES=\", paste(times, collapse=\",\"), \"\\n\", sep=\"\")",
    sep = "\n"), reps, n, seed)
  output <- system2("Rscript", c("--vanilla", "-e", shQuote(code)),
                    stdout = TRUE, stderr = TRUE)
  version <- grep("^INSTALLED_GAU_PACKAGE=", output, value = TRUE)
  times <- grep("^INSTALLED_GAU_TIMES=", output, value = TRUE)
  if (length(times) != 1L) stop("Could not parse installed TMB benchmark output.")
  ans <- as.numeric(strsplit(sub("^INSTALLED_GAU_TIMES=", "", times), ",")[[1L]])
  attr(ans, "package_version") <- sub("^INSTALLED_GAU_PACKAGE=", "", version[1L])
  ans
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package before running this benchmark.")
}
pkg_dir <- if (file.exists("DESCRIPTION")) "." else "glmmTMB"
pkgload::load_all(pkg_dir, quiet = TRUE)
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

benchmark_case <- function(d, label, seed) {
  installed <- installed_tmb_times(nrow(d), seed)
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
                   backend = rep(c("installed TMB", "local TMB", "local RTMB"),
                                 each = reps),
                   run = rep(seq_len(reps), times = 3L),
                   elapsed_seconds = c(installed, tmb, rtmb)), row.names = FALSE)
  data.frame(case = label, n = nrow(d),
             installed_tmb_version = attr(installed, "package_version"),
             installed_tmb_mean_seconds = mean(installed),
             local_tmb_mean_seconds = mean(tmb),
             local_rtmb_mean_seconds = mean(rtmb),
             local_rtmb_vs_installed_tmb_ratio = mean(rtmb) / mean(installed),
             local_rtmb_vs_local_tmb_ratio = mean(rtmb) / mean(tmb))
}

low <- make_data(100L, 1L)
high <- make_data(high_n, 2L)
results <- rbind(benchmark_case(low, "low-dimensional volcano example", 1L),
                 benchmark_case(high, "high-dimensional spatial block", 2L))
print(results, row.names = FALSE)
