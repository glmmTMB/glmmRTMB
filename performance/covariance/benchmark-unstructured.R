## Benchmark covariance structures that use the unstructured-correlation path.
## This compares installed-package TMB from a clean Rscript process against
## local-source RTMB after pkgload::load_all().
##
## These structures scale quadratically in the block dimension, so the default
## high dimension is intentionally much smaller than the AR(1)/spatial examples.

## Pin BLAS/OpenMP threading to 1 so that parallel workers (spawned below via
## parallel::mclapply, and the installed-TMB Rscript subprocesses they call)
## don't each try to claim every core, which oversubscribes the machine.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")

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

## Run the independent installed-TMB Rscript launches in parallel, using up
## to half of the available cores.
n_cores <- max(1L, round(parallel::detectCores() / 2))
run_parallel <- function(X, FUN) {
  if (.Platform$OS.type == "windows" || n_cores <= 1L) {
    return(lapply(X, FUN))
  }
  result <- parallel::mclapply(X, FUN, mc.cores = min(length(X), n_cores))
  failed <- vapply(result, function(x) inherits(x, "try-error"), logical(1L))
  if (any(failed)) {
    messages <- vapply(result[failed], function(x) {
      cond <- attr(x, "condition")
      if (is.null(cond)) as.character(x) else conditionMessage(cond)
    }, character(1L))
    stop("Parallel benchmark run failed: ", paste(messages, collapse = " | "))
  }
  result
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

fit_expression <- function(structure, n) {
  if (structure == "us") {
    "glmmTMB(z ~ 1 + us(0 + id | group), data = d)"
  } else if (structure == "propto") {
    "glmmTMB(z ~ 1 + propto(0 + id | group, V), data = d)"
  } else if (structure == "equalto") {
    "glmmTMB(z ~ 1 + equalto(0 + id | group, V), data = d)"
  } else {
    stop("Unknown covariance structure: ", structure)
  }
}

installed_tmb_times <- function(structure, n, seed) {
  code <- sprintf(
    paste(
      "library(glmmTMB)",
      "reps <- %dL",
      "n <- %dL",
      "seed <- %dL",
      "make_unstructured_data <- function(n, seed) {",
      "  set.seed(seed)",
      "  id <- factor(seq_len(n), levels = seq_len(n))",
      "  data.frame(",
      "    z = 1 + rnorm(n, sd = 1),",
      "    id = id,",
      "    group = factor(rep(1, n))",
      "  )",
      "}",
      "make_covariance_matrix <- function(n) {",
      "  nm <- paste0(\"id\", seq_len(n))",
      "  V <- diag(n)",
      "  dimnames(V) <- list(nm, nm)",
      "  V",
      "}",
      "d <- make_unstructured_data(n, seed)",
      "V <- make_covariance_matrix(n)",
      "fit <- function() %s",
      "fit()",
      "times <- numeric(reps)",
      "for (i in seq_len(reps)) times[i] <- system.time(fit())[[\"elapsed\"]]",
      "cat(\"INSTALLED_TMB_PACKAGE=\", as.character(packageVersion(\"glmmTMB\")), \"\\n\", sep = \"\")",
      "cat(\"INSTALLED_TMB_TIMES=\", paste(times, collapse = \",\"), \"\\n\", sep = \"\")",
      sep = "\n"
    ),
    reps, n, seed, fit_expression(structure, n)
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

pkgload::load_all("glmmTMB", quiet = TRUE)

old_use_rtmb <- glmmTMB::useRTMB()
on.exit(glmmTMB::useRTMB(old_use_rtmb), add = TRUE)

benchmark_case <- function(structure, n, seed, label, installed_tmb) {
  d <- make_unstructured_data(n, seed)
  V <- make_covariance_matrix(n)

  fit_rtmb <- switch(
    structure,
    us = function() {
      glmmTMB::useRTMB(TRUE)
      stopifnot(identical(glmmTMB::useRTMB(), TRUE))
      glmmTMB(z ~ 1 + us(0 + id | group), data = d)
    },
    propto = function() {
      glmmTMB::useRTMB(TRUE)
      stopifnot(identical(glmmTMB::useRTMB(), TRUE))
      glmmTMB(z ~ 1 + propto(0 + id | group, V), data = d)
    },
    equalto = function() {
      glmmTMB::useRTMB(TRUE)
      stopifnot(identical(glmmTMB::useRTMB(), TRUE))
      glmmTMB(z ~ 1 + equalto(0 + id | group, V), data = d)
    },
    stop("Unknown covariance structure: ", structure)
  )

  ## Build the local RTMB objective before timing, as in issue #9.
  fit_rtmb()

  rtmb_times <- numeric(reps)
  for (i in seq_len(reps)) {
    rtmb_times[i] <- system.time(fit_rtmb())[["elapsed"]]
  }

  cat("\n", label, "\n", sep = "")
  print(data.frame(
    structure = structure,
    backend = rep(c("installed TMB", "local RTMB"), each = reps),
    run = rep(seq_len(reps), times = 2L),
    elapsed_seconds = c(installed_tmb, rtmb_times),
    row.names = NULL
  ), row.names = FALSE)

  result <- data.frame(
    structure = structure,
    case = label,
    n = n,
    installed_tmb_version = attr(installed_tmb, "package_version"),
    installed_tmb_mean_seconds = mean(installed_tmb),
    local_rtmb_mean_seconds = mean(rtmb_times),
    local_rtmb_vs_installed_tmb_ratio = mean(rtmb_times) / mean(installed_tmb),
    row.names = NULL
  )
  print(result, row.names = FALSE)
  result
}

structures <- c("us", "propto", "equalto")

cases <- do.call(c, lapply(structures, function(structure) {
  list(
    list(structure = structure, n = low_n, seed = 1L,
         label = paste("low-dimensional", structure, "block")),
    list(structure = structure, n = high_n, seed = 2L,
         label = paste("high-dimensional", structure, "block"))
  )
}))

installed_all <- run_parallel(cases, function(case) {
  installed_tmb_times(case$structure, case$n, case$seed)
})

results <- do.call(rbind, Map(
  function(case, installed_tmb) {
    benchmark_case(case$structure, case$n, case$seed, case$label,
                   installed_tmb)
  },
  cases, installed_all
))

print(results, row.names = FALSE)

results_file <- "performance/covariance/benchmark-unstructured-results.rds"
saveRDS(results, results_file)
cat("\nSaved results table to ", results_file, "\n", sep = "")
