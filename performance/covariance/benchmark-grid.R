## Controlled TMB vs RTMB benchmark grid for families and covariance structures.
##
## Run from the repository root:
##   Rscript performance/covariance/benchmark-grid.R
##
## Useful smaller runs:
##   RTMB_GRID_FAMILIES=gaussian,poisson RTMB_GRID_COVSTRUCTS=cs,ou,exp \
##     Rscript performance/covariance/benchmark-grid.R
##
## Optional controls:
##   RTMB_GRID_FAMILIES=all
##   RTMB_GRID_COVSTRUCTS=all
##   RTMB_GRID_BLOCK_N=5,30
##   RTMB_GRID_GROUPS=4
##   RTMB_GRID_REPS=1
##   RTMB_GRID_TOP=30
##   RTMB_GRID_EVAL_MAX=50
##   RTMB_GRID_RATIO_MIN_TMB=0.05
##   RTMB_GRID_PROGRESS_EVERY=10

## Pin BLAS/OpenMP threading to 1 so that the parallel TMB/RTMB backend
## subprocesses spawned below (via parallel::mclapply) don't each try to
## claim every core, which oversubscribes the machine. This also applies to
## the child processes themselves, since they inherit this environment.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package before running this benchmark.")
}
if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  stop("Install the 'RhpcBLASctl' package before running this benchmark.")
}

## Sys.setenv(OPENBLAS_NUM_THREADS = ...) above only reaches this process's
## own threading if set before OpenBLAS's thread pool is first created,
## which for THIS process already happened at library-load time -- it does
## still correctly propagate to freshly spawned child processes, though.
## RhpcBLASctl calls the underlying C APIs directly and can resize the
## thread pool at any point, so each process (parent and children) also
## pins itself explicitly with it below.
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

script_path <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
  } else {
    normalizePath("performance/covariance/benchmark-grid.R", mustWork = TRUE)
  }
})

repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                           mustWork = TRUE)
pkg_dir <- file.path(repo_root, "glmmTMB")

all_families <- c(
  "gaussian", "poisson", "truncated_poisson", "binomial", "betabinomial",
  "combinomial", "beta", "ordbeta", "Gamma", "nbinom1", "nbinom2",
  "nbinom12", "truncated_nbinom1", "truncated_nbinom2", "genpois",
  "truncated_genpois", "compois", "truncated_compois", "t", "tweedie",
  "lognormal", "skewnormal", "bell"
)

all_covstructs <- c(
  "diag", "homdiag", "us", "cs", "homcs", "toep", "homtoep",
  "ar1", "hetar1", "ou", "exp", "gau", "mat", "rr", "propto", "equalto"
)

split_env <- function(var, default) {
  value <- Sys.getenv(var, default)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

parse_ints <- function(var, default) {
  out <- as.integer(split_env(var, default))
  if (any(is.na(out)) || any(out < 1L)) {
    stop(var, " must be a comma-separated list of positive integers.")
  }
  out
}

families <- split_env("RTMB_GRID_FAMILIES", "all")
if (identical(families, "all")) {
  families <- all_families
}
covstructs <- split_env("RTMB_GRID_COVSTRUCTS", "all")
if (identical(covstructs, "all")) {
  covstructs <- all_covstructs
}

bad_families <- setdiff(families, all_families)
bad_covstructs <- setdiff(covstructs, all_covstructs)
if (length(bad_families) > 0L) {
  stop("Unknown families: ", paste(bad_families, collapse = ", "))
}
if (length(bad_covstructs) > 0L) {
  stop("Unknown covariance structures: ", paste(bad_covstructs, collapse = ", "))
}

block_sizes <- parse_ints("RTMB_GRID_BLOCK_N", "5,30")
group_counts <- parse_ints("RTMB_GRID_GROUPS", "4")
reps <- as.integer(Sys.getenv("RTMB_GRID_REPS", "1"))
top_n <- as.integer(Sys.getenv("RTMB_GRID_TOP", "30"))
eval_max <- as.integer(Sys.getenv("RTMB_GRID_EVAL_MAX", "50"))
ratio_min_tmb <- as.numeric(Sys.getenv("RTMB_GRID_RATIO_MIN_TMB", "0.05"))
progress_every <- as.integer(Sys.getenv("RTMB_GRID_PROGRESS_EVERY", "10"))
if (any(is.na(c(reps, top_n, eval_max))) ||
    reps < 1L || top_n < 1L || eval_max < 1L) {
  stop("RTMB_GRID_REPS, RTMB_GRID_TOP, and RTMB_GRID_EVAL_MAX must be positive.")
}
if (is.na(ratio_min_tmb) || ratio_min_tmb < 0) {
  stop("RTMB_GRID_RATIO_MIN_TMB must be non-negative.")
}
if (is.na(progress_every) || progress_every < 0L) {
  stop("RTMB_GRID_PROGRESS_EVERY must be a non-negative integer.")
}

## Run the independent per-backend Rscript launches in parallel, using up to
## half of the available cores.
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

family_spec <- function(family_name) {
  switch(
    family_name,
    gaussian = list(family = stats::gaussian(), response = "y"),
    poisson = list(family = stats::poisson(), response = "y"),
    truncated_poisson = list(
      family = glmmTMB::truncated_poisson(),
      response = "y"
    ),
    binomial = list(family = stats::binomial(), response = "cbind(success, failure)"),
    betabinomial = list(
      family = glmmTMB::betabinomial(link = "logit"),
      response = "cbind(success, failure)"
    ),
    combinomial = list(
      family = glmmTMB::combinomial(),
      response = "cbind(success, failure)"
    ),
    beta = list(family = glmmTMB::beta_family(link = "logit"), response = "y"),
    ordbeta = list(
      family = glmmTMB::ordbeta(),
      response = "y",
      start = list(psi = c(-1, 1)),
      map = list(psi = factor(c(NA, NA)))
    ),
    Gamma = list(family = stats::Gamma(link = "log"), response = "y"),
    nbinom1 = list(family = glmmTMB::nbinom1(), response = "y"),
    nbinom2 = list(family = glmmTMB::nbinom2(), response = "y"),
    nbinom12 = list(
      family = glmmTMB::nbinom12(),
      response = "y",
      start = list(psi = 0.1),
      map = list(psi = factor(NA))
    ),
    truncated_nbinom1 = list(
      family = glmmTMB::truncated_nbinom1(),
      response = "y"
    ),
    truncated_nbinom2 = list(
      family = glmmTMB::truncated_nbinom2(),
      response = "y"
    ),
    genpois = list(family = glmmTMB::genpois(), response = "y"),
    truncated_genpois = list(
      family = glmmTMB::truncated_genpois(),
      response = "y"
    ),
    compois = list(family = glmmTMB::compois(), response = "y"),
    truncated_compois = list(
      family = glmmTMB::truncated_compois(),
      response = "y"
    ),
    t = list(
      family = glmmTMB::t_family(),
      response = "y",
      start = list(psi = log(6)),
      map = list(psi = factor(NA))
    ),
    tweedie = list(
      family = glmmTMB::tweedie(),
      response = "y",
      start = list(psi = stats::qlogis(0.5)),
      map = list(psi = factor(NA))
    ),
    lognormal = list(family = glmmTMB::lognormal(link = "log"), response = "y"),
    skewnormal = list(
      family = glmmTMB::skewnormal(),
      response = "y",
      start = list(psi = -2),
      map = list(psi = factor(NA))
    ),
    bell = list(family = glmmTMB::bell(), response = "y")
  )
}

merge_named <- function(x, y) {
  if (is.null(x)) {
    x <- list()
  }
  if (is.null(y)) {
    y <- list()
  }
  modifyList(x, y)
}

cov_theta <- function(covstruct, n) {
  corr_one <- 0.2
  logsd <- log(rep(0.35, n))
  switch(
    covstruct,
    diag = logsd,
    homdiag = log(0.35),
    us = c(logsd, rep(0.05, n * (n - 1L) / 2L)),
    cs = c(logsd, corr_one),
    homcs = c(log(0.35), corr_one),
    toep = c(logsd, rep(0.05, n - 1L)),
    homtoep = c(log(0.35), rep(0.05, n - 1L)),
    ar1 = c(log(0.35), corr_one),
    hetar1 = c(logsd, corr_one),
    ou = c(log(0.35), log(0.5)),
    exp = c(log(0.35), log(max(n / 3, 1))),
    gau = c(log(0.35), log(max(n / 3, 1))),
    mat = c(log(0.35), log(max(n / 3, 1)), log(1.2)),
    rr = {
      rank <- min(2L, n)
      c(rep(0.25, rank), rep(0.02, n * rank - rank * (rank + 1L) / 2L))
    },
    propto = c(logsd, rep(0.05, n * (n - 1L) / 2L), 0),
    equalto = c(logsd, rep(0.05, n * (n - 1L) / 2L))
  )
}

make_data <- function(family_name, covstruct, n, groups, seed) {
  set.seed(seed)
  level <- rep(seq_len(n), times = groups)
  group <- rep(seq_len(groups), each = n)
  total_n <- length(level)

  d <- data.frame(
    x = stats::rnorm(total_n),
    z = stats::rnorm(total_n),
    level = factor(level, levels = seq_len(n), ordered = TRUE),
    group = factor(group)
  )

  grid_side <- ceiling(sqrt(n))
  coords <- expand.grid(sx = seq_len(grid_side), sy = seq_len(grid_side))
  coords <- coords[seq_len(n), , drop = FALSE]
  d$sx <- coords$sx[level]
  d$sy <- coords$sy[level]
  d$pos <- glmmTMB::numFactor(d$sx, d$sy)

  time_values <- cumsum(stats::runif(n, min = 0.5, max = 1.5))
  d$time <- glmmTMB::numFactor(time_values[level])

  eta <- -0.25 + 0.25 * d$x
  mu <- exp(pmin(eta, 1.5))
  prob <- stats::plogis(eta)
  d$size <- 10L
  d$success <- stats::rbinom(total_n, size = d$size, prob = prob)
  d$failure <- d$size - d$success

  d$y <- switch(
    family_name,
    gaussian = eta + stats::rnorm(total_n, sd = 0.7),
    poisson = stats::rpois(total_n, lambda = mu),
    truncated_poisson = pmax(1L, stats::rpois(total_n, lambda = mu)),
    binomial = d$success,
    betabinomial = d$success,
    combinomial = d$success,
    beta = pmin(pmax(stats::plogis(eta + stats::rnorm(total_n, sd = 0.3)),
                    0.02), 0.98),
    ordbeta = {
      y <- pmin(pmax(stats::plogis(eta + stats::rnorm(total_n, sd = 0.5)),
                    0.02), 0.98)
      y[seq(1L, total_n, by = 11L)] <- 0
      y[seq(6L, total_n, by = 13L)] <- 1
      y
    },
    Gamma = stats::rgamma(total_n, shape = 3, scale = mu / 3),
    nbinom1 = stats::rnbinom(total_n, size = 3, mu = mu),
    nbinom2 = stats::rnbinom(total_n, size = 3, mu = mu),
    nbinom12 = stats::rnbinom(total_n, size = 3, mu = mu),
    truncated_nbinom1 = pmax(1L, stats::rnbinom(total_n, size = 3, mu = mu)),
    truncated_nbinom2 = pmax(1L, stats::rnbinom(total_n, size = 3, mu = mu)),
    genpois = stats::rpois(total_n, lambda = mu),
    truncated_genpois = pmax(1L, stats::rpois(total_n, lambda = mu)),
    compois = stats::rpois(total_n, lambda = mu),
    truncated_compois = pmax(1L, stats::rpois(total_n, lambda = mu)),
    t = eta + 0.7 * stats::rt(total_n, df = 6),
    tweedie = stats::rgamma(total_n, shape = 2, scale = mu / 2),
    lognormal = exp(eta + stats::rnorm(total_n, sd = 0.4)),
    skewnormal = eta + stats::rnorm(total_n, sd = 0.7),
    bell = stats::rpois(total_n, lambda = mu)
  )

  d
}

cov_term <- function(covstruct, n) {
  switch(
    covstruct,
    ou = "ou(0 + time | group)",
    exp = "exp(0 + pos | group)",
    gau = "gau(0 + pos | group)",
    mat = "mat(0 + pos | group)",
    rr = paste0("rr(0 + level | group, d = ", min(2L, n), ")"),
    propto = "propto(0 + level | group, fixed_covariance)",
    equalto = "equalto(0 + level | group, fixed_covariance)",
    paste0(covstruct, "(0 + level | group)")
  )
}

make_formula <- function(response, covstruct, n) {
  stats::as.formula(
    paste(response, "~ x +", cov_term(covstruct, n)),
    env = parent.frame()
  )
}

fixed_covariance_matrix <- function(n) {
  level_names <- paste0("level", seq_len(n))
  ans <- diag(n)
  dimnames(ans) <- list(level_names, level_names)
  ans
}

fit_case <- function(case, backend) {
  glmmTMB::useRTMB(backend == "RTMB")
  spec <- family_spec(case$family)
  data <- make_data(case$family, case$covstruct, case$block_n,
                    case$groups, case$seed)
  theta <- cov_theta(case$covstruct, case$block_n)
  start <- merge_named(spec$start, list(theta = theta))
  map <- merge_named(spec$map, list(theta = factor(rep(NA, length(theta)))))
  fixed_covariance <- fixed_covariance_matrix(case$block_n)
  formula <- make_formula(spec$response, case$covstruct, case$block_n)
  control <- glmmTMB::glmmTMBControl(
    optCtrl = list(iter.max = eval_max, eval.max = eval_max)
  )

  warnings <- character()
  fit <- withCallingHandlers(
    tryCatch(
      glmmTMB::glmmTMB(
        formula,
        family = spec$family,
        data = data,
        start = start,
        map = map,
        control = control,
        se = FALSE
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    ok = !inherits(fit, "error"),
    message = if (inherits(fit, "error")) conditionMessage(fit) else "",
    warnings = paste(unique(warnings), collapse = " | ")
  )
}

time_case <- function(case, backend) {
  elapsed <- system.time(result <- fit_case(case, backend))[["elapsed"]]
  data.frame(
    backend = backend,
    family = case$family,
    covstruct = case$covstruct,
    block_n = case$block_n,
    groups = case$groups,
    observations = case$block_n * case$groups,
    rep = case$rep,
    seconds = elapsed,
    status = if (result$ok) "OK" else "ERROR",
    message = result$message,
    warnings = result$warnings,
    stringsAsFactors = FALSE
  )
}

run_backend <- function(backend) {
  ## pkgload::load_all() defaults to debug = TRUE, which appends "-O0" after
  ## R's own "-O2" -- gcc/g++ take the last -O flag, so this silently builds
  ## an unoptimized binary. debug = FALSE gives a build comparable to a
  ## normal R CMD INSTALL (and to CRAN).
  pkgload::load_all(pkg_dir, quiet = TRUE, debug = FALSE)
  ## Re-pin after load_all(), in case compiling/loading the package reset
  ## the thread pool.
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  grid <- expand.grid(
    family = families,
    covstruct = covstructs,
    block_n = block_sizes,
    groups = group_counts,
    rep = seq_len(reps),
    stringsAsFactors = FALSE
  )
  grid$seed <- seq_len(nrow(grid)) + 9000L
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    if (progress_every > 0L &&
        (i == 1L || i == nrow(grid) || i %% progress_every == 0L)) {
      cat(backend, ": case ", i, " / ", nrow(grid), "\n", sep = "")
    }
    rows[[i]] <- time_case(grid[i, ], backend)
  }
  do.call(rbind, rows)
}

run_child <- function() {
  backend <- Sys.getenv("RTMB_GRID_CHILD_BACKEND")
  output <- Sys.getenv("RTMB_GRID_CHILD_OUTPUT")
  if (!backend %in% c("TMB", "RTMB")) {
    stop("Child backend must be TMB or RTMB.")
  }
  if (!nzchar(output)) {
    stop("Child output path is missing.")
  }
  saveRDS(run_backend(backend), output)
}

format_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

clip <- function(x, width) {
  x <- as.character(x)
  too_wide <- nchar(x, type = "width") > width
  x[too_wide] <- paste0(substr(x[too_wide], 1L, width - 1L), "~")
  x
}

print_table <- function(x, columns, widths) {
  if (nrow(x) == 0L) {
    cat("(none)\n")
    return(invisible(NULL))
  }
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

summarize_backend <- function(x) {
  aggregate(
    seconds ~ family + covstruct + block_n + groups + observations + backend,
    x,
    mean
  )
}

compare_backends <- function(tmb, rtmb) {
  tmb_sum <- summarize_backend(tmb)
  rtmb_sum <- summarize_backend(rtmb)
  names(tmb_sum)[names(tmb_sum) == "seconds"] <- "tmb_sec"
  names(rtmb_sum)[names(rtmb_sum) == "seconds"] <- "rtmb_sec"

  tmb_status <- aggregate(
    status ~ family + covstruct + block_n + groups + observations,
    tmb,
    function(z) if (all(z == "OK")) "OK" else "ERROR"
  )
  rtmb_status <- aggregate(
    status ~ family + covstruct + block_n + groups + observations,
    rtmb,
    function(z) if (all(z == "OK")) "OK" else "ERROR"
  )
  names(tmb_status)[names(tmb_status) == "status"] <- "tmb_status"
  names(rtmb_status)[names(rtmb_status) == "status"] <- "rtmb_status"

  out <- merge(
    tmb_sum,
    rtmb_sum,
    by = c("family", "covstruct", "block_n", "groups", "observations"),
    all = TRUE
  )
  out <- merge(out, tmb_status,
               by = c("family", "covstruct", "block_n", "groups",
                      "observations"), all.x = TRUE)
  out <- merge(out, rtmb_status,
               by = c("family", "covstruct", "block_n", "groups",
                      "observations"), all.x = TRUE)
  out$extra_sec <- out$rtmb_sec - out$tmb_sec
  out$ratio <- out$rtmb_sec / out$tmb_sec
  out$status <- ifelse(out$tmb_status == "OK" & out$rtmb_status == "OK",
                       "OK", paste(out$tmb_status, out$rtmb_status, sep = "/"))
  out
}

print_results <- function(comparison, raw) {
  ok <- comparison[comparison$status == "OK", ]
  slow_abs <- ok[order(-ok$extra_sec), ]
  slow_ratio <- ok[
    is.finite(ok$ratio) & ok$tmb_sec >= ratio_min_tmb,
  ]
  slow_ratio <- slow_ratio[order(-slow_ratio$ratio), ]
  errors <- comparison[comparison$status != "OK", ]

  cat("\nControlled RTMB covariance/family benchmark grid\n")
  cat("===============================================\n")
  cat("Families: ", paste(families, collapse = ", "), "\n", sep = "")
  cat("Covariance structures: ", paste(covstructs, collapse = ", "), "\n",
      sep = "")
  cat("Block sizes: ", paste(block_sizes, collapse = ", "), "\n", sep = "")
  cat("Groups: ", paste(group_counts, collapse = ", "), "\n", sep = "")
  cat("Repeats per backend/case: ", reps, "\n", sep = "")
  cat("Optimizer eval.max/iter.max: ", eval_max, "\n", sep = "")
  cat("Raw fits timed: ", nrow(raw), "\n", sep = "")
  cat("Comparable OK cases: ", nrow(ok), " / ", nrow(comparison), "\n\n",
      sep = "")

  display <- head(slow_abs, top_n)
  display$rank <- seq_len(nrow(display))
  display$tmb_sec <- format_num(display$tmb_sec)
  display$rtmb_sec <- format_num(display$rtmb_sec)
  display$extra_sec <- format_num(display$extra_sec)
  display$ratio <- paste0(format_num(display$ratio, 2L), "x")

  cat("Highest-priority slowdowns, ranked by extra RTMB seconds\n")
  print_table(
    display,
    c("rank", "family", "covstruct", "block_n", "groups", "observations",
      "tmb_sec", "rtmb_sec", "extra_sec", "ratio"),
    c(rank = 4L, family = 18L, covstruct = 10L, block_n = 7L, groups = 6L,
      observations = 12L, tmb_sec = 8L, rtmb_sec = 8L, extra_sec = 9L,
      ratio = 7L)
  )

  display <- head(slow_ratio, top_n)
  display$rank <- seq_len(nrow(display))
  display$tmb_sec <- format_num(display$tmb_sec)
  display$rtmb_sec <- format_num(display$rtmb_sec)
  display$extra_sec <- format_num(display$extra_sec)
  display$ratio <- paste0(format_num(display$ratio, 2L), "x")

  cat("\nLargest RTMB/TMB ratios, excluding TMB cases below ",
      ratio_min_tmb, " sec\n", sep = "")
  print_table(
    display,
    c("rank", "family", "covstruct", "block_n", "groups", "observations",
      "tmb_sec", "rtmb_sec", "extra_sec", "ratio"),
    c(rank = 4L, family = 18L, covstruct = 10L, block_n = 7L, groups = 6L,
      observations = 12L, tmb_sec = 8L, rtmb_sec = 8L, extra_sec = 9L,
      ratio = 7L)
  )

  if (nrow(errors) > 0L) {
    cat("\nCases that did not fit cleanly in both backends\n")
    error_messages <- aggregate(
      message ~ family + covstruct + block_n + groups + observations,
      raw[raw$status != "OK", ],
      function(z) paste(unique(z[nzchar(z)]), collapse = " | ")
    )
    errors <- merge(
      errors,
      error_messages,
      by = c("family", "covstruct", "block_n", "groups", "observations"),
      all.x = TRUE
    )
    print_table(
      errors,
      c("family", "covstruct", "block_n", "groups", "status", "message"),
      c(family = 18L, covstruct = 10L, block_n = 7L, groups = 6L,
        status = 12L, message = 52L)
    )
  }
}

run_parent <- function() {
  if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    stop("Run this script from the repository root.")
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install the 'pkgload' package before running this benchmark.")
  }
  ## Compile once, here in the parent, before the TMB and RTMB child
  ## processes are launched in parallel below -- otherwise both children
  ## could try to compile the same source concurrently and race.
  cat("Compiling local source (debug = FALSE) before timing...\n")
  pkgload::load_all(pkg_dir, quiet = TRUE, debug = FALSE)

  tmb_file <- tempfile(fileext = ".rds")
  rtmb_file <- tempfile(fileext = ".rds")
  on.exit(unlink(c(tmb_file, rtmb_file)), add = TRUE)

  backend_output <- list(TMB = tmb_file, RTMB = rtmb_file)

  run_backend_child <- function(backend) {
    output <- backend_output[[backend]]
    cat("Timing ", backend, " grid...\n", sep = "")
    status <- system2(
      "Rscript",
      c("--vanilla", script_path),
      env = c(
        paste0("RTMB_GRID_CHILD_BACKEND=", backend),
        paste0("RTMB_GRID_CHILD_OUTPUT=", output),
        paste0("RTMB_GRID_FAMILIES=", paste(families, collapse = ",")),
        paste0("RTMB_GRID_COVSTRUCTS=", paste(covstructs, collapse = ",")),
        paste0("RTMB_GRID_BLOCK_N=", paste(block_sizes, collapse = ",")),
        paste0("RTMB_GRID_GROUPS=", paste(group_counts, collapse = ",")),
        paste0("RTMB_GRID_REPS=", reps),
        paste0("RTMB_GRID_EVAL_MAX=", eval_max),
        paste0("RTMB_GRID_PROGRESS_EVERY=", progress_every)
      )
    )
    if (!identical(status, 0L)) {
      stop(backend, " benchmark failed with exit status ", status)
    }
    invisible(NULL)
  }

  run_parallel(c("TMB", "RTMB"), run_backend_child)

  tmb <- readRDS(tmb_file)
  rtmb <- readRDS(rtmb_file)
  raw <- rbind(tmb, rtmb)
  comparison <- compare_backends(tmb, rtmb)
  print_results(comparison, raw)

  results_file <- file.path(repo_root, "performance", "covariance",
                            "benchmark-grid-results.rds")
  saveRDS(list(comparison = comparison, raw = raw), results_file)
  cat("\nSaved results table to ", results_file, "\n", sep = "")

  invisible(comparison)
}

if (nzchar(Sys.getenv("RTMB_GRID_CHILD_BACKEND"))) {
  run_child()
} else {
  run_parent()
}
