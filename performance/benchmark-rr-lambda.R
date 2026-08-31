if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("Package 'microbenchmark' is required for this benchmark")
}

## Benchmark ways to fill the reduced-rank Lambda matrix.
##
## This checks whether the explicit nested loop in rtmb_tpl.R can be replaced
## by a simpler matrix-indexing assignment without changing the result.

fill_lambda_loop <- function(n, rr_rank, lam_diag, lam_lower) {
  Lambda <- matrix(0, n, rr_rank)
  lower_index <- 1L

  for (j in seq_len(rr_rank)) {
    for (i in seq_len(n)) {
      if (j == i) {
        Lambda[i, j] <- lam_diag[j]
      } else if (j < i) {
        Lambda[i, j] <- lam_lower[lower_index]
        lower_index <- lower_index + 1L
      }
    }
  }

  Lambda
}

fill_lambda_mask <- function(n, rr_rank, lam_diag, lam_lower) {
  Lambda <- matrix(0, n, rr_rank)
  Lambda[row(Lambda) == col(Lambda)] <- lam_diag
  Lambda[row(Lambda) > col(Lambda)] <- lam_lower
  Lambda
}

bench_lambda <- function(n = 20L, rr_rank = 5L, times = 1000L, run = NA_integer_) {
  stopifnot(rr_rank <= n)

  set.seed(1)
  lam_diag <- rnorm(rr_rank)
  lam_lower <- rnorm(n * rr_rank - rr_rank * (rr_rank + 1L) / 2L)

  loop_value <- fill_lambda_loop(n, rr_rank, lam_diag, lam_lower)
  mask_value <- fill_lambda_mask(n, rr_rank, lam_diag, lam_lower)

  stopifnot(identical(loop_value, mask_value))

  mb <- microbenchmark::microbenchmark(
    loop = fill_lambda_loop(n, rr_rank, lam_diag, lam_lower),
    mask = fill_lambda_mask(n, rr_rank, lam_diag, lam_lower),
    times = times
  )
  timings <- tapply(mb$time, mb$expr, mean) / 1e9

  data.frame(
    run = run,
    n = n,
    rr_rank = rr_rank,
    method = names(timings),
    seconds_per_eval = unname(timings),
    relative_to_loop = unname(timings / timings[["loop"]]),
    row.names = NULL
  )
}

bench_all <- function(run = NA_integer_) {
  rbind(
    bench_lambda(n = 5L, rr_rank = 2L, times = 5000L, run = run),
    bench_lambda(n = 20L, rr_rank = 5L, times = 2000L, run = run),
    bench_lambda(n = 100L, rr_rank = 10L, times = 500L, run = run),
    bench_lambda(n = 500L, rr_rank = 20L, times = 100L, run = run)
  )
}

bench_repeated <- function(iterations = 20L) {
  raw <- do.call(rbind, lapply(seq_len(iterations), bench_all))
  summary <- aggregate(
    cbind(seconds_per_eval, relative_to_loop) ~ n + rr_rank + method,
    data = raw,
    FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
  )

  summary <- do.call(data.frame, summary)
  names(summary) <- sub("\\.mean$", "_mean", names(summary))
  names(summary) <- sub("\\.sd$", "_sd", names(summary))

  list(raw = raw, summary = summary)
}

results <- bench_repeated(iterations = 20L)
print(results$summary)
