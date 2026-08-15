cmb <- function(f, d) function(p) f(p, d)

osa_keep <- function(x) {
  if (inherits(x, "osa")) {
    as.vector(x@keep[, 1L])
  } else {
    rep(1, length(x))
  }
}

osa_value <- function(x) {
  if (inherits(x, "osa")) x@x else x
}

## Translated from logit_inverse_linkfun(), glmmTMB.cpp:213-232.
## The binomial likelihood uses logit(probability), not necessarily eta.
logit_inverse_linkfun_rtmb <- function(eta, link) {
  switch(
    names(link),
    logit = eta,
    probit = RTMB::pnorm(eta, log.p = TRUE) -
      RTMB::pnorm(eta, lower.tail = FALSE, log.p = TRUE),
    cloglog = RTMB::logspace_sub(exp(eta), 0),
    {
      mu <- switch(
        names(link),
        log = exp(eta),
        identity = eta,
        sqrt = eta * eta,
        inverse = 1 / eta,
        lambertW = exp(eta) * exp(exp(eta)),
        stop("link not yet implemented for binomial: ", names(link))
      )
      log(mu) - log(1 - mu)
    }
  )
}

## Translated from log_inverse_linkfun(), glmmTMB.cpp:234-249.
## Negative-binomial likelihoods use log(mu) for robust density evaluation.
log_inverse_linkfun_rtmb <- function(eta, link) {
  switch(
    names(link),
    log = eta,
    logit = -RTMB::logspace_add(0, -eta),
    {
      mu <- switch(
        names(link),
        probit = RTMB::pnorm(eta),
        cloglog = 1 - exp(-exp(eta)),
        identity = eta,
        sqrt = eta * eta,
        inverse = 1 / eta,
        lambertW = exp(eta) * exp(exp(eta)),
        stop("link not yet implemented for log inverse-link: ", names(link))
      )
      log(mu)
    }
  )
}

logit_mu_rtmb <- function(eta, link) {
  logit_inverse_linkfun_rtmb(eta, link)
}

log_mu_rtmb <- function(eta, link) {
  log_inverse_linkfun_rtmb(eta, link)
}

log_var_minus_mu_rtmb <- function(family_name, log_mu, etadisp, psi) {
  switch(
    family_name,
    nbinom1 = log_mu + etadisp,
    truncated_nbinom1 = log_mu + etadisp,
    nbinom2 = 2 * log_mu - etadisp,
    truncated_nbinom2 = 2 * log_mu - etadisp,
    nbinom12 = {
      log_mu_vec <- log_mu[seq_along(log_mu)]
      etadisp_vec <- etadisp[seq_along(etadisp)]
      log_mu_vec + RTMB::logspace_add(etadisp_vec, log_mu_vec - psi[1L])
    },
    stop("log(var - mu) not defined for distribution: ", family_name)
  )
}

#' Simulate from a zero-inflated density wrapper
#'
#' This helper implements the simulation branch used by [dZI()].  The
#' zero-inflation predictor `eta_zi` is on the logit scale, so the structural
#' zero probability is `p_zi = 1 / (1 + exp(-eta_zi))` and the complementary
#' conditional probability is `1 - p_zi = 1 / (1 + exp(eta_zi))`.  Simulated
#' values therefore come from the same mixture represented in the likelihood:
#' either the structural-zero component returns zero, or the wrapped
#' conditional density simulates the response.  Keeping this logic outside
#' [dZI()] makes the likelihood code easier to read and preserves RTMB's
#' convention that simulation is triggered by evaluating log-density functions
#' on `simref` objects.
#'
#' @param density A log-density function such as `RTMB::dpois()` or
#'   `RTMB::dnorm()` that also supports RTMB simulation through `simref`
#'   objects.
#' @param x The response vector, normally an RTMB `simref` object during
#'   simulation.
#' @param ... Distribution-specific arguments passed to `density`.
#' @param eta_zi Zero-inflation linear predictor on the logit scale.
#'
#' @return A zero vector on the log-density scale, after mutating the `simref`
#'   response object with simulated values.
#'
#' @noRd
simZI <- function(density, x, ..., eta_zi) {
  prob_nonzero <- 1 / (1 + exp(eta_zi))

  nonzero <- as.logical(stats::rbinom(length(x), 1, prob_nonzero))
  density_args <- list(...)
  density_args <- lapply(density_args, function(arg) {
    if (length(arg) == length(x)) arg[nonzero] else arg
  })

  if (any(nonzero)) {
    do.call(
      density,
      c(list(x = x[nonzero]), density_args, list(log = TRUE))
    )
  }
  if (any(!nonzero)) {
    structural_zero <- x[!nonzero]
    structural_zero[] <- 0
  }
  rep(0, length(x))
}

#' Add zero inflation to a density function
#'
#' `dZI()` takes an ordinary density function and
#' returns a new density function with glmmTMB-style zero-inflation behavior.
#' When `eta_zi` is `NULL`, the wrapper deliberately reduces to the original
#' density, so non-zero-inflated models use the same likelihood path.  When
#' `eta_zi` is supplied, the likelihood is the standard zero-inflated mixture.
#' A nonzero observation can only come from the conditional distribution, so it
#' contributes `log(1 - p_zi) + log f(y)`.  An observed zero can come either
#' from the structural-zero component or from the conditional distribution, so
#' it contributes `log(p_zi + (1 - p_zi) * f(0))`.  The wrapper evaluates these
#' terms on the log scale, using `RTMB::logspace_add()` for the zero case so
#' the two possible zero sources are combined stably.  The separate `is_zero`
#' argument tells the wrapper which observations should use the zero-mixture
#' formula; this is especially important for truncated or hurdle-like families
#' where the conditional density should not be evaluated at zero.
#'
#' @param density A log-density function to wrap.
#'
#' @return A function with the same distribution-specific arguments as
#'   `density`, plus `eta_zi`, `log`, and `is_zero`.
#'
#' @noRd
dZI <- function(density) {
  force(density)

  function(x, ..., eta_zi = NULL, log = FALSE, is_zero = NULL) {
    if (inherits(x, "simref") && !is.null(eta_zi)) {
      return(simZI(density, x, ..., eta_zi = eta_zi))
    }

    x <- osa_value(x)
    if (is.null(eta_zi)) {
      loglik <- density(x, ..., log = TRUE)
      return(if (log) loglik else exp(loglik))
    }

    if (is.null(is_zero)) {
      is_zero <- x == 0
    }
    has_zero <- any(is_zero)

    loglik <- density(x, ..., log = TRUE)
    log_1mpz <- -RTMB::logspace_add(0, eta_zi)
    ans <- log_1mpz + loglik

    if (has_zero) {
      log_pz <- -RTMB::logspace_add(0, -eta_zi)
      ans[is_zero] <- RTMB::logspace_add(
        log_pz[is_zero],
        ans[is_zero]
      )
    }
    if (log) ans else exp(ans)
  }
}

## Fitting uses RTMB::dbinom_robust() to match glmmTMB.cpp:981.
## Simulation follows glmmTMB.cpp:982 but delegates to RTMB::dbinom(),
## because RTMB::dbinom_robust() does not provide rbinom_robust().
dbinom_robust_rtmb <- function(x, size, logit_p, log = FALSE) {
  if (inherits(x, "simref")) {
    prob <- 1 / (1 + exp(-logit_p))
    return(RTMB::dbinom(x, size = size, prob = prob, log = log))
  }
  RTMB::dbinom_robust(x, size = size, logit_p = logit_p, log = log)
}

## Translated from glmmtmb::logspace_gamma(), distrib.h:27-47.
## This avoids lgamma(exp(x)) underflow for very small x, matching the C++
## stable beta-binomial density helper.
logspace_gamma_rtmb <- RTMB::Vectorize(
  function(x) {
    `if` <- RTMB::ADoverload("if")
    if (x < -150) -x else lgamma(exp(x))
  },
  vectorize.args = "x"
)

## Translated from glmmtmb::dbetabinom_robust(), distrib.h:49-62, and
## the betabinomial_family case in glmmTMB.cpp:1037-1047.
dbetabinom_robust_rtmb <- function(x, log_shape1, log_shape2, size,
                                   log = FALSE) {
  if (inherits(x, "simref")) {
    log_shape1 <- rep(as.vector(log_shape1), length.out = length(x))
    log_shape2 <- rep(as.vector(log_shape2), length.out = length(x))
    size <- rep(as.vector(size), length.out = length(x))
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      prob <- stats::rbeta(1L, exp(log_shape1[i]), exp(log_shape2[i]))
      ans[i] <- stats::rbinom(1L, size = size[i], prob = prob)
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  log_x <- log(x)
  log_size_minus_x <- log(size - x)
  ans <-
    lgamma(size + 1) - lgamma(x + 1) - lgamma(size - x + 1) +
    logspace_gamma_rtmb(RTMB::logspace_add(log_x, log_shape1)) +
    logspace_gamma_rtmb(RTMB::logspace_add(log_size_minus_x, log_shape2)) -
    lgamma(size + exp(log_shape1) + exp(log_shape2)) +
    lgamma(exp(log_shape1) + exp(log_shape2)) -
    logspace_gamma_rtmb(log_shape1) - logspace_gamma_rtmb(log_shape2)
  if (log) ans else exp(ans)
}

## Mean-parameterized Conway-Maxwell-Binomial density; translated from
## dcombinom2() and combinom_utils in TMB's distributions_R.hpp:660-695 and
## tiny_ad/compois/combinom.hpp:1-96.
combinom_logZ_rtmb <- function(logit_p, nu, size) {
  ans <- -Inf
  for (k in 0:size) {
    ans <- RTMB::logspace_add(
      ans,
      nu * lchoose(size, k) + k * logit_p
    )
  }
  ans
}

combinom_moments_rtmb <- function(logit_p, nu, size) {
  log_z <- combinom_logZ_rtmb(logit_p, nu, size)
  mean <- second_moment <- 0
  for (k in 0:size) {
    probability <- exp(
      nu * lchoose(size, k) + k * logit_p - log_z
    )
    mean <- mean + k * probability
    second_moment <- second_moment + k * k * probability
  }
  list(mean = mean, variance = second_moment - mean * mean)
}

combinom_logitp_rtmb <- function(mean, nu, size) {
  "if" <- RTMB::ADoverload("if")
  lower <- -30
  upper <- 30
  for (i in seq_len(10L)) {
    midpoint <- (lower + upper) / 2
    midpoint_mean <- combinom_moments_rtmb(midpoint, nu, size)$mean
    new_lower <- if (midpoint_mean > mean) lower else midpoint
    new_upper <- if (midpoint_mean > mean) midpoint else upper
    lower <- new_lower
    upper <- new_upper
  }

  logit_p <- (lower + upper) / 2
  for (i in seq_len(5L)) {
    moments <- combinom_moments_rtmb(logit_p, nu, size)
    logit_p <- logit_p - (moments$mean - mean) / moments$variance
  }
  logit_p
}

rcombinom2_rtmb <- function(mean, nu, size) {
  logit_p <- combinom_logitp_rtmb(mean, nu, size)
  log_weights <- vapply(
    0:size,
    function(k) nu * lchoose(size, k) + k * logit_p,
    numeric(1)
  )
  weights <- exp(log_weights - max(log_weights))
  random_number <- stats::runif(1L) * sum(weights)
  which(cumsum(weights) >= random_number)[1L] - 1L
}

dcombinom2_rtmb <- function(x, size, mean, nu, log = FALSE) {
  if (inherits(x, "simref")) {
    mean <- as.vector(mean)
    nu <- as.vector(nu)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rcombinom2_rtmb(mean[i], nu[i], size[i])
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  "[<-" <- RTMB::ADoverload("[<-")
  ans <- mean * 0
  for (i in seq_along(x)) {
    logit_p <- combinom_logitp_rtmb(mean[i], nu[i], size[i])
    ans[i] <- nu[i] * lchoose(size[i], x[i]) + x[i] * logit_p -
      combinom_logZ_rtmb(logit_p, nu[i], size[i])
  }
  if (log) ans else exp(ans)
}

## Translated from the Gamma_family case in glmmTMB.cpp:991-996 and
## zt_lik_zero(), glmmTMB.cpp:959. Gamma is mean/shape parameterized here:
## shape = phi and scale = mu / phi. Exact zeros are excluded from the
## conditional density so zero-inflated Gamma behaves as a hurdle model.
dgamma_rtmb <- function(x, mean, shape, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- stats::rgamma(
      length(x),
      shape = as.vector(shape),
      scale = as.vector(mean / shape)
    )
    return(rep(0, length(x)))
  }

  is_zero <- x == 0
  if (!any(is_zero)) {
    ans <- RTMB::dgamma(x, shape = shape, scale = mean / shape, log = TRUE)
  } else {
    not_zero <- !is_zero
    if (!any(not_zero)) {
      return(rep(if (log) -Inf else 0, length(x)))
    }

    "[<-" <- RTMB::ADoverload("[<-")
    subset_arg <- function(arg) {
      if (length(arg) == length(x)) arg[not_zero] else arg
    }
    mean_nz <- subset_arg(mean)
    shape_nz <- subset_arg(shape)
    loglik_nz <- RTMB::dgamma(
      x[not_zero],
      shape = shape_nz,
      scale = mean_nz / shape_nz,
      log = TRUE
    )
    ans <- rep(loglik_nz[1L] * 0 - Inf, length(x))
    ans[not_zero] <- loglik_nz
  }
  if (log) ans else exp(ans)
}

## Translated from the beta_family case in glmmTMB.cpp:997-1002 and
## zt_lik_zero(), glmmTMB.cpp:959. The beta family uses Ferrari-Cribari-Neto
## parameterization: shape1 = mu * phi and shape2 = (1 - mu) * phi. Exact zeros
## are excluded from the conditional density so zero-inflated beta behaves as a
## hurdle model.
dbeta_rtmb <- function(x, mean, phi, log = FALSE) {
  shape1 <- mean * phi
  shape2 <- (1 - mean) * phi

  if (inherits(x, "simref")) {
    x[] <- stats::rbeta(
      length(x),
      shape1 = as.vector(shape1),
      shape2 = as.vector(shape2)
    )
    return(rep(0, length(x)))
  }

  is_zero <- x == 0
  if (!any(is_zero)) {
    ans <- RTMB::dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE)
  } else {
    not_zero <- !is_zero
    if (!any(not_zero)) {
      return(rep(if (log) -Inf else 0, length(x)))
    }

    "[<-" <- RTMB::ADoverload("[<-")
    subset_arg <- function(arg) {
      if (length(arg) == length(x)) arg[not_zero] else arg
    }
    shape1_nz <- subset_arg(shape1)
    shape2_nz <- subset_arg(shape2)
    loglik_nz <- RTMB::dbeta(
      x[not_zero],
      shape1 = shape1_nz,
      shape2 = shape2_nz,
      log = TRUE
    )
    ans <- rep(loglik_nz[1L] * 0 - Inf, length(x))
    ans[not_zero] <- loglik_nz
  }
  if (log) ans else exp(ans)
}

## Translated from the ordbeta_family case in glmmTMB.cpp:1004-1031.
dordbeta_rtmb <- function(x, eta, mean, phi, cutpoints, log = FALSE) {
  if (inherits(x, "simref")) {
    eta <- as.vector(eta)
    mean <- as.vector(mean)
    phi <- as.vector(phi)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      if (stats::runif(1L) < stats::plogis(cutpoints[1L] - eta[i])) {
        ans[i] <- 0
      } else if (stats::runif(1L) < stats::plogis(eta[i] - cutpoints[2L])) {
        ans[i] <- 1
      } else {
        ans[i] <- stats::rbeta(
          1L,
          shape1 = mean[i] * phi[i],
          shape2 = (1 - mean[i]) * phi[i]
        )
      }
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  "[<-" <- RTMB::ADoverload("[<-")
  is_zero <- x == 0
  is_one <- x == 1
  is_interior <- !is_zero & !is_one
  ans <- eta * 0

  if (any(is_zero)) {
    ans[is_zero] <- -RTMB::logspace_add(
      0,
      eta[is_zero] - cutpoints[1L]
    )
  }
  if (any(is_one)) {
    ans[is_one] <- -RTMB::logspace_add(
      0,
      cutpoints[2L] - eta[is_one]
    )
  }
  if (any(is_interior)) {
    log_lower <- -RTMB::logspace_add(
      0,
      cutpoints[1L] - eta[is_interior]
    )
    log_upper <- -RTMB::logspace_add(
      0,
      cutpoints[2L] - eta[is_interior]
    )
    log_middle <- RTMB::logspace_sub(log_lower, log_upper)
    ans[is_interior] <- log_middle + RTMB::dbeta(
      x[is_interior],
      shape1 = mean[is_interior] * phi[is_interior],
      shape2 = (1 - mean[is_interior]) * phi[is_interior],
      log = TRUE
    )
  }
  if (log) ans else exp(ans)
}

## Translated from the lognormal_family case in glmmTMB.cpp:1164-1179.
## The lognormal family is parameterized by mean and SD on the data scale.
dlognormal_rtmb <- function(x, mean, sd, log = FALSE) {
  log_var <- RTMB::logspace_add(2 * (log(sd) - log(mean)), 0)
  meanlog <- log(mean) - log_var / 2
  sdlog <- sqrt(log_var)

  if (inherits(x, "simref")) {
    x[] <- stats::rlnorm(length(x), meanlog = as.vector(meanlog),
                         sdlog = as.vector(sdlog))
    return(rep(0, length(x)))
  }

  is_zero <- x == 0
  if (!any(is_zero)) {
    ans <- RTMB::dnorm(log(x), meanlog, sdlog, log = TRUE) - log(x)
  } else {
    not_zero <- !is_zero
    if (!any(not_zero)) {
      return(rep(if (log) -Inf else 0, length(x)))
    }

    "[<-" <- RTMB::ADoverload("[<-")
    subset_arg <- function(arg) {
      if (length(arg) == length(x)) arg[not_zero] else arg
    }
    meanlog_nz <- subset_arg(meanlog)
    sdlog_nz <- subset_arg(sdlog)
    loglik_nz <- RTMB::dnorm(
      log(x[not_zero]),
      meanlog_nz,
      sdlog_nz,
      log = TRUE
    ) - log(x[not_zero])
    ans <- rep(loglik_nz[1L] * 0 - Inf, length(x))
    ans[not_zero] <- loglik_nz
  }
  if (log) ans else exp(ans)
}

## Translated from the t_family case in glmmTMB.cpp:1182-1190.
## The response is standardized by the fitted scale phi, so the log-density
## subtracts log(phi), represented by etadisp in the C++ code.
dt_rtmb <- function(x, mean, scale, df, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- as.vector(mean) + as.vector(scale) * stats::rt(length(x), df)
    return(rep(0, length(x)))
  }

  ans <- RTMB::dt((x - mean) / scale, df = df, log = TRUE) - log(scale)
  if (log) ans else exp(ans)
}

## Translated from glmmtmb::dskewnorm(), distrib.h:78-88, and the
## skewnormal_family case in glmmTMB.cpp:975-980.
dskewnormal_rtmb <- function(x, mean, sd, alpha, log = FALSE) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sd / sqrt(1 - 2 / pi * delta^2)
  xi <- mean - omega * delta * sqrt(2 / pi)

  if (inherits(x, "simref")) {
    n <- length(x)
    xi <- rep(as.vector(xi), length.out = n)
    omega <- rep(as.vector(omega), length.out = n)
    delta <- rep(as.vector(delta), length.out = n)
    ans <- numeric(n)
    for (i in seq_len(n)) {
      chi <- abs(stats::rnorm(1L))
      nrv <- stats::rnorm(1L)
      z <- delta[i] * chi + sqrt(1 - delta[i]^2) * nrv
      ans[i] <- xi[i] + omega[i] * z
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  z <- (x - xi) / omega
  ans <- log(2) - log(omega) +
    RTMB::dnorm(z, 0, 1, log = TRUE) +
    log(RTMB::pnorm(alpha * z))
  if (log) ans else exp(ans)
}

## Fitting translates glmmTMB.cpp:1042-1075 for nbinom1/nbinom2:
## both families use dnbinom_robust(log_mu, log_var_minus_mu). Simulation
## follows the same mean/variance by converting back to size/mu.
dnbinom_robust_rtmb <- function(x, log_mu, log_var_minus_mu, log = FALSE) {
  if (inherits(x, "simref")) {
    mu <- exp(as.vector(log_mu))
    size <- exp(as.vector(2 * log_mu - log_var_minus_mu))
    x[] <- stats::rnbinom(length(x), size = size, mu = mu)
    return(rep(0, length(x)))
  }
  RTMB::dnbinom_robust(x, log_mu = log_mu, log_var_minus_mu = log_var_minus_mu,
                       log = log)
}

## Translated from glmmtmb::dgenpois(); distrib.h:64-78.
dgenpois_log_rtmb <- RTMB::Vectorize(
  function(x, theta, lambda) {
    term <- theta + lambda * x
    log(theta) + (x - 1) * log(term) - term - lgamma(x + 1)
  },
  vectorize.args = c("x", "theta", "lambda")
)

## Translated from glmmtmb::rgenpois(); distrib.h:172-183.
rgenpois_rtmb <- function(theta, lambda) {
  ans <- 0
  random_number <- stats::runif(1L)
  kum <- exp(dgenpois_log_rtmb(0, theta, lambda))
  while (random_number > kum) {
    ans <- ans + 1
    kum <- kum + exp(dgenpois_log_rtmb(ans, theta, lambda))
  }
  ans
}

## Translated from glmmtmb::rtruncated_genpois(); distrib.h:186-197.
rtruncated_genpois_rtmb <- function(theta, lambda) {
  nloop <- 10000L
  counter <- 0L
  ans <- rgenpois_rtmb(theta, lambda)
  while (ans < 1 && counter < nloop) {
    ans <- rgenpois_rtmb(theta, lambda)
    counter <- counter + 1L
  }
  if (ans < 1) {
    warning(
      "Zeros in simulation of zero-truncated data. ",
      "Possibly due to low estimated mean.",
      call. = FALSE
    )
  }
  ans
}

## Translated from the genpois_family case in glmmTMB.cpp:1128-1133.
dgenpois_rtmb <- function(x, theta, lambda, log = FALSE) {
  if (inherits(x, "simref")) {
    theta <- rep(as.vector(theta), length.out = length(x))
    lambda <- rep(as.vector(lambda), length.out = length(x))
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rgenpois_rtmb(theta[i], lambda[i])
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  ans <- dgenpois_log_rtmb(x, theta, lambda)
  if (log) ans else exp(ans)
}

## AD-compatible Lambert W transformation used by the Bell family; translated
## from glmmtmb::LambertW(), distrib.h:486-521.
lambertW_rtmb <- RTMB::Vectorize(
  function(x) {
    y <- log1p(x)
    for (i in seq_len(12L)) {
      y <- y - (y - x * exp(-y)) / (1 + y)
    }
    y
  },
  vectorize.args = "x"
)

## Translation of glmmtmb::Bell(), distrib.h:442-464.
bell_number_rtmb <- function(n) {
  if (n < 2L) {
    return(1)
  }

  bell <- bell_new <- numeric(n)
  bell[1L] <- 1
  for (i in seq_len(n - 1L)) {
    bell_new[1L] <- bell[i]
    for (j in seq_len(i)) {
      bell_new[j + 1L] <- bell[j] + bell_new[j]
    }
    bell <- bell_new
  }
  bell_new[n]
}

## Translation of glmmtmb::rbell() and glmmtmb::dbell(),
## distrib.h:417-432 and 466-478.
dbell_rtmb <- function(x, mean, log = FALSE) {
  theta <- lambertW_rtmb(mean)
  if (inherits(x, "simref")) {
    theta <- as.vector(theta)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      n_compound <- stats::rpois(1L, expm1(theta[i]))
      if (n_compound > 0L) {
        ans[i] <- sum(vapply(
          seq_len(n_compound),
          function(j) rtruncated_poisson_rtmb(theta[i]),
          numeric(1)
        ))
      }
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  bell_number <- vapply(as.integer(x), bell_number_rtmb, numeric(1))
  ans <- x * log(theta) - exp(theta) + 1 +
    log(bell_number) - lgamma(x + 1)
  if (log) ans else exp(ans)
}

## Translated from calc_log_nzprob(), glmmTMB.cpp:286-291.
log_nzprob_truncated_genpois_rtmb <- RTMB::Vectorize(
  function(theta) RTMB::logspace_sub(0, -theta),
  vectorize.args = "theta"
)

## Translated from the truncated_genpois_family case,
## glmmTMB.cpp:1134-1139.
dtruncated_genpois_rtmb <- function(x, theta, lambda, log = FALSE) {
  if (inherits(x, "simref")) {
    theta <- rep(as.vector(theta), length.out = length(x))
    lambda <- rep(as.vector(lambda), length.out = length(x))
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rtruncated_genpois_rtmb(theta[i], lambda[i])
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  log_nzprob <- log_nzprob_truncated_genpois_rtmb(theta)
  ans <- dgenpois_log_rtmb(x, theta, lambda) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## Translated from the compois_family case in glmmTMB.cpp:1115-1119.
## Fitting uses RTMB::dcompois2(mean, nu); simulation uses RTMB's internal
## rcompois2() because the exported density is the fitting interface.
dcompois2_rtmb <- function(x, mean, nu, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- get("rcompois2", envir = asNamespace("RTMB"))(
      length(x),
      mean = as.vector(mean),
      nu = as.vector(nu)
    )
    return(rep(0, length(x)))
  }
  RTMB::dcompois2(x, mean = mean, nu = nu, log = log)
}

## Translation of glmmtmb::rtruncated_compois2(); distrib.h:280-289.
rtruncated_compois2_rtmb <- function(n, mean, nu) {
  rcompois2 <- get("rcompois2", envir = asNamespace("RTMB"))
  mean <- rep(mean, length.out = n)
  nu <- rep(nu, length.out = n)
  ans <- rcompois2(n, mean = mean, nu = nu)

  nloop <- 10000L
  counter <- 0L
  while (any(ans < 1) && counter < nloop) {
    zero <- ans < 1
    ans[zero] <- rcompois2(sum(zero), mean = mean[zero], nu = nu[zero])
    counter <- counter + 1L
  }
  if (any(ans < 1)) {
    warning(
      "Zeros in simulation of zero-truncated data. ",
      "Possibly due to low estimated mean.",
      call. = FALSE
    )
  }
  ans
}

## translation of glmmtmb::rtruncated_nbinom(); distrib.h:130-168
rtruncated_nbinom_rtmb <- function(n, size, k = 0L, mu) {
  ans <- numeric(n)
  size <- rep(size, length.out = n)
  mu <- rep(mu, length.out = n)

  for (i in seq_len(n)) {
    if (size[i] <= 0) {
      stop("non-positive size in k-truncated-neg-bin simulator")
    }
    if (mu[i] <= 0) {
      stop("non-positive mu in k-truncated-neg-bin simulator")
    }
    if (k < 0) {
      stop("negative k in k-truncated-neg-bin simulator")
    }

    p <- size[i] / (mu[i] + size[i])
    q <- mu[i] / (mu[i] + size[i])
    m <- ceiling(max((k + 1) * p - size[i] * q, 0))

    repeat {
      x <- stats::rnbinom(1L, size = size[i] + m, prob = p) + m
      if (m > 0) {
        a <- 1
        u <- stats::runif(1L)
        for (j in seq_len(m)) {
          a <- a * (k + 2 - j) / (x - j + 1)
        }
        if (u < a && x > k) {
          break
        }
      } else if (x > k) {
        break
      }
    }
    ans[i] <- x
  }
  ans
}

## Scalar formulas vectorized to mirror the observation loop in
## calc_log_nzprob(), glmmTMB.cpp:270-290
log_nzprob_truncated_poisson_rtmb <- RTMB::Vectorize(
  function(mu) RTMB::logspace_sub(0, -mu),
  vectorize.args = "mu"
)

log_nzprob_truncated_nbinom1_rtmb <- RTMB::Vectorize(
  function(mu, log_phi) {
    log_phi_plus_one <- RTMB::logspace_add(0, log_phi)
    RTMB::logspace_sub(0, -mu / exp(log_phi) * log_phi_plus_one)
  },
  vectorize.args = c("mu", "log_phi")
)

log_nzprob_truncated_nbinom2_rtmb <- RTMB::Vectorize(
  function(log_mu, log_size) {
    log_ratio_plus_one <- RTMB::logspace_add(0, log_mu - log_size)
    RTMB::logspace_sub(0, -exp(log_size) * log_ratio_plus_one)
  },
  vectorize.args = c("log_mu", "log_size")
)

log_nzprob_truncated_compois_rtmb <- RTMB::Vectorize(
  function(mean, nu) {
    RTMB::logspace_sub(
      0,
      RTMB::dcompois2(0, mean = mean, nu = nu, log = TRUE)
    )
  },
  vectorize.args = c("mean", "nu")
)

## Translation of glmmtmb::rtruncated_poisson(); distrib.h:94-128.
rtruncated_poisson_rtmb <- function(mu, k = 0L) {
  if (mu <= 0) {
    stop("non-positive mu in k-truncated-poisson simulator")
  }
  if (k < 0) {
    stop("negative k in k-truncated-poisson simulator")
  }

  mdoub <- max(k + 1 - mu, 0)
  m <- ceiling(mdoub)

  repeat {
    x <- stats::rpois(1L, mu) + m
    if (m > 0) {
      a <- 1
      u <- stats::runif(1L)
      for (j in seq_len(m) - 1L) {
        a <- a * (k + 1 - j) / (x - j)
      }
      if (u < a && x > k) {
        return(x)
      }
    } else if (x > k) {
      return(x)
    }
  }
}

## Translated from calc_log_nzprob() and truncated_compois_family,
## glmmTMB.cpp:293-294 and 1121-1127.
dtruncated_compois2_rtmb <- function(x, mean, nu, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- rtruncated_compois2_rtmb(
      length(x),
      mean = as.vector(mean),
      nu = as.vector(nu)
    )
    return(rep(0, length(x)))
  }

  log_nzprob <- log_nzprob_truncated_compois_rtmb(mean, nu)
  ans <- RTMB::dcompois2(x, mean = mean, nu = nu, log = TRUE) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## Translated from calc_log_nzprob() and the truncated_nbinom1_family
## likelihood case, glmmTMB.cpp:274-277 and 1042-1064
dtruncated_nbinom1_rtmb <- function(x, log_mu, log_var_minus_mu, log_phi,
                                    log = FALSE) {
  if (inherits(x, "simref")) {
    sim_mu <- exp(as.vector(log_mu))
    sim_phi <- exp(as.vector(log_phi))
    x[] <- rtruncated_nbinom_rtmb(
      length(x),
      size = sim_mu / sim_phi,
      k = 0L,
      mu = sim_mu
    )
    return(rep(0, length(x)))
  }

  mu <- exp(log_mu)
  log_nzprob <- log_nzprob_truncated_nbinom1_rtmb(mu, log_phi)
  ans <- RTMB::dnbinom_robust(x, log_mu = log_mu,
                              log_var_minus_mu = log_var_minus_mu,
                              log = TRUE) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## Translated from calc_log_nzprob() and the truncated_nbinom2_family
## likelihood case, glmmTMB.cpp:278-283 and 1066-1081
dtruncated_nbinom2_rtmb <- function(x, log_mu, log_var_minus_mu, log_size,
                                    log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- rtruncated_nbinom_rtmb(
      length(x),
      size = exp(as.vector(log_size)),
      k = 0L,
      mu = exp(as.vector(log_mu))
    )
    return(rep(0, length(x)))
  }

  log_nzprob <- log_nzprob_truncated_nbinom2_rtmb(log_mu, log_size)
  ans <- RTMB::dnbinom_robust(x, log_mu = log_mu,
                              log_var_minus_mu = log_var_minus_mu,
                              log = TRUE) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## zero-truncated poisson density
dtruncated_poisson_rtmb <- function(x, lambda, log = FALSE) {
  if (inherits(x, "simref")) {
    lambda <- rep(as.vector(lambda), length.out = length(x))
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rtruncated_poisson_rtmb(lambda[i], k = 0L)
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  log_nzprob <- RTMB::logspace_sub(0, -lambda)
  ans <- RTMB::dpois(x, lambda = lambda, log = TRUE) - log_nzprob

  ## the conditional distribution has strictly positive support
  is_zero <- x < 0.001
  if (any(is_zero)) {
    ## return -Inf to let dZI() treat observed zeros as structural
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

apply_zi_prediction <- function(mu, eta, etazi, ziPredictCode) {
  if (ziPredictCode == .valid_zipredictcode[["corrected"]]) {
    pz <- 1 / (1 + exp(-etazi))
    mu <- mu * (1 - pz)
  } else if (ziPredictCode == .valid_zipredictcode[["uncorrected"]]) {
    ## leave mu and eta unchanged
  } else if (ziPredictCode == .valid_zipredictcode[["prob"]]) {
    mu <- 1 / (1 + exp(-etazi))
    eta <- etazi
  } else if (ziPredictCode == .valid_zipredictcode[["disp"]]) {
    ## handled separately by caller
  } else {
    stop("Invalid ziPredictCode: ", ziPredictCode)
  }

  list(mu = mu, eta = eta)
}

linkfun_rtmb <- function(mu, link) {
  switch(
    names(link),
    log = log(mu),
    identity = mu,
    sqrt = sqrt(mu),
    logit = log(mu / (1 - mu)),
    probit = stats::qnorm(mu),
    cloglog = log(-log(1 - mu)),
    inverse = 1 / mu,
    lambertW = stop("linkfun for lambertW not yet implemented"),
    stop("link not yet implemented for prediction aggregation: ", names(link))
  )
}

family_name_rtmb <- function(family) {
  if (is.list(family) && !is.null(family$family)) {
    return(family$family)
  }

  family_name <- names(family)
  if (length(family_name) == 0L) {
    family_name <- names(.valid_family)[match(family, .valid_family)]
  }
  family_name
}

dcauchy_rtmb <- function(x, location, scale, log = FALSE) {
  resid <- (x - location) / scale
  ans <- -log(pi) - log(scale) - log1p(resid * resid)
  if (log) ans else exp(ans)
}

dlkj_rtmb <- function(x, eta, log = FALSE) {
  "[<-" <- RTMB::ADoverload("[<-")

  len <- length(x)
  if (len == 0) {
    return(if (log) 0 else 1)
  }

  n <- (1 + sqrt(1 + 8 * len)) / 2
  if (abs(n - round(n)) > sqrt(.Machine$double.eps)) {
    stop("Invalid number of LKJ correlation parameters: ", len)
  }
  n <- as.integer(round(n))

  L <- diag(n)
  k <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        L[i, j] <- x[k]
        k <- k + 1L
      }
    }
  }

  row_sums <- L * L
  log_det_x <- 0
  for (i in seq_len(n)) {
    log_det_x <- log_det_x - log(sum(row_sums[i, ]))
  }
  ans <- (eta - 1) * log_det_x
  if (log) ans else exp(ans)
}

prior_nll <- function(beta, betazi, betadisp, theta, thetazi, psi,
                      prior_distrib, prior_whichpar, prior_elstart,
                      prior_elend, prior_npar, prior_params) {
  nll <- 0
  par_ind <- 1L

  for (i in seq_along(prior_distrib)) {
    parvec <- switch(
      as.character(prior_whichpar[i]),
      "0" = beta,
      "1" = betazi,
      "2" = betadisp,
      "10" = theta,
      "20" = thetazi,
      "30" = psi,
      stop("Unknown prior parameter vector code: ", prior_whichpar[i])
    )

    par_start <- prior_elstart[i] + 1L
    par_end <- prior_elend[i] + 1L
    if (par_start > par_end) {
      par_ind <- par_ind + prior_npar[i]
      next
    }
    if (par_start < 1L || par_end > length(parvec)) {
      stop(
        "Bad prior index for prior ", i, ": requested elements ",
        prior_elstart[i], ":", prior_elend[i],
        " in a parameter vector of length ", length(parvec)
      )
    }

    if (prior_distrib[i] == .valid_prior[["lkj"]]) {
      corpars <- parvec[par_start:par_end]
      nll <- nll - dlkj_rtmb(
        corpars,
        prior_params[par_ind],
        log = TRUE
      )
    } else {
      for (j in par_start:par_end) {
        parval <- parvec[j]
        logpriorval <- switch(
          as.character(prior_distrib[i]),
          "0" = RTMB::dnorm(
            parval,
            mean = prior_params[par_ind],
            sd = prior_params[par_ind + 1L],
            log = TRUE
          ),
          "1" = {
            location <- prior_params[par_ind]
            scale <- prior_params[par_ind + 1L]
            df <- prior_params[par_ind + 2L]
            RTMB::dt((parval - location) / scale, df = df, log = TRUE) -
              log(scale)
          },
          "2" = dcauchy_rtmb(
            parval,
            location = prior_params[par_ind],
            scale = prior_params[par_ind + 1L],
            log = TRUE
          ),
          "10" = {
            shape <- prior_params[par_ind + 1L]
            scale <- prior_params[par_ind] / prior_params[par_ind + 1L]
            RTMB::dgamma(exp(parval), shape = shape, scale = scale,
                         log = TRUE)
          },
          stop("Prior distribution not implemented: ", prior_distrib[i])
        )
        nll <- nll - logpriorval
      }
    }
    par_ind <- par_ind + prior_npar[i]
  }
  nll
}

## Variables injected into rtmb_tpl() by RTMB::getAll()
utils::globalVariables(c(
  "X", "XS", "Z", "offset", "terms", "family", "link", "weights", "size",
  "beta", "b", "theta",
  "Xzi", "XziS", "Zzi", "zioffset", "termszi",
  "betazi", "bzi", "thetazi",
  "Xdisp", "XdispS", "Zdisp", "dispoffset", "termsdisp",
  "betadisp", "bdisp", "thetadisp",
  "psi", "combinom_disp_link", "ziPredictCode", "doPredict",
  "whichPredict", "aggregate",
  "prior_distrib", "prior_whichpar", "prior_elstart", "prior_elend",
  "prior_npar", "prior_params"
))

rtmb_tpl <- function(parameters, data) {
  RTMB::getAll(data, parameters)
  family_name <- family_name_rtmb(family)
  ## Keep the original response for NA and structural-zero checks; OBS() may
  ## replace yobs with a simulation or OSA reference. During OSA calculations
  ## yobs is moved from data into parameters, so data$yobs may be absent.
  yobs_obs <- if (!is.null(data$yobs)) data$yobs else osa_value(yobs)
  yobs <- RTMB::OBS(yobs)

  nll <- 0

  ## Random-effects contribution; translated from glmmTMB.cpp:900-903
  cond_re <- allterms_nll(b, theta, terms)
  zi_re <- allterms_nll(bzi, thetazi, termszi)
  disp_re <- allterms_nll(bdisp, thetadisp, termsdisp)
  nll <- nll + cond_re$nll + zi_re$nll + disp_re$nll
  b <- cond_re$u
  bzi <- zi_re$u
  bdisp <- disp_re$u

  ## Conditional linear predictor and inverse link; adapted from
  ## glmmTMB.cpp:833, 911-918, and 934-937
  sparseX <- nrow(X) == 0 && ncol(X) == 0
  Xc <- if (sparseX) XS else X
  eta <- Xc %*% beta + Z %*% b + offset
  eta <- as.vector(eta)

  mu <- switch(
    names(link),
    log = exp(eta),
    identity = eta,
    sqrt = eta * eta,
    logit = 1 / (1 + exp(-eta)),
    probit = RTMB::pnorm(eta),
    cloglog = 1 - exp(-exp(eta)),
    inverse = 1 / eta,
    lambertW = exp(eta) * exp(exp(eta)),
    stop(
      "link not yet implemented: ", names(link),
      "; implemented links are: log, identity, sqrt, logit, probit, ",
      "cloglog, inverse, lambertW"
    )
  )

  ## Zero-inflation linear predictor; adapted from
  ## glmmTMB.cpp:836, 880, and 919-925
  has_zi <- length(betazi) > 0 || length(bzi) > 0
  if (has_zi) {
    sparseXzi <- nrow(Xzi) == 0 && ncol(Xzi) == 0
    Xzic <- if (sparseXzi) XziS else Xzi
    etazi <- Xzic %*% betazi + Zzi %*% bzi + zioffset
    etazi <- as.vector(etazi)
  }

  ## Dispersion linear predictor; adapted from
  ## glmmTMB.cpp:839, 926-932, and 939
  sparseXdisp <- nrow(Xdisp) == 0 && ncol(Xdisp) == 0
  Xdispc <- if (sparseXdisp) XdispS else Xdisp
  etadisp <- Xdispc %*% betadisp + Zdisp %*% bdisp + dispoffset
  etadisp <- as.vector(etadisp)
  phi <- exp(etadisp)
  if (family_name == "combinomial" && combinom_disp_link == 1L) {
    phi <- etadisp
  }

  ## Observation likelihoods; adapted from glmmTMB.cpp:961-978,
  ## 1095-1101, and 1180-1199
  i <- !is.na(yobs_obs) | inherits(yobs, "simref")
  yobs_i <- yobs[i]
  keep <- osa_keep(yobs_i)
  eta_zi <- if (has_zi) etazi[i] else NULL
  logit_mu <- function() logit_mu_rtmb(eta, link)
  log_mu <- function() log_mu_rtmb(eta, link)
  log_var_minus_mu <- function() {
    log_var_minus_mu_rtmb(family_name, log_mu(), etadisp, psi)
  }

  tmp_loglik <- switch(
    family_name,
    poisson = dZI(RTMB::dpois)(yobs_i, lambda = mu[i], eta_zi = eta_zi, log = TRUE,
                               is_zero = yobs_obs[i] == 0),
    truncated_poisson = dZI(dtruncated_poisson_rtmb)(
      yobs_i, lambda = mu[i], eta_zi = eta_zi, log = TRUE,
      is_zero = yobs_obs[i] == 0
    ),
    gaussian = dZI(RTMB::dnorm)(
      yobs_i, mean = mu[i], sd = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the Gamma_family case in glmmTMB.cpp:991-996.
    Gamma = dZI(dgamma_rtmb)(
      yobs_i, mean = mu[i], shape = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the beta_family case in glmmTMB.cpp:997-1002.
    beta = dZI(dbeta_rtmb)(
      yobs_i, mean = mu[i], phi = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the ordbeta_family case in glmmTMB.cpp:1004-1031.
    ordbeta = dZI(dordbeta_rtmb)(
      yobs_i, eta = eta[i], mean = mu[i], phi = phi[i], cutpoints = psi,
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the lognormal_family case in glmmTMB.cpp:1164-1179.
    lognormal = dZI(dlognormal_rtmb)(
      yobs_i, mean = mu[i], sd = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the t_family case in glmmTMB.cpp:1182-1190.
    t = dZI(dt_rtmb)(
      yobs_i, mean = mu[i], scale = phi[i], df = exp(psi[1L]),
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the skewnormal_family case in glmmTMB.cpp:975-980.
    skewnormal = dZI(dskewnormal_rtmb)(
      yobs_i, mean = mu[i], sd = phi[i], alpha = psi[1L],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the tweedie_family case in glmmTMB.cpp:1155-1163.
    tweedie = dZI(RTMB::dtweedie)(
      yobs_i, mu = mu[i], phi = phi[i],
      p = 1 / (1 + exp(-psi[1L])) + 1,
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the bell_family case in glmmTMB.cpp:1191-1200.
    bell = dZI(dbell_rtmb)(
      yobs_i, mean = mu[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the binomial_family case in glmmTMB.cpp:979-983.
    binomial = dZI(dbinom_robust_rtmb)(
      yobs_i, size = size[i], logit_p = logit_mu()[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from betabinomial_family, glmmTMB.cpp:1037-1047.
    betabinomial = {
      logit_p <- logit_mu()[i]
      dZI(dbetabinom_robust_rtmb)(
        yobs_i,
        log_shape1 = -RTMB::logspace_add(0, -logit_p) + etadisp[i],
        log_shape2 = -RTMB::logspace_add(0, logit_p) + etadisp[i],
        size = size[i], eta_zi = eta_zi, log = TRUE,
        is_zero = yobs_obs[i] == 0
      )
    },
    ## Translated from combinomial_family, glmmTMB.cpp:1049-1067.
    combinomial = dZI(dcombinom2_rtmb)(
      yobs_i, size = size[i], mean = mu[i] * size[i], nu = phi[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the nbinom1_family case in glmmTMB.cpp:1042-1056.
    nbinom1 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_nbinom1_family, glmmTMB.cpp:1042-1064.
    truncated_nbinom1 = dZI(dtruncated_nbinom1_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      log_phi = etadisp[i], eta_zi = eta_zi, log = TRUE,
      is_zero = yobs_obs[i] == 0),
    ## Translated from the nbinom2_family case in glmmTMB.cpp:1066-1075.
    nbinom2 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the nbinom12_family case in glmmTMB.cpp:1084-1094.
    nbinom12 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the genpois_family case in glmmTMB.cpp:1128-1133.
    genpois = dZI(dgenpois_rtmb)(
      yobs_i, theta = mu[i] / sqrt(phi[i]),
      lambda = 1 - 1 / sqrt(phi[i]), eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_genpois_family, glmmTMB.cpp:1134-1139.
    truncated_genpois = dZI(dtruncated_genpois_rtmb)(
      yobs_i, theta = mu[i] / sqrt(phi[i]),
      lambda = 1 - 1 / sqrt(phi[i]), eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the compois_family case in glmmTMB.cpp:1115-1119.
    compois = dZI(dcompois2_rtmb)(
      yobs_i, mean = mu[i], nu = 1 / phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_compois_family, glmmTMB.cpp:1121-1127.
    truncated_compois = dZI(dtruncated_compois2_rtmb)(
      yobs_i, mean = mu[i], nu = 1 / phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_nbinom2_family, glmmTMB.cpp:1066-1081.
    truncated_nbinom2 = dZI(dtruncated_nbinom2_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      log_size = etadisp[i], eta_zi = eta_zi, log = TRUE,
      is_zero = yobs_obs[i] == 0),
    stop(
      "distribution not implemented yet for use with RTMB backend: ",
      family_name
    )
  )

  nll <- nll - sum(keep * weights[i] * tmp_loglik)

  ## Prior contribution; translated from glmmTMB.cpp:1203-1267
  nll <- nll + prior_nll(
    beta = beta,
    betazi = betazi,
    betadisp = betadisp,
    theta = theta,
    thetazi = thetazi,
    psi = psi,
    prior_distrib = prior_distrib,
    prior_whichpar = prior_whichpar,
    prior_elstart = prior_elstart,
    prior_elend = prior_elend,
    prior_npar = prior_npar,
    prior_params = prior_params
  )

  ## Prediction output; translated from glmmTMB.cpp:1353-1379
  mu_pred_all <- mu
  eta_pred_all <- eta

  ## Convert untruncated mean to the conditional mean of truncated distribution
  ## translated from glmmTMB.cpp:1331-1334
  if (family_name == "truncated_poisson") {
    mu_vector <- mu[seq_along(mu)]
    log_nzprob_pred <- log_nzprob_truncated_poisson_rtmb(mu_vector)
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_nbinom1") {
    mu_vector <- mu[seq_along(mu)]
    etadisp_vector <- etadisp[seq_along(etadisp)]
    log_nzprob_pred <- log_nzprob_truncated_nbinom1_rtmb(
      mu_vector,
      etadisp_vector
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_nbinom2") {
    log_mu_value <- log_mu()
    log_mu_vector <- log_mu_value[seq_along(log_mu_value)]
    etadisp_vector <- etadisp[seq_along(etadisp)]
    log_nzprob_pred <- log_nzprob_truncated_nbinom2_rtmb(
      log_mu_vector,
      etadisp_vector
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_genpois") {
    mu_vector <- mu[seq_along(mu)]
    phi_vector <- phi[seq_along(phi)]
    log_nzprob_pred <- log_nzprob_truncated_genpois_rtmb(
      mu_vector / sqrt(phi_vector)
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_compois") {
    mu_vector <- mu[seq_along(mu)]
    phi_vector <- phi[seq_along(phi)]
    log_nzprob_pred <- log_nzprob_truncated_compois_rtmb(
      mu_vector,
      1 / phi_vector
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  }

  if (has_zi || ziPredictCode == .valid_zipredictcode[["prob"]]) {
    zi_pred <- apply_zi_prediction(
      mu = mu_pred_all,
      eta = eta_pred_all,
      etazi = etazi,
      ziPredictCode = ziPredictCode
    )
    mu_pred_all <- zi_pred$mu
    eta_pred_all <- zi_pred$eta
  }

  if (ziPredictCode == .valid_zipredictcode[["disp"]]) {
    mu_pred_all <- if (family_name == "Gamma") 1 / sqrt(phi) else phi
    eta_pred_all <- etadisp
  }

  mu_predict <- mu_pred_all[whichPredict]
  eta_predict <- eta_pred_all[whichPredict]

  if (length(aggregate) > 0) {
    if (length(aggregate) != length(mu_predict)) {
      stop(
        "'aggregate' wrong size; got length ", length(aggregate),
        " but prediction length is ", length(mu_predict)
      )
    }

    "[<-" <- RTMB::ADoverload("[<-")
    n_aggregate <- max(as.integer(aggregate))
    tmp <- rep(mu_predict[1L] * 0, n_aggregate)
    for (j in seq_along(mu_predict)) {
      tmp[as.integer(aggregate[j])] <- tmp[as.integer(aggregate[j])] +
        mu_predict[j]
    }

    mu_predict <- tmp
    eta_predict <- linkfun_rtmb(mu_predict, link)
  }

  corr <- cond_re$corr
  sd <- cond_re$sd
  corrzi <- zi_re$corr
  sdzi <- zi_re$sd
  corrdisp <- disp_re$corr
  sddisp <- disp_re$sd
  fact_load <- cond_re$fact_load

  REPORT(corr)
  REPORT(sd)
  REPORT(corrzi)
  REPORT(sdzi)
  REPORT(corrdisp)
  REPORT(sddisp)
  REPORT(fact_load)
  REPORT(b)
  REPORT(bzi)
  REPORT(bdisp)
  REPORT(mu_predict)
  REPORT(eta_predict)

  if (doPredict == 1) {
    ADREPORT(mu_predict)
  } else if (doPredict == 2) {
    ADREPORT(eta_predict)
  } else if (doPredict == 3) {
    ADREPORT(b)
    ADREPORT(bzi)
    ADREPORT(bdisp)
  }

  nll
}
## Partition the concatenated random effects and covariance parameters by term
## Term slicing is translated from allterms_nll() in glmmTMB.cpp:803-826
allterms_nll <- function(u, theta, terms) {
  "[<-" <- RTMB::ADoverload("[<-")

  nll <- 0
  corr <- vector("list", length(terms))
  sd <- vector("list", length(terms))
  fact_load <- vector("list", length(terms))
  names(corr) <- names(terms)
  names(sd) <- names(terms)
  names(fact_load) <- names(terms)

  if (length(terms) == 0) {
    output_u <- if (inherits(u, "simref")) u$value else u
    return(list(
      nll = nll, corr = corr, sd = sd, fact_load = fact_load, u = output_u
    ))
  }

  transformed_u <- u
  upointer <- 0L
  tpointer <- 0L
  np <- 0L

  for (i in seq_along(terms)) {
    term <- terms[[i]]
    nr <- term$blockSize * term$blockReps
    ## A zero-length theta block reuses the prev term's covariance parameters
    emptyTheta <- term$blockNumTheta == 0

    if (!emptyTheta) {
      np <- term$blockNumTheta
      theta_start <- tpointer + 1L
    } else {
      theta_start <- tpointer - np + 1L
    }

    useg <- u[(upointer + 1L):(upointer + nr)]

    if (np > 0) {
        tseg <- theta[theta_start:(theta_start + np - 1L)]
    } else {
      tseg <- numeric(0)
    }

    ans <- termwise_nll(useg, tseg, term)
    nll <- nll + ans$nll
    corr[[i]] <- ans$corr
    sd[[i]] <- ans$sd
    fact_load[[i]] <- ans$fact_load
    if (!inherits(transformed_u, "simref")) {
      transformed_u[(upointer + 1L):(upointer + nr)] <- ans$u
    }

    upointer <- upointer + nr
    tpointer <- tpointer + term$blockNumTheta
  }

  output_u <- if (inherits(transformed_u, "simref")) {
    transformed_u$value
  } else {
    transformed_u
  }

  list(
    nll = nll,
    corr = corr,
    sd = sd,
    fact_load = fact_load,
    u = output_u
  )
}

## Construct the correlation matrix used by TMB's
## density::UNSTRUCTURED_CORR_t. TMB fills the lower triangle row-wise,
## whereas matrix lower-triangle assignment in R fills it column-wise
## leading to a different theta ordering for dim >= 4
tmb_unstructured_corr <- function(n, theta) {
  "[<-" <- RTMB::ADoverload("[<-")

  expected <- n * (n - 1L) / 2L
  if (length(theta) != expected) {
    stop(
      "Expected ", expected, " correlation parameters for unstructured ",
      n, " by ", n, " correlation matrix, got ", length(theta)
    )
  }

  L <- diag(n)
  k <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        L[i, j] <- theta[k]
        k <- k + 1L
      }
    }
  }

  llt <- L %*% t(L)
  corr <- llt
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      corr[i, j] <- llt[i, j] / sqrt(llt[i, i] * llt[j, j])
    }
  }
  corr
}

## Simulation factor for density::UNSTRUCTURED_CORR_t, used in
## glmmTMB.cpp:407-425. This preserves the C++ simulation order for us()
## random effects instead of delegating simulation to RTMB::dmvnorm().
tmb_unstructured_sim_factor <- function(n, theta) {
  expected <- n * (n - 1L) / 2L
  if (length(theta) != expected) {
    stop(
      "Expected ", expected, " correlation parameters for unstructured ",
      n, " by ", n, " simulation factor, got ", length(theta)
    )
  }

  L <- diag(n)
  k <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        L[i, j] <- theta[k]
        k <- k + 1L
      }
    }
  }

  L / sqrt(rowSums(L * L))
}

simulate_tmb_unstructured <- function(U, sd, corr_par) {
  n <- nrow(U)
  reps <- ncol(U)
  sim_factor <- tmb_unstructured_sim_factor(n, as.vector(corr_par))
  sim_sd <- as.vector(sd)

  for (j in seq_len(reps)) {
    U_col <- U[, j]
    U_col[] <- sim_sd * as.vector(sim_factor %*% stats::rnorm(n))
  }
}

## Evaluate one random-effects term under its covariance structure
## Translation of the currently supported cases in
## termwise_nll(), glmmTMB.cpp:358-799
termwise_nll <- function(U, theta, term) {
  ## Preserve automatic differentiation when filling correlation matrices
  "[<-" <- RTMB::ADoverload("[<-")

  block_code <- term$blockCode
  name <- if (is.character(block_code) && length(block_code) == 1L) {
    block_code
  } else {
    block_name <- names(block_code)
    if (length(block_name) == 0L) {
      names(.valid_covstruct)[match(block_code, .valid_covstruct)]
    } else {
      block_name[1L]
    }
  }
  supported <- c(
    "diag", "homdiag", "us", "cs", "homcs", "toep", "homtoep",
    "ar1", "hetar1", "ou", "exp", "gau", "mat", "rr", "propto", "equalto"
  )

  if (!name %in% supported) {
    stop(
      "covariance structure not yet implemented: ", name,
      "; implemented structures are: ", paste(supported, collapse = ", ")
    )
  }

  n <- term$blockSize
  reps <- term$blockReps
  dim(U) <- c(n, reps)

  rr_rank <- NA_integer_
  if (name == "rr") {
    ntheta <- length(theta)
    rank_discriminant <- (2 * n + 1)^2 - 8 * ntheta
    if (rank_discriminant < 0) {
      stop(
        "Invalid covariance parameter count for 'rr': ", ntheta,
        "; rank discriminant is ", rank_discriminant,
        ", so no real-valued rank can be inferred for block size ", n
      )
    }
    rank_value <- (
      2 * n + 1 - sqrt(rank_discriminant)
    ) / 2
    rr_rank <- as.integer(round(rank_value))
    valid_rank <- is.finite(rank_value) &&
      abs(rank_value - rr_rank) < sqrt(.Machine$double.eps) &&
      rr_rank >= 1L &&
      rr_rank <= n

    if (!valid_rank) {
      stop(
        "Invalid covariance parameter count for 'rr': ", ntheta,
        "; inferred rank value is ", rank_value,
        ", rounded rank is ", rr_rank,
        ", valid ranks are integers from 1 to ", n
      )
    }
  }

  expected_num_theta <- switch(
    name,
    diag = n,
    homdiag = 1L,
    us = n * (n + 1L) / 2L,
    cs = n + 1L,
    homcs = 2L,
    toep = 2L * n - 1L,
    homtoep = n,
    ar1 = 2L,
    hetar1 = n + 1L,
    ou = 2L,
    exp = 2L,
    gau = 2L,
    mat = 3L,
    rr = n * rr_rank - (rr_rank - 1L) * rr_rank / 2L,
    propto = n * (n + 1L) / 2L + 1L,
    equalto = n * (n + 1L) / 2L
  )
  if (length(theta) != expected_num_theta) {
    stop(
      "Expected ", expected_num_theta, " covariance parameters for '",
      name, "', got ", length(theta)
    )
  }

  if (name == "rr") {
    ## Reduced-rank covariance; glmmTMB.cpp:698-761. The optimized random
    ## effects are spherical, while the linear predictor uses Lambda %*% u.
    nll <- 0
    simulation <- inherits(U, "simref")

    if (simulation && !term$simCode %in% .valid_simcode) {
      stop(
        "unknown simCode for rr covariance structure: ", term$simCode,
        "; known simCodes are: ",
        paste(names(.valid_simcode), .valid_simcode, sep = "=",
              collapse = ", ")
      )
    }

    if (!simulation || term$simCode == .valid_simcode[["random"]]) {
      for (j in seq_len(reps)) {
        nll <- nll - sum(RTMB::dnorm(U[, j], 0, 1, log = TRUE))
      }
    } else if (term$simCode == .valid_simcode[["zero"]]) {
      U[] <- 0
    } else {
      U[] <- U$getOrig(seq_along(U))
    }

    Lambda <- matrix(0, n, rr_rank)
    lam_diag <- head(theta, rr_rank)
    lam_lower <- utils::tail(theta, length(theta) - rr_rank)

    Lambda[row(Lambda) == col(Lambda)] <- lam_diag
    Lambda[row(Lambda) > col(Lambda)] <- lam_lower

    if (term$simCode != .valid_simcode[["fix"]]) {
      for (j in seq_len(reps)) {
        transformed_column <- Lambda %*% U[seq_len(rr_rank), j]
        if (simulation) {
          U_column <- U[, j]
          U_column[] <- transformed_column
        } else {
          U[, j] <- transformed_column
        }
      }
    }

    report_corr <- matrix(numeric(0), 0, 0)
    report_sd <- numeric(0)
    if (term$fullCor == 1L) {
      covariance <- Lambda %*% t(Lambda)
      report_sd <- sqrt(diag(covariance))
      report_corr <- covariance /
        (report_sd %*% t(report_sd))
    }

    return(list(
      nll = nll,
      corr = report_corr,
      sd = report_sd,
      fact_load = Lambda,
      u = if (simulation) NULL else as.vector(U)
    ))
  }

  ## Homogeneous structures use one standard-deviation parameter;
  ## heterogeneous structures use one parameter per term component.
  homogeneous <- c(
    "homdiag", "homcs", "homtoep", "ar1", "ou", "exp", "gau", "mat"
  )
  hetvar <- !name %in% homogeneous
  n_sd_par <- if (hetvar) n else 1L

  logsd <- if (hetvar) {
    head(theta, n)
  } else {
    rep(theta[1L], n)
  }

  sd <- exp(logsd)
  corr_par <- theta[-seq_len(n_sd_par)]

  ## propto uses an unstructured correlation matrix with an additional
  ## parameter that proportionally scales the covariance matrix.
  if (name == "propto") {
    loglambda <- utils::tail(corr_par, 1L)
    corr_par <- head(corr_par, -1L)
    sd <- exp(logsd + loglambda / 2)
  }

  ## Remove the "hom" prefix because homogeneous and heterogeneous
  ## variants differ only in their standard-deviation parameterization.
  cov_structure <- sub("^hom", "", name)

  ## propto and equalto use the unstructured correlation parameterization.
  density_structure <- if (cov_structure %in% c("propto", "equalto")) {
    "us"
  } else if (cov_structure == "hetar1") {
    "ar1"
  } else {
    cov_structure
  }

  C <- switch(
    density_structure,

    ## Diagonal covariance; glmmTMB.cpp:358-405
    diag = {
      matrix(numeric(0), 0, 0)
    },

    ## Unstructured covariance; glmmTMB.cpp:407-440
    us = {
      tmb_unstructured_corr(n, corr_par)
    },

    ## Compound-symmetry covariance; glmmTMB.cpp:441-473
    cs = {
      a <- 1 / (n - 1)
      rho <- (1 / (1 + exp(-corr_par[1L]))) * (1 + a) - a
      corr <- diag(n)
      corr[row(corr) != col(corr)] <- rho
      corr
    },

    ## Toeplitz covariance; glmmTMB.cpp:474-506
    toep = {
      corr_params <- corr_par / sqrt(1 + corr_par^2)
      corr <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          corr[i, j] <- if (i == j) {
            1
          } else {
            corr_params[abs(i - j)]
          }
        }
      }

      corr
    },

    ## Homogeneous AR(1) covariance; glmmTMB.cpp:507-590
    ar1 = {
      phi <- corr_par[1L] / sqrt(1 + corr_par[1L]^2)
      corr <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          corr[i, j] <- phi^abs(i - j)
        }
      }
      corr
    },

    ## OU covariance; glmmTMB.cpp:593-650
    ou = {
      times <- term$times
      if (length(times) != n) {
        stop(
          "OU time vector length must equal block size; got length(times)=",
          length(times), " and blockSize=", n
        )
      }
      decay <- exp(corr_par[1L])
      corr <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          corr[i, j] <- exp(
            -decay * abs(times[i] - times[j])
          )
        }
      }
      corr
    },

    ## Exponential spatial covariance; glmmTMB.cpp:653-700
    exp = {
      spatial_dist <- term$dist
      spatial_dim <- dim(spatial_dist)
      if (length(spatial_dim) != 2L || any(spatial_dim != n)) {
        stop(
          "Dimension of distance matrix must equal block size for ", name,
          "; got dim(dist)=",
          paste(spatial_dim, collapse = " x "),
          " and blockSize=", n
        )
      }
      corr <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          corr[i, j] <- if (i == j) {
            1
          } else {
            exp(-spatial_dist[i, j] * exp(-corr_par[1L]))
          }
        }
      }
      corr
    },

    ## Gaussian spatial covariance; glmmTMB.cpp:653-700
    gau = {
      spatial_dist <- term$dist
      spatial_dim <- dim(spatial_dist)
      if (length(spatial_dim) != 2L || any(spatial_dim != n)) {
        stop(
          "Dimension of distance matrix must equal block size for ", name,
          "; got dim(dist)=",
          paste(spatial_dim, collapse = " x "),
          " and blockSize=", n
        )
      }
      corr <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          corr[i, j] <- if (i == j) {
            1
          } else {
            exp(
              -(spatial_dist[i, j]^2) * exp(-2 * corr_par[1L])
            )
          }
        }
      }
      corr
    },

    ## Matern covariance; glmmTMB.cpp:653-700
    mat = {
      spatial_dist <- term$dist
      spatial_dim <- dim(spatial_dist)
      if (length(spatial_dim) != 2L || any(spatial_dim != n)) {
        stop(
          "Dimension of distance matrix must equal block size for ", name,
          "; got dim(dist)=",
          paste(spatial_dim, collapse = " x "),
          " and blockSize=", n
        )
      }
      range <- exp(corr_par[1L])
      smoothness <- exp(corr_par[2L])
      corr <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i == j) {
            corr[i, j] <- 1
          } else {
            scaled_dist <- spatial_dist[i, j] / range
            corr[i, j] <-
              scaled_dist^smoothness * RTMB::besselK(
                scaled_dist,
                smoothness
              ) /
              (exp(lgamma(smoothness)) * 2^(smoothness - 1))
          }
        }
      }
      corr
    },
    stop(
      "covariance structure not yet implemented: ", name,
      "; implemented density structures are: diag, us, cs, toep, ar1, ",
      "ou, exp, gau, mat"
    )
  )

  simulation <- inherits(U, "simref")
  simulate_density <- TRUE
  if (simulation) {
    if (!term$simCode %in% .valid_simcode) {
      stop(
        "unknown simCode for ", name, " covariance structure: ",
        term$simCode,
        "; known simCodes are: ",
        paste(names(.valid_simcode), .valid_simcode, sep = "=",
              collapse = ", ")
      )
    }

    flexible_simulation <- c("diag", "us", "ar1", "hetar1", "ou")
    random_only_simulation <- c(
      "homdiag", "cs", "homcs", "toep", "homtoep",
      "exp", "gau", "mat"
    )

    if (name %in% flexible_simulation) {
      if (term$simCode == .valid_simcode[["zero"]]) {
        U[] <- 0
        simulate_density <- FALSE
      } else if (term$simCode == .valid_simcode[["fix"]]) {
        U[] <- U$getOrig(seq_along(U))
        simulate_density <- FALSE
      } else if (
        name == "us" && term$simCode == .valid_simcode[["random"]]
      ) {
        simulate_tmb_unstructured(U, sd, corr_par)
        simulate_density <- FALSE
      }
    } else if (
      name %in% random_only_simulation &&
      term$simCode != .valid_simcode[["random"]]
    ) {
      stop(
        "simCode '",
        names(.valid_simcode)[match(term$simCode, .valid_simcode)],
        "' is not implemented for ", name,
        " covariance structure; only random simulation is currently supported"
      )
    }
  }

  ## Diagonal structures factor into univariate normal densities;
  ## correlated structures use a scaled multivariate normal density.
  if (!simulate_density) {
    nll <- 0
  } else if (density_structure == "diag") {
    nll <- 0

    for (k in seq_len(n)) {
      nll <- nll - sum(RTMB::dnorm(U[k, ], 0, sd[k], log = TRUE))
    }
  } else {
    ## Keep scale dimensions identical to t(U). A bare vector is ambiguous to
    ## RTMB::dmvnorm when there is exactly one block repetition.
    scale_matrix <- rep(sd, reps)
    dim(scale_matrix) <- c(n, reps)
    scale_matrix <- t(scale_matrix)
    nll <- -sum(RTMB::dmvnorm(t(U), Sigma = C, log = TRUE, scale = scale_matrix)
    )
  }

  ## Match C++ full-correlation reporting; equalto always reports its matrix.
  report_corr <- C
  if (name %in% c("ar1", "hetar1") && term$fullCor == 0) {
    report_corr <- matrix(phi, 1L, 1L)
  }
  if (name == "ou" && term$fullCor == 0) {
    report_corr <- matrix(decay, 1L, 1L)
  }
  if (name %in% c("exp", "gau", "mat") && term$fullCor == 0) {
    report_corr <- matrix(numeric(0), 0, 0)
  }

  conditional_full_cor <- c(
    "us", "cs", "homcs", "toep", "homtoep", "propto"
  )
  if (name %in% conditional_full_cor && term$fullCor == 0) {
    report_corr <- matrix(NaN, 1L, 1L)
  }

  report_sd <- if (name == "ar1" || (name == "ou" && term$fullCor == 0)) {
    sd[1L]
  } else {
    sd
  }
  list(
    nll = nll,
    corr = report_corr,
    sd = report_sd,
    fact_load = matrix(numeric(0), 0, 0),
    u = if (inherits(U, "simref")) NULL else as.vector(U)
  )
}
