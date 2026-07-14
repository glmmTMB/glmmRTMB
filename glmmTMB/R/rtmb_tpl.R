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
    {
      mu <- switch(
        names(link),
        probit = {
          probit_eta <- if (inherits(eta, "simref")) eta$value else eta
          if (inherits(probit_eta, "Matrix")) {
            probit_eta <- as.matrix(probit_eta)
          }
          RTMB::pnorm(probit_eta)
        },
        cloglog = 1 - exp(-exp(eta)),
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
        probit = {
          probit_eta <- if (inherits(eta, "simref")) eta$value else eta
          if (inherits(probit_eta, "Matrix")) {
            probit_eta <- as.matrix(probit_eta)
          }
          RTMB::pnorm(probit_eta)
        },
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
#'   `dnorm_tmb()` that also supports RTMB simulation through `simref` objects.
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
  if (inherits(prob_nonzero, "simref")) {
    prob_nonzero <- prob_nonzero$value
  }

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

## matches tmb's dnorm() arithmetic ordering
dnorm_tmb <- local({
  log_density <- RTMB::Vectorize(
    function(x, mean, sd) {
      z <- (x - mean) / sd
      -log(sqrt(2 * pi)) - log(sd) - 0.5 * z * z
    },
    vectorize.args = c("x", "mean", "sd")
  )

  function(x, mean = 0, sd = 1, log = FALSE) {
    if (inherits(x, "simref")) {
      if (inherits(mean, "simref")) {
        mean <- mean$value
      }
      if (inherits(sd, "simref")) {
        sd <- sd$value
      }
      if (inherits(mean, "Matrix")) {
        mean <- as.matrix(mean)
      }
      if (inherits(sd, "Matrix")) {
        sd <- as.matrix(sd)
      }
      x[] <- stats::rnorm(length(x), mean = as.vector(mean), sd = as.vector(sd))
      return(rep(0, length(x)))
    }
    ans <- log_density(x, mean, sd)
    if (log) ans else exp(ans)
  }
})

## Fitting uses RTMB::dbinom_robust() to match glmmTMB.cpp:981.
## Simulation follows glmmTMB.cpp:982 but delegates to RTMB::dbinom(),
## because RTMB::dbinom_robust() does not provide rbinom_robust().
dbinom_robust_rtmb <- function(x, size, logit_p, log = FALSE) {
  if (inherits(x, "simref")) {
    prob <- 1 / (1 + exp(-logit_p))
    if (inherits(size, "simref")) {
      size <- size$value
    }
    if (inherits(prob, "simref")) {
      prob <- prob$value
    }
    return(RTMB::dbinom(x, size = size, prob = prob, log = log))
  }
  RTMB::dbinom_robust(x, size = size, logit_p = logit_p, log = log)
}

## Fitting translates glmmTMB.cpp:1042-1075 for nbinom1/nbinom2:
## both families use dnbinom_robust(log_mu, log_var_minus_mu). Simulation
## follows the same mean/variance by converting back to size/mu.
dnbinom_robust_rtmb <- function(x, log_mu, log_var_minus_mu, log = FALSE) {
  if (inherits(x, "simref")) {
    if (inherits(log_mu, "simref")) {
      log_mu <- log_mu$value
    }
    if (inherits(log_var_minus_mu, "simref")) {
      log_var_minus_mu <- log_var_minus_mu$value
    }
    if (inherits(log_mu, "Matrix")) {
      log_mu <- as.matrix(log_mu)
    }
    if (inherits(log_var_minus_mu, "Matrix")) {
      log_var_minus_mu <- as.matrix(log_var_minus_mu)
    }
    mu <- exp(as.vector(log_mu))
    size <- exp(as.vector(2 * log_mu - log_var_minus_mu))
    x[] <- stats::rnbinom(length(x), size = size, mu = mu)
    return(rep(0, length(x)))
  }
  RTMB::dnbinom_robust(x, log_mu = log_mu, log_var_minus_mu = log_var_minus_mu,
                       log = log)
}

## zero-truncated poisson density
dtruncated_poisson_rtmb <- function(x, lambda, log = FALSE) {
  if (inherits(x, "simref")) {
    if (inherits(lambda, "simref")) {
      lambda <- lambda$value
    }
    ## Draw from the conditional upper tail
    # expm1() helps for accuracy when lambda is close to zero.
    upper_prob <- stats::runif(length(x)) * (-expm1(-lambda))
    x[] <- stats::qpois(upper_prob, lambda = lambda, lower.tail = FALSE)
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
          "0" = dnorm_tmb(
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
  "psi", "ziPredictCode", "doPredict", "whichPredict", "aggregate",
  "prior_distrib", "prior_whichpar", "prior_elstart", "prior_elend",
  "prior_npar", "prior_params"
))

rtmb_tpl <- function(parameters, data) {
  RTMB::getAll(data, parameters)
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

  mu <- switch(
    names(link),
    log = exp(eta),
    identity = eta,
    sqrt = eta * eta,
    logit = 1 / (1 + exp(-eta)),
    probit = {
      probit_eta <- if (inherits(eta, "simref")) eta$value else eta
      if (inherits(probit_eta, "Matrix")) {
        probit_eta <- as.matrix(probit_eta)
      }
      RTMB::pnorm(probit_eta)
    },
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
  }

  ## Dispersion linear predictor; adapted from
  ## glmmTMB.cpp:839, 926-932, and 939
  sparseXdisp <- nrow(Xdisp) == 0 && ncol(Xdisp) == 0
  Xdispc <- if (sparseXdisp) XdispS else Xdisp
  etadisp <- Xdispc %*% betadisp + Zdisp %*% bdisp + dispoffset
  phi <- exp(etadisp)

  ## Observation likelihoods; adapted from glmmTMB.cpp:961-978,
  ## 1095-1101, and 1180-1199
  i <- !is.na(yobs_obs) | inherits(yobs, "simref")
  yobs_i <- yobs[i]
  keep <- osa_keep(yobs_i)
  eta_zi <- if (has_zi) etazi[i] else NULL
  logit_mu <- if (names(family) == "binomial") {
    logit_inverse_linkfun_rtmb(eta, link)
  } else {
    NULL
  }
  log_mu <- if (names(family) %in% c("nbinom1", "nbinom2")) {
    log_inverse_linkfun_rtmb(eta, link)
  } else {
    NULL
  }
  log_var_minus_mu <- switch(
    names(family),
    nbinom1 = log_mu + etadisp,
    nbinom2 = 2 * log_mu - etadisp,
    NULL
  )

  tmp_loglik <- switch(
    names(family),
    poisson = dZI(RTMB::dpois)(yobs_i, lambda = mu[i], eta_zi = eta_zi, log = TRUE,
                               is_zero = yobs_obs[i] == 0),
    truncated_poisson = dZI(dtruncated_poisson_rtmb)(
      yobs_i, lambda = mu[i], eta_zi = eta_zi, log = TRUE,
      is_zero = yobs_obs[i] == 0
    ),
    gaussian = dZI(dnorm_tmb)(yobs_i, mean = mu[i], sd = phi[i], eta_zi = eta_zi,
                              log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the binomial_family case in glmmTMB.cpp:979-983.
    binomial = dZI(dbinom_robust_rtmb)(
      yobs_i, size = size[i], logit_p = logit_mu[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the nbinom1_family case in glmmTMB.cpp:1042-1056.
    nbinom1 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu[i], log_var_minus_mu = log_var_minus_mu[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the nbinom2_family case in glmmTMB.cpp:1066-1075.
    nbinom2 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu[i], log_var_minus_mu = log_var_minus_mu[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    stop(
      "family not yet implemented: ", names(family),
      "; implemented families are: poisson, truncated_poisson, gaussian, ",
      "binomial, nbinom1, nbinom2"
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
    mu_pred_all <- phi
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
        nll <- nll - sum(dnorm_tmb(U[, j], 0, 1, log = TRUE))
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
