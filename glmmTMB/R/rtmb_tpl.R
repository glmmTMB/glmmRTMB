cmb <- function(f, d) function(p) f(p, d)

## Wrap a density function with zero-inflation likelihood and simulation support
dZI <- function(density) {
  force(density)

  function(x, ..., zi = NULL, log = FALSE, is_zero = NULL) {
    if (inherits(x, "simref") && !is.null(zi)) {
      prob_nonzero <- 1 / (1 + exp(zi))
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
      return(rep(0, length(x)))
    }

    loglik <- density(x, ..., log = TRUE)
    if (is.null(zi)) {
      return(if (log) loglik else exp(loglik))
    }

    if (is.null(is_zero)) {
      is_zero <- x == 0
    }

    log_pz <- -RTMB::logspace_add(0, -zi)
    log_1mpz <- -RTMB::logspace_add(0, zi)
    ans <- log_1mpz + loglik

    if (any(is_zero)) {
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
      return(RTMB::dnorm(x, mean, sd, log))
    }
    ans <- log_density(x, mean, sd)
    if (log) ans else exp(ans)
  }
})

## Variables injected into rtmb_tpl() by RTMB::getAll()
utils::globalVariables(c(
  "X", "XS", "Z", "offset", "terms", "family", "link", "weights",
  "beta", "b", "theta",
  "Xzi", "XziS", "Zzi", "zioffset", "termszi",
  "betazi", "bzi", "thetazi",
  "Xdisp", "XdispS", "Zdisp", "dispoffset", "termsdisp",
  "betadisp", "bdisp", "thetadisp"
))

rtmb_tpl <- function(parameters, data) {
  ## Keep the original response for NA and structural-zero checks; OBS() may
  ## replace yobs with a simulation reference
  yobs_obs <- data$yobs
  RTMB::getAll(data, parameters)
  yobs <- RTMB::OBS(yobs)

  nll <- 0

  ## Random-effects contribution; translated from glmmTMB.cpp:900-903
  cond_re <- allterms_nll(b, theta, terms)
  zi_re <- allterms_nll(bzi, thetazi, termszi)
  disp_re <- allterms_nll(bdisp, thetadisp, termsdisp)
  nll <- nll + cond_re$nll + zi_re$nll + disp_re$nll

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
    probit = RTMB::pnorm(eta),
    cloglog = 1 - exp(-exp(eta)),
    inverse = 1 / eta,
    lambertW = exp(eta) * exp(exp(eta)),
    stop("link not yet implemented: ", names(link))
  )

  ## Zero-inflation linear predictor; adapted from
  ## glmmTMB.cpp:836, 880, and 919-925
  # has_zi <- length(betazi) > 0 || length(bzi) > 0
  has_zi <- length(betazi) > 0
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

  ## Gaussian and Poisson observation likelihoods; adapted from
  ## glmmTMB.cpp:961-967, 975-978, and 1180-1199
  i <- !is.na(yobs_obs) | inherits(yobs, "simref")
  zi <- if (has_zi) etazi[i] else NULL

  tmp_loglik <- switch(
    names(family),
    poisson = dZI(RTMB::dpois)(yobs[i], lambda = mu[i], zi = zi, log = TRUE,
                               is_zero = yobs_obs[i] == 0),
    gaussian = dZI(dnorm_tmb)(yobs[i], mean = mu[i], sd = phi[i], zi = zi,
                              log = TRUE, is_zero = yobs_obs[i] == 0),
    stop("family not yet implemented: ", names(family))
  )

  nll <- nll - sum(weights[i] * tmp_loglik)
  corr <- cond_re$corr
  sd <- cond_re$sd
  corrzi <- zi_re$corr
  sdzi <- zi_re$sd
  corrdisp <- disp_re$corr
  sddisp <- disp_re$sd

  REPORT(corr)
  REPORT(sd)
  REPORT(corrzi)
  REPORT(sdzi)
  REPORT(corrdisp)
  REPORT(sddisp)
  REPORT(b)
  REPORT(bzi)
  REPORT(bdisp)

  nll
}
## Partition the concatenated random effects and covariance parameters by term
## Term slicing is translated from allterms_nll() in glmmTMB.cpp:803-826
allterms_nll <- function(u, theta, terms) {
  nll <- 0
  corr <- vector("list", length(terms))
  sd <- vector("list", length(terms))
  names(corr) <- names(terms)
  names(sd) <- names(terms)

  if (length(terms) == 0) {
    return(list(nll = nll, corr = corr, sd = sd))
  }

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

    upointer <- upointer + nr
    tpointer <- tpointer + term$blockNumTheta
  }

  list(nll = nll, corr = corr, sd = sd)
}

## Construct the correlation matrix used by TMB's
## density::UNSTRUCTURED_CORR_t. TMB fills the lower triangle row-wise,
## whereas matrix lower-triangle assignment in R fills it column-wise
## leading to a different theta ordering for dim >= 4
tmb_unstructured_corr <- function(n, theta) {
  "[<-" <- RTMB::ADoverload("[<-")

  expected <- n * (n - 1L) / 2L
  if (length(theta) != expected) {
    stop("Expected ", expected, " correlation parameters")
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
## termwise_nll(), glmmTMB.cpp:358-650
termwise_nll <- function(U, theta, term) {
  ## Preserve automatic differentiation when filling correlation matrices
  "[<-" <- RTMB::ADoverload("[<-")

  name <- names(term$blockCode)
  supported <- c(
    "diag", "homdiag", "us", "cs", "homcs", "toep", "homtoep",
    "ar1", "hetar1", "ou", "propto", "equalto"
  )

  if (!name %in% supported) {
    stop("covariance structure not yet implemented: ", name)
  }

  n <- term$blockSize
  reps <- term$blockReps
  dim(U) <- c(n, reps)

  ## Homogeneous structures use one standard-deviation parameter;
  ## heterogeneous structures use one parameter per term component.
  homogeneous <- c("homdiag", "homcs", "homtoep", "ar1", "ou")
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
        stop("OU time vector length must equal block size")
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
    stop("covariance structure not yet implemented: ", name)
  )

  ## Diagonal structures factor into univariate normal densities;
  ## correlated structures use a scaled multivariate normal density.
  if (density_structure == "diag") {
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
  list(nll = nll, corr = report_corr, sd = report_sd)
}
