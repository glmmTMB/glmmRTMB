cmb <- function(f, d) function(p) f(p, d)

rtmb_tpl <- function(parameters, data) {
  ## Keep the original response for NA and structural-zero checks; OBS() may
  ## replace yobs with a simulation reference
  yobs_obs <- data$yobs
  RTMB::getAll(data, parameters) ##but R will complain about visible bindings...
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

  if (names(link) == "log") {
    mu <- exp(eta)
  } else if(names(link) == "identity"){
    mu <- eta
  } else {
    stop("not yet implemented")
  }

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
  if (names(family) == "gaussian") {
    sparseXdisp <- nrow(Xdisp) == 0 && ncol(Xdisp) == 0
    Xdispc <- if (sparseXdisp) XdispS else Xdisp
    etadisp <- Xdispc %*% betadisp + Zdisp %*% bdisp + dispoffset

    sigma <- exp(etadisp)
  }

  ## Gaussian and Poisson observation likelihoods; adapted from
  ## glmmTMB.cpp:961-967, 975-978, and 1196-1199
  for(j in seq_along(yobs)){
    if(!is.na(yobs_obs[j]) || inherits(yobs, "simref")){
      if (names(family) == "poisson") {
        tmp_loglik <- RTMB::dpois(yobs[j], mu[j], log=TRUE)
      } else if(names(family) == "gaussian"){
        if (inherits(yobs, "simref")) {
          tmp_loglik <- RTMB::dnorm(yobs[j], mu[j], sd=sigma[j], log=TRUE)
        } else {
          z <- (yobs[j] - mu[j]) / sigma[j]
          tmp_loglik <- -(log(sigma[j]) + 0.5 * log(2 * pi) + 0.5 * z * z)
        }
      } else {
        stop("not yet implemented")
      }

      if(has_zi){
        ## Compute log(p_zero) and log(1 - p_zero) directly from the logit-scale
        ## predictor then combine the structural-zero and conditional components
        ## Adapted from glmmTMB.cpp:1180-1193

        log_pz <- -RTMB::logspace_add(0, -etazi[j])
        log_1mpz <- -RTMB::logspace_add(0, etazi[j])
        if(yobs_obs[j] == 0) {
          tmp_loglik <- RTMB::logspace_add(log_pz, log_1mpz + tmp_loglik)
        } else {
          tmp_loglik <- log_1mpz + tmp_loglik
        }
      }
      nll <- nll - weights[j] * tmp_loglik
    }
  }
  corr <- cond_re$corr
  sd <- cond_re$sd
  corrzi <- zi_re$corr
  sdzi <- zi_re$sd
  corrdisp <- disp_re$corr
  sddisp <- disp_re$sd

  RTMB::REPORT(corr)
  RTMB::REPORT(sd)
  RTMB::REPORT(corrzi)
  RTMB::REPORT(sdzi)
  RTMB::REPORT(corrdisp)
  RTMB::REPORT(sddisp)
  RTMB::REPORT(b)
  RTMB::REPORT(bzi)
  RTMB::REPORT(bdisp)

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

## Evaluate one random-effects term under its covariance structure
## Partial translation of termwise_nll() in glmmTMB.cpp:353-800
termwise_nll <- function(U, theta, term) {
  ## Preserve automatic differentiation when filling matrices by assignment
  "[<-" <- RTMB::ADoverload("[<-")
  nll <- 0
  name <- names(term$blockCode)
  if (name == "us") {
    ## Unstructured covariance; translated from
    ## glmmTMB.cpp:407-416 and 436-440
    n <- term$blockSize
    reps <- term$blockReps
    logsd <- head(theta, n)
    corr_transf <- theta[-(seq_len(n))]
    sd <- exp(logsd)
    us <- RTMB::unstructured(n)
    C <- us$corr(corr_transf)
    dim(U) <- c(n, reps)
    nll <- nll - sum(RTMB::dmvnorm(t(U), Sigma=C, log=TRUE, scale=sd))
    return(list(nll = nll, corr = C, sd = sd))
  } else if (name == "diag") {
    ## Heterogeneous diagonal covariance; translated from
    ## glmmTMB.cpp:358-362 and 382
    n <- term$blockSize
    reps <- term$blockReps
    logsd <- head(theta, n)
    sd <- exp(logsd)
    dim(U) <- c(n, reps)
    for (k in seq_len(n)) {
      nll <- nll - sum(RTMB::dnorm(U[k, ], 0, sd[k], log = TRUE))
    }
    return(list(nll = nll, corr = matrix(numeric(0), 0, 0), sd = sd))
  } else if (name == "homdiag") {
    ## Homogeneous diagonal covariance; translated from
    ## glmmTMB.cpp:384-389 and 398-405
    n <- term$blockSize
    reps <- term$blockReps
    sd <- rep(exp(theta[1]), n)
    dim(U) <- c(n, reps)
    for (k in seq_len(n)) {
      nll <- nll - sum(RTMB::dnorm(U[k, ], 0, sd[k], log = TRUE))
    }
    return(list(nll = nll, corr = matrix(numeric(0), 0, 0), sd = sd))
  } else if (name == "cs" || name == "homcs") {
    ## Compound-symmetry covariance; translated from
    ## glmmTMB.cpp:441-463 and 471-473
    n <- term$blockSize
    reps <- term$blockReps
    if (name == "cs") {
      logsd <- theta[seq_len(n)]
      corr_transf <- theta[n + 1L]
    } else {
      logsd <- rep(theta[1], n)
      corr_transf <- theta[2]
    }

    sd <- exp(logsd)
    a <- 1 / (n-1)
    rho <- (1 / (1 + exp(-corr_transf))) * (1 + a) - a
    C <- diag(n)
    C[row(C) != col(C)] <- rho

    dim(U) <- c(n, reps)
    nll <- nll - sum(RTMB::dmvnorm(t(U), Sigma = C, log = TRUE, scale = sd))
    return(list(nll = nll, corr = C, sd = sd))

  } else if (name == "toep" || name == "homtoep") {
    ## Toeplitz covariance; translated from
    ## glmmTMB.cpp:474-496 and 504-506
    n <- term$blockSize
    reps <- term$blockReps
    if (name == "toep") {
      logsd <- theta[seq_len(n)]
      corr_raw <- theta[-seq_len(n)]
    } else {
      logsd <- rep(theta[1], n)
      corr_raw <- theta[-1]
    }

    sd <- exp(logsd)
    corr_params <- corr_raw / sqrt(1 + corr_raw^2)
    C <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        C[i, j] <- if (i == j) 1 else corr_params[abs(i - j)]
      }
    }

    dim(U) <- c(n, reps)
    nll <- nll - sum(RTMB::dmvnorm(t(U), Sigma = C, log = TRUE, scale = sd))

    return(list(nll = nll, corr = C, sd = sd))

  } else {
    stop("not yet implemented")
  }
}
