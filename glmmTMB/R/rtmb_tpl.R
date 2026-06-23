cmb <- function(f, d) function(p) f(p, d)

logspace_add <- function(a, b) {
  log(exp(a) + exp(b))
}


rtmb_tpl <- function(parameters, data) {
  yobs_obs <- data$yobs
  RTMB::getAll(data, parameters) ## but R will complain about visible bindings...
  yobs <- RTMB::OBS(yobs)

  nll <- 0

  ## Random effect likelihood
  cond_re <- allterms_nll(b, theta, terms)
  zi_re <- allterms_nll(bzi, thetazi, termszi)
  disp_re <- allterms_nll(bdisp, thetadisp, termsdisp)
  nll <- nll + cond_re$nll + zi_re$nll + disp_re$nll
  
  ## Linear predictor (ignoring 'zi', 'disp')
  sparseX <- nrow(X)==0 && ncol(X)==0
  Xc <- if (sparseX) XS else X
  eta <- Xc %*% beta + Z %*% b + offset

  ## Apply link
  if (names(link) == "log") {
    mu <- exp(eta)
  } else if(names(link) == "identity"){
    mu <- eta
  } else {
    stop("not yet implemented")
  }

  ## ZI Linear Predictor
  #has_zi <- length(betazi) > 0 || length(bzi) > 0
  has_zi <- length(betazi) > 0
  if (has_zi) {
    sparseXzi <- nrow(Xzi)==0 && ncol(Xzi)==0
    Xzic <- if(sparseXzi) XziS else Xzi
    etazi <- Xzic %*% betazi + Zzi %*% bzi + zioffset
  }

  ## Dispersion Linear Predictor
  if (names(family) == "gaussian") {
    sparseXdisp <- nrow(Xdisp) == 0 && ncol(Xdisp) == 0
    Xdispc <- if (sparseXdisp) XdispS else Xdisp
    etadisp <- Xdispc %*% betadisp + Zdisp %*% bdisp + dispoffset

    sigma <- exp(etadisp)
  }

  ## Data likelihood
  for(j in seq_along(yobs)){
    if(!is.na(yobs_obs[j]) || inherits(yobs, "simref")){
      if (names(family) == "poisson") {
        tmp_loglik <- RTMB::dpois(yobs[j], mu[j], log=TRUE)
      } else if(names(family) == "gaussian"){
        tmp_loglik <- RTMB::dnorm(yobs[j], mu[j], sd=sigma[j], log=TRUE)
      } else {
        stop("not yet implemented")
      }

      if(has_zi){
        #observation is a structural zero
        log_pz <- -logspace_add(0, -etazi[j])
        #observation is not a structural zero; drawn from Normal(u,o)
        log_1mpz <- -logspace_add(0, etazi[j])
        if(yobs_obs[j] == 0) {
          tmp_loglik <- logspace_add(log_pz, log_1mpz + tmp_loglik)
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

termwise_nll <- function(U, theta, term) {
  nll <- 0
  name <- names(term$blockCode)
  if (name == "us") {
    ## Direct translation of glmmTMB.cpp:407-440
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
  } else if(name == "diag"){
    n <- term$blockSize
    reps <- term$blockReps
    logsd <- head(theta, n)
    sd <- exp(logsd)
    dim(U) <- c(n, reps)
    for (k in seq_len(n)) {
      nll <- nll - sum(RTMB::dnorm(U[k, ], 0, sd[k], log = TRUE))
    }
    return(list(nll = nll, corr = matrix(numeric(0), 0, 0), sd = sd))
  } else if(name == "homdiag"){
    n <- term$blockSize
    reps <- term$blockReps
    sd <- rep(exp(theta[1]), n)
    dim(U) <- c(n, reps)
    for (k in seq_len(n)) {
      nll <- nll - sum(RTMB::dnorm(U[k, ], 0, sd[k], log = TRUE))
    }
    return(list(nll = nll, corr = matrix(numeric(0), 0, 0), sd = sd))
  } else {
    stop("not yet implemented")
  }
}
