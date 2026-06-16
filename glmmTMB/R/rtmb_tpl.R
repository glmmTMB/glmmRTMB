cmb <- function(f, d) function(p) f(p, d)

rtmb_tpl <- function(parameters, data) {
  RTMB::getAll(data, parameters) ## but R will complain about visible bindings...
  yobs <- RTMB::OBS(yobs)

  nll <- 0

  ## Random effect likelihood (ignoring 'zi', 'disp')
  nll <- nll + allterms_nll(b, theta, terms)

  ## Linear predictor (ignoring 'zi', 'disp')
  sparseX <- nrow(X)==0 && ncol(X)==0
  if (sparseX) X <- XS
  eta <- X %*% beta + Z %*% b + offset

  ## Apply link
  if (names(link) == "log") {
    mu <- exp(eta)
  } else if(names(link) == "identity"){
    mu <- eta
  } else {
    stop("not yet implemented")
  }

  ## Data likelihood
  i <- !is.na(yobs) | inherits(yobs, "simref")
  if (names(family) == "poisson") {
    nll <- nll - sum(RTMB::dpois(yobs[i], mu[i], log=TRUE))
  } else if(names(family) == "gaussian"){
    #dispersion linear predictor (fixed effects only)
    sparseXdist <- nrow(Xdist)==0 && ncol(Xdist)==0
    if (sparseXdist) Xdist <- XdistS
    etadisp <- Xdist %*% betadisp + offsetdisp

    sigma <- exp(etadisp)    
    nll <- nll - sum(RTMB::dnorm(yobs[i], mu[i], sd=sigma[i], log=TRUE))
  } else {
    stop("not yet implemented")
  }

  nll
}

allterms_nll <- function(u, theta, terms) {
  nll <- 0
  if (length(terms) > 1) stop("not yet implemented")
  for (term in terms) { ## TODO: Get segments as allterms_nll
    useg <- seq_along(u)
    tseg <- seq_along(theta)
    nll <- nll + termwise_nll(u[useg], theta[tseg], term)
  }
  nll
}

termwise_nll <- function(U, theta, term) {
  nll <- 0
  name <- names(term$blockCode)
  if (name == "us") {
    ## Direct translation of glmmTMB.cpp:407-440
    n <- term$blockSize
    reps <- term$blockReps
    logsd <- head(theta, n)
    corr_transf <- tail(theta, -n)
    sd <- exp(logsd)
    us <- RTMB::unstructured(n)
    C <- us$corr(corr_transf)
    dim(U) <- c(n, reps)
    nll <- nll - sum(RTMB::dmvnorm(t(U), Sigma=C, log=TRUE, scale=sd))
  } else {
    stop("not yet implemented")
  }
  nll
}
