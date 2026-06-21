cmb <- function(f, d) function(p) f(p, d)

logspace_add <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

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

  ## ZI Linear Predictor
  has_zi <- length(betazi) > 0
  if (has_zi) {
    sparseXzi <- nrow(Xzi)==0 && ncol(Xzi)==0
    if (sparseXzi) Xzi <- XziS
    etazi <- Xzi %*% betazi + offsetzi
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

    if(!has_zi){
      nll <- nll - sum(RTMB::dnorm(yobs[i], mu[i], sd=sigma[i], log=TRUE))
    } else{
      is_zero <- yobs[i] == 0
      log_ll_cont <- RTMB::dnorm(yobs[i], mu[i], sd = sigma[i], log = TRUE)
      #observation is a structural zero
      log_pz <- -log1p(exp(-etazi[i]))
      #observation is not a structural zero; drawn from Normal(u,o)
      log_1mpz <- -log1p(exp( etazi[i]))
      ll <- ifelse(is_zero, logspace_add(log_pz, log_1mpz + log_ll_cont), log_1mpz + log_ll_cont)      
      nll <- nll - sum(ll)

    }    
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
    corr_transf <- theta[-(seq_len(n))]
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
