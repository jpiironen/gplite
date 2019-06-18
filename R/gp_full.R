
library('rstan')
#rstan_options(auto_write = TRUE)
stanmodels <- list(
   gaussian=stan_model('~/Aalto/codes/R/rcodes/experim/gp_hilbert/gp_full_gaussian.stan'),
   binomial_logit=stan_model('~/Aalto/codes/R/rcodes/experim/gp_hilbert/gp_full_binomial_logit.stan'),
   binomial_probit=stan_model('~/Aalto/codes/R/rcodes/experim/gp_hilbert/gp_full_binomial_probit.stan')
  )

# TODO: add matern kernel

get_param <- function (object, ...) {
  UseMethod("get_param", object)
}

get_param.cf_const <- function(object, ...) {
  log(object$magn)
}

get_param.cf_sexp <- function(object, ...) {
  log(c(object$lscale, object$magn))
}

get_param.lik_gaussian <- function(object, ...) {
  log(object$sigma)
}

get_param.lik_binomial <- function(object, ...) {
  c()
}

set_param <- function (object, ...) {
  UseMethod("set_param", object)
}

set_param.cf_const <- function(object, param, ...) {
  object$magn <- exp(param[1])
  object
}

set_param.cf_sexp <- function(object, param, ...) {
  object$lscale <- exp(param[1])
  object$magn <- exp(param[2])
  object
}

set_param.lik_gaussian <- function(object, param, ...) {
  object$sigma <- exp(param[1])
  object
}

set_param.lik_binomial <- function(object, param, ...) {
  object
}

get_standata <- function(object, ...) {
  UseMethod("get_standata", object)
}

get_standata.lik_gaussian <- function(object, ...) {
  list(sigma=object$sigma)
}

get_standata.lik_binomial <- function(object, ...) {
  args <- list(...)
  if (is.null(args$trials))
    stop('trials must be provided for the binomial likelihood.')
  list(trials=args$trials)
}

get_response <- function(object, ...) {
  UseMethod("get_response", object)
}

get_response.gp <- function(object, ...) {
  get_response(object$lik, ...)
}

get_response.lik_gaussian <- function(object, f, ...) {
  f
}

get_response.lik_binomial <- function(object, f, ...) {
  if (object$link == 'probit')
    return(pnorm(f))
  else
    return(1/(1+exp(-f)))
}


eval_cf <- function (object, ...) {
  UseMethod("eval_cf", object)
}

eval_cf.list <- function(object, x1, x2, ...) {
  K <- 0
  for (k in seq_along(object)) {
    K <- K + eval_cf(object[[k]], x1, x2, ...)
  }
  K
}

eval_cf.cf_sexp <- function(object, x1, x2, ...) {
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  n1 <- NROW(x1)
  n2 <- NROW(x2)
  K <- matrix(nrow=n1,ncol=n2)
  lscale <- object$lscale
  magn <- object$magn
  ind <- object$ind
  if (is.null(ind))
    ind <- c(1:NCOL(x1))
    
  for (i in 1:n1) {
    for (j in 1:n2) {
      K[i,j] <- magn^2*exp(-sum((x1[i,]-x2[j,])^2/lscale^2))
    }
  }
  return(K)
}

eval_cf.cf_const <- function(object, x1, x2, ...) {
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  n1 <- NROW(x1)
  n2 <- NROW(x2)
  K <- matrix(object$magn^2, nrow=n1,ncol=n2)
  return(K)
}


gpcf_const <- function(magn=1.0) {
  cf <- list()
  cf$magn <- magn
  class(cf) <- 'cf_const'
  cf
}

gpcf_sexp <- function(ind=NULL, lscale=0.5, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_sexp'
  cf
}

lik_gaussian <- function(sigma=1.0) {
  lik <- list()
  lik$sigma <- sigma
  lik$stanmodel <- stanmodels$gaussian
  class(lik) <- 'lik_gaussian'
  lik
}

lik_binomial <- function(link='logit') {
  lik <- list()
  lik$link <- link
  if (link == 'probit')
    lik$stanmodel <- stanmodels$binomial_probit
  else
    lik$stanmodel <- stanmodels$binomial_logit
  class(lik) <- 'lik_binomial'
  lik
}


gp_init <- function(cfs=gpcf_sexp(), lik=lik_gaussian()) {
  gp <- list()
  if (class(cfs) != 'list')
    cfs <- list(cfs)
  gp$cfs <- cfs
  gp$lik <- lik
  class(gp) <- 'gp'
  gp
}

get_param.gp <- function(object, ...) {
  param <- c()
  # covariance parameters
  for (k in seq_along(object$cfs))
    param <- c(param, get_param(object$cfs[[k]]))
  # likelihood parameters
  param <- c(param, get_param(object$lik))
  param
}

set_param.gp <- function(object, param, ...) {
  # covariance parameters
  j <- 1
  for (k in seq_along(object$cfs)) {
    np <- length(get_param(object$cfs[[k]]))
    object$cfs[[k]] <- set_param(object$cfs[[k]], param[j:(j+np)])
    j <- j + np
  }
  # likelihood parameters
  object$lik <- set_param(object$lik, param[j:length(param)])
  object
}



gp_fit <- function(gp, x,y, trials=NULL, jitter=1e-3, ...) {
  x <- as.matrix(x)
  n <- length(y)
  K <- eval_cf(gp$cf, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(list(n=n,K=K,y=y), get_standata(gp$lik, trials=trials))
  gp$fit <- optimizing(gp$lik$stanmodel, data=data, hessian=T, as_vector=F, ...)
  gp$fmean <- gp$fit$par$f
  gp$fprec_chol <- t(chol(-as.matrix(gp$fit$hessian))) # cholesky of precision
  gp$log_evidence <- gp$fit$value + 0.5*n*log(2*pi) - sum(log(diag(gp$fprec_chol)))
  return(gp)
}


gp_sample <- function(gp, x,y, trials=NULL, jitter=1e-3, ...) {
  x <- as.matrix(x)
  n <- length(y)
  K <- eval_cf(gp$cf, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(list(n=n,K=K,y=y), get_standata(gp$lik, trials=trials))
  gp$fit <- sampling(gp$lik$stanmodel, data=data, ...)
  gp$fsample <- t(extract(gp$fit)$f)
  gp
}


gp_optim <- function(gp, x,y, method='Nelder-Mead', tol=1e-3, verbose=T, ...) {
  energy <- function(param) {
    gp <- set_param(gp, param)
    gp <- gp_fit(gp,x,y, ...)
    if (verbose)
      print(paste0('Energy: ', -gp$log_evidence))
    -gp$log_evidence
  }
  param0 <- get_param(gp)
  res <- optim(param0, energy, method=method, control = list(reltol=tol))
  param <- res$par
  gp <- set_param(gp, param)
  gp <- gp_fit(gp,x,y, ...)
  gp
}



gp_pred <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=1e-3) {
  # TODO: implement predictions for the sampling gp
  
  if (is.null(gp$fsample)) {
    # predictions based on analytical approximation
    pred <- gp_pred_analytic(gp, xt, var=var, draws=draws, transform=transform, jitter=jitter)
  } else {
    # prediction based on (markov chain) monte carlo draws from the posterior
    pred <- gp_pred_sampling(gp, xt, draws=draws, transform=transform, jitter=jitter)
  }
  return(pred)
}


gp_pred_analytic <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=1e-3) {
  
  # compute the latent mean first
  K_chol <- gp$K_chol
  Kt <- eval_cf(gp$cf, xt, gp$x)
  Ktt <- eval_cf(gp$cf, xt, xt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, gp$fmean))
  pred_mean <- as.vector(pred_mean)
  nt <- length(pred_mean)
  
  if (var == T || !is.null(draws)) {
    # covariance of the latent function
    aux1 <- solve(K_chol, t(Kt))
    aux2 <- solve(gp$fprec_chol, solve(t(K_chol), solve(K_chol, t(Kt))))
    pred_cov <- Ktt - t(aux1) %*% aux1 + t(aux2) %*% aux2 + jitter*diag(nt)
    
    if (is.null(draws))
      return(list(mean=pred_mean, var=diag(pred_cov)))
    else {
      # sample from the predictive distribution
      sample <- t(chol(pred_cov)) %*% matrix(rnorm(nt*draws), nrow=nt) + pred_mean
      if (transform)
        sample <- get_response(gp$lik, sample)
      return(sample)
    }
  }
  return(pred_mean)
}

gp_pred_sampling <- function(gp, xt, draws=NULL, transform=T, jitter=1e-3) {
  
  draws <- NCOL(gp$fsample)
  nt <- NROW(xt)
  Kt <- eval_cf(gp$cf, xt, gp$x)
  Ktt <- eval_cf(gp$cf, xt, xt)
  K_chol <- gp$K_chol
  aux <- solve(K_chol, t(Kt))
  pred_cov <- Ktt - t(aux) %*% aux + jitter*diag(nt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, gp$fsample))
  sample <- t(chol(pred_cov)) %*% matrix(rnorm(nt*draws), nrow=nt) + pred_mean
  if (transform)
    sample <- get_response(gp$lik, sample)
  return(sample)
}



