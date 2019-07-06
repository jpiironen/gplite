

#' @rdname gp_fit
#' @export
gp_sample <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  if (gp$method == 'full') {
    gp <- gp_mcmc_full(gp, x, y, trials=trials, jitter=jitter, ...)
  } else if (gp$method == 'rf') {
    num_inputs <- NCOL(x)
    featuremap <- rf_featmap(gp$cfs, num_inputs, gp$num_basis, seed=gp$rf_seed)
    gp <- gp_mcmc_linearized(gp, x, y, featuremap, trials=trials, ...)
  } else {
    stop('Unknown method: ', gp$method)
  }
  return(gp)
}


gp_mcmc_full <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cf, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(list(n=n,L=gp$K_chol,y=y), get_standata(gp$lik, trials=trials))
  model <- get_stanmodel(gp$lik, gp$method)
  gp$fit <- rstan::sampling(model, data=data, ...)
  gp$fsample <- t(rstan::extract(gp$fit)$f)
  return(gp)
}

gp_mcmc_linearized <- function(gp, x, y, featuremap, trials=NULL, ...) {
  gp$featuremap <- featuremap
  z <- gp$featuremap(x)
  n <- length(y)
  data <- c(list(n=n,m=ncol(z),Z=z,y=y), get_standata(gp$lik, trials=trials))
  model <- get_stanmodel(gp$lik, gp$method)
  gp$fit <- rstan::sampling(model, data=data, ...)
  gp$wsample <- t(rstan::extract(gp$fit)$w)
  return(gp)
}









