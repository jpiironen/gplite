

#' @rdname gp_fit
#' @export
gp_mcmc <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  if (gp$method == 'full') {
    gp <- gp_mcmc_full(gp, x, y, trials=trials, jitter=jitter, ...)
  } else if (gp$method == 'rf') {
    gp <- gp_mcmc_linearized(gp, x, y, trials=trials, ...)
  } else {
    stop('Unknown method: ', gp$method)
  }
  gp$fitted <- TRUE
  return(gp)
}


gp_mcmc_full <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cfs, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(list(n=n,L=gp$K_chol,y=as.array(y)), get_standata(gp$lik, trials=trials))
  model <- get_stanmodel(gp$lik, gp$method)
  gp$fit <- rstan::sampling(model, data=data, ...)
  gp$fsample <- t(rstan::extract(gp$fit)$f)
  return(gp)
}

gp_mcmc_linearized <- function(gp, x, y, featuremap, trials=NULL, ...) {
  num_inputs <- NCOL(x)
  featuremap <- get_featuremap(gp, num_inputs)
  z <- featuremap(x)
  n <- length(y)
  data <- c(list(n=n,m=ncol(z),Z=z,y=as.array(y)), get_standata(gp$lik, trials=trials))
  model <- get_stanmodel(gp$lik, gp$method)
  gp$fit <- rstan::sampling(model, data=data, ...)
  gp$wsample <- t(rstan::extract(gp$fit)$w)
  return(gp)
}









