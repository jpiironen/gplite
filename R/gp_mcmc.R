

#' @rdname gp_fit
#' @export
gp_mcmc <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  gp <- learn_scales(gp, x)
  gp <- fit_mcmc(gp, x, y, trials=trials, jitter=jitter, ...)
  gp$fitted <- TRUE
  return(gp)
}

fit_mcmc <- function(object, ...) {
  UseMethod('fit_mcmc', object)
}

fit_mcmc.gp <- function(object, ...) {
  fit_mcmc(object$approx, object, ...)
}

fit_mcmc.approx_full <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cfs, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(
    list(n=n, L=gp$K_chol, y=as.array(y)), 
    get_standata(gp$lik, trials=trials, n=n)
  )
  model <- get_stanmodel(gp$lik, gp$approx$name)
  gp$fit <- rstan::sampling(model, data=data, ...)
  gp$fsample <- t(rstan::extract(gp$fit)$f)
  return(gp)
}

fit_mcmc.approx_rf <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  num_inputs <- NCOL(x)
  featuremap <- get_featuremap(gp, num_inputs)
  z <- featuremap(x)
  n <- length(y)
  data <- c(
    list(n=n, m=ncol(z), Z=z, y=as.array(y)), 
    get_standata(gp$lik, trials=trials, n=n)
  )
  model <- get_stanmodel(gp$lik, gp$approx$name)
  gp$fit <- rstan::sampling(model, data=data, ...)
  gp$wsample <- t(rstan::extract(gp$fit)$w)
  return(gp)
}









