

#' @rdname gp_fit
#' @export
gp_mcmc <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  gp <- learn_scales(gp, x)
  gp <- fit_mcmc(gp, x, y, trials=trials, jitter=jitter, ...)
  gp$fit$type <- 'mcmc'
  gp$fitted <- TRUE
  return(gp)
}

fit_mcmc <- function(object, ...) {
  UseMethod('fit_mcmc', object)
}

fit_mcmc.gp <- function(object, ...) {
  fit_mcmc(object$approx, object, ...)
}

fit_mcmc.approx <- function(object, ...) {
  stop("MCMC for ", class(object)[1], " not implemented yet.")
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
  stanfit <- rstan::sampling(model, data=data, ...)
  draws <- rstan::extract(stanfit)
  gp$fit <- list(fsample=t(draws$f))
  return(gp)
}

fit_mcmc.approx_fitc <- function(object, gp, x,y, trials=NULL, jitter, ...) {
  x <- as.matrix(x)
  n <- length(y)
  z <- get_inducing(gp, x)
  jitter <- get_jitter(gp,jitter)
  Kz <- eval_cf(gp$cfs, z, z) + jitter*diag(gp$approx$num_inducing)
  Kxz <- eval_cf(gp$cfs, x, z)
  Kz_chol <- t(chol(Kz))
  Kxz_U_inv <- t(forwardsolve(Kz_chol, t(Kxz)))
  D <- eval_cf(gp$cfs, x, x, diag_only=T) + jitter
  D <- D - colSums(forwardsolve(Kz_chol, t(Kxz))^2)
  gp$x <- x
  gp$x_inducing <- z
  data <- c(
    list(n=n, m=NROW(z), Kxz_U_inv=Kxz_U_inv, D=D, y=as.array(y)), 
    get_standata(gp$lik, trials=trials, n=n)
  )
  model <- get_stanmodel(gp$lik, gp$approx$name)
  stanfit <- rstan::sampling(model, data=data, ...)
  draws <- rstan::extract(stanfit)
  u <- Kz_chol %*% t(draws$u_white)
  gp$fit <- list(usample=u, Kz_chol=Kz_chol)
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
  stanfit <- rstan::sampling(model, data=data, ...)
  draws <- rstan::extract(stanfit)
  gp$fit <- list(wsample=t(draws$w))
  return(gp)
}










