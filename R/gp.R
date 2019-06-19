
# main gp functions


# TODO: 
# - export main functions
# - add documentation
# - add function like gp_evidence(gp, param=NULL) which evaluates the log marginal likelihood
# at the given hyperparameter values, or at the current values if no new values are given










#' Initialize a GP model
#'
#' Initializes a GP model with given covariance function(s) and likelihood. The model can then be fitted using \code{\link{gp_fit}} or \code{\link{gp_sample}}. For hyperparameter optimization, see \code{\link{gp_optim}}
#' 
#' @param cf The covariance function(s). Either a single covariance function or a list of them. See \code{\link{cf}}.
#' @param lik Likelihood (observation model). See \code{\link{lik}}.
#'
#'
#' @return A GP model object that can be passed to other functions, for example when optimizing the hyperparameters or making predictions.
#' 
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' ### Basic usage (single covariance function)
#' cf <- gpcf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x ,y, trials)
#' 
#' }
#'

#' @export
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



