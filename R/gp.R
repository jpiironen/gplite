
# main gp functions





#' Initialize a GP model
#'
#' Initializes a GP model with given covariance function(s) and likelihood. The model can then be fitted using \code{\link{gp_fit}} or \code{\link{gp_sample}}. For hyperparameter optimization, see \code{\link{gp_optim}}
#' 
#' @param cfs The covariance function(s). Either a single covariance function or a list of them. See \code{\link{cf}}.
#' @param lik Likelihood (observation model). See \code{\link{lik}}.
#' 
#' @return A GP model object that can be passed to other functions, for example when optimizing the hyperparameters or making predictions.
#' 
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' # Basic usage (single covariance function)
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


#' Energy of a GP model
#'
#' Returns the energy (negative log marginal likelihood) of a fitted GP model with the current hyperparameters. The result is exact for the Gaussian likelihood and based on Laplace approximation for other cases.
#' 
#' @param gp The fitted GP model.
#' 
#' @return The energy value (negative log marginal likelihood).
#' 
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' # Basic usage (single covariance function)
#' cf <- gpcf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' e <- gp_energy(gp)
#' 
#' }
#'

#' @export
gp_energy <- function(gp) {
  if (!is_fitted(gp, type='analytic'))
    stop('The GP must be fitted. Call gp_fit first.')
  -gp$log_evidence
}

#' @export
get_param.gp <- function(object, ...) {
  param <- c()
  # covariance parameters
  for (k in seq_along(object$cfs))
    param <- c(param, get_param(object$cfs[[k]]))
  # likelihood parameters
  param <- c(param, get_param(object$lik))
  param
}

#' @export
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

is_fitted.gp <- function(object, type, ...) {
  if (type=='analytic')
    ifelse(is.null(object$fmean), F, T)
  else if (type=='sampling')
    ifelse(is.null(object$fsample), F, T)
}




#' Fit a GP model
#'
#' Function \code{gp_fit} fits a GP model with the current hyperparameters. Notice that this function 
#' does not optimize the hyperparameters in any way, but only finds the Laplace approximation (or the analytical 
#' true posterior in the case of Gaussian likelihood) to the latent values. Function \code{gp_sample} 
#' draws from the posterior of the latent values given the current hyperparameter estimates.
#' 
#' @name gp_fit
#' 
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input dimension). 
#' Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param trials Vecton of length n giving the number of trials for each observation in binomial 
#' (and beta binomial) model.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability. Default is 1e-4 for Gaussian and 1e-2 for other likelihoods.
#' @param ... Further arguments to be passed to \link{rstan}'s function \link{optimizing} (if \code{gp_fit} or
#' \code{gp_optim} was called) or \link{sampling} (when \code{gp_sample} was called).
#'
#'
#' @return An updated GP model object.
#'  
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' # Analytic approximation
#' cf <- gpcf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_fit(gp, x, y)
#' 
#' # MCMC solution
#' gp <- gp_sample(gp, x, y, trials=trials, chains=2, iter=1000)
#' }
#'
NULL

#' @rdname gp_fit
#' @export
gp_fit <- function(gp, x,y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cf, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(list(n=n,K=K,y=y), get_standata(gp$lik, trials=trials))
  gp$fit <- rstan::optimizing(gp$lik$stanmodel, data=data, hessian=T, as_vector=F, init=0, ...)
  gp$fmean <- gp$fit$par$f
  gp$fprec_chol <- t(chol(-as.matrix(gp$fit$hessian))) # cholesky of precision
  gp$log_evidence <- gp$fit$value + 0.5*n*log(2*pi) - sum(log(diag(gp$fprec_chol)))
  return(gp)
}


#' @rdname gp_fit
#' @export
gp_sample <- function(gp, x,y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cf, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(list(n=n,K=K,y=y), get_standata(gp$lik, trials=trials))
  gp$fit <- rstan::sampling(gp$lik$stanmodel, data=data, ...)
  gp$fsample <- t(rstan::extract(gp$fit)$f)
  gp
}

#' Optimize hyperparameters of a GP model
#' 
#' This function can be used to optimize the hyperparameters of the model to the maximum marginal
#' likelihood solution (type-II maximum likelihood) based on the Laplace approximation.
#' 
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input dimension). 
#' Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param method Optimization method that will be passed to \link{optim} function.
#' @param tol Relative change in the objective function value below the optimization is terminated. 
#' @param verbose If TRUE, then some information about the progress of the optimization is printed to the console.
#' @param ... Further arguments to be passed to \link{gp_fit} that are needed in the fitting process, for
#' example \code{trials} in the case of binomial likelihood.
#'
#'
#' @return An updated GP model object.
#'  
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' # Basic usage (single covariance function)
#' cf <- gpcf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x, y, trials=trials)
#' 
#' }
#'
#'
#' 
#' @export
gp_optim <- function(gp, x,y, method='Nelder-Mead', tol=1e-3, verbose=T, ...) {
  energy <- function(param) {
    gp <- set_param(gp, param)
    gp <- gp_fit(gp,x,y, ...)
    if (verbose)
      print(paste0('Energy: ', -gp$log_evidence))
    -gp$log_evidence
  }
  param0 <- get_param(gp)
  res <- stats::optim(param0, energy, method=method, control = list(reltol=tol))
  param <- res$par
  gp <- set_param(gp, param)
  gp <- gp_fit(gp,x,y, ...)
  gp
}














#' Make predictions with a GP model
#' 
#' This function can be used to make predictions with a fitted model.
#' 
#' @param gp A fitted GP model object.
#' @param xt N-by-d matrix of input values (N is the number of test points and d the input dimension). 
#' Can also be a vector of length N if the model has only a single input.
#' @param var Whether to compute the predictive variances along with predictive mean
#' @param draws Number of draws to generate from the predictive distribution for the latent values. 
#' @param transform Whether to transform the draws of latent values to the same scale as the target y.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability. Default is 1e-4 for Gaussian and 1e-2 for other likelihoods.
#'
#'
#' @return If draws is given (sampling from the predictive distribution), then returns an N-by-draws
#' matrix of random draws from the predictive distribution. If draws is not set, then returns either only
#' the predictive means (if \code{var=F}) or a list with both the predictive means and variances (if \code{var=T}).
#'  
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' # fit GP
#' gp <- gp_init()
#' gp <- gp_optim(gp, x, y, trials=trials)
#' 
#' # analytic prediction
#' pred <- gp_pred(gp, xt, var=T)
#' pred$mean 
#' pred$var 
#' 
#' # draws based on the analytic posterior approximation
#' draws <- gp_pred(gp, xt, draws=1000) 
#' 
#' # fit using mcmc and predict then
#' gp <- gp_sample(gp, x, y, trials=trials, iter=1000, chains=2)
#' draws <- gp_pred(gp, xt, draws=1000)
#' }
#'
#'
#' @export
gp_pred <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=NULL) {
  
  if (!is_fitted(gp, 'sampling')) {
    # predictions based on analytical approximation
    if (!is_fitted(gp, 'analytic')) 
      stop('The given model appears not to be fitted yet.')
    pred <- gp_pred_analytic(gp, xt, var=var, draws=draws, transform=transform, jitter=jitter)
  } else {
    # prediction based on (markov chain) monte carlo draws from the posterior
    pred <- gp_pred_sampling(gp, xt, draws=draws, transform=transform, jitter=jitter)
  }
  return(pred)
}


gp_pred_analytic <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=NULL) {
  
  # compute the latent mean first
  K_chol <- gp$K_chol
  Kt <- eval_cf(gp$cf, xt, gp$x)
  Ktt <- eval_cf(gp$cf, xt, xt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, gp$fmean))
  pred_mean <- as.vector(pred_mean)
  nt <- length(pred_mean)
  jitter <- get_jitter(gp,jitter)
  
  if (var == T || !is.null(draws)) {
    # covariance of the latent function
    aux1 <- solve(K_chol, t(Kt))
    aux2 <- solve(gp$fprec_chol, solve(t(K_chol), solve(K_chol, t(Kt))))
    pred_cov <- Ktt - t(aux1) %*% aux1 + t(aux2) %*% aux2 + jitter*diag(nt)
    
    if (is.null(draws))
      return(list(mean=pred_mean, var=diag(pred_cov)))
    else {
      # sample from the predictive distribution
      sample <- t(chol(pred_cov)) %*% matrix(stats::rnorm(nt*draws), nrow=nt) + pred_mean
      if (transform)
        sample <- get_response(gp$lik, sample)
      return(sample)
    }
  }
  return(pred_mean)
}

gp_pred_sampling <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {
  
  draws <- NCOL(gp$fsample)
  nt <- NROW(xt)
  jitter <- get_jitter(gp,jitter)
  Kt <- eval_cf(gp$cf, xt, gp$x)
  Ktt <- eval_cf(gp$cf, xt, xt)
  K_chol <- gp$K_chol
  aux <- solve(K_chol, t(Kt))
  pred_cov <- Ktt - t(aux) %*% aux + jitter*diag(nt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, gp$fsample))
  sample <- t(chol(pred_cov)) %*% matrix(stats::rnorm(nt*draws), nrow=nt) + pred_mean
  if (transform)
    sample <- get_response(gp$lik, sample)
  return(sample)
}






