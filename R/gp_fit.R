#' Fit a GP model
#'
#' Function \code{gp_fit} fits a GP model with the current hyperparameters. 
#' Notice that this function does not optimize the hyperparameters in any way, 
#' but only finds the Laplace approximation (or the analytical 
#' true posterior in the case of Gaussian likelihood) to the latent values. 
#' Function \code{gp_mcmc} draws from the posterior of the latent values 
#' given the current hyperparameter estimates using MCMC. For optimizing the hyperparameter
#' values, see \code{gp_optim}.
#' 
#' @name gp_fit
#' 
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input dimension). 
#' Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param trials Vecton of length n giving the number of trials for each observation in binomial 
#' (and beta binomial) model.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability.
#'  Default is 1e-6.
#' @param ... Further arguments to be passed to \link{rstan}'s function 
#' \code{\link[rstan]{optimizing}} (if \code{gp_fit} was called) or 
#' \code{\link[rstan]{sampling}} (if \code{gp_mcmc} was called).
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
#' cf <- cf_sexp(lscale=0.3, magn=0.5)
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_fit(gp, x, y)
#' 
#' # MCMC solution
#' gpmc <- gp_mcmc(gp, x, y, trials=trials, chains=2, iter=1000)
#' }
#'
NULL

#' @rdname gp_fit
#' @export
gp_fit <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  if (gp$method == 'full') {
    gp <- gp_laplace_full(gp, x, y, trials=trials, jitter=jitter, ...)
  } else if (gp$method == 'rf') {
    gp <- gp_laplace_linearized(gp, x, y, trials=trials, jitter=jitter, ...)
  } else
    stop('Unknown method: ', gp$method)
  gp$fitted <- TRUE
  return(gp)
}


gp_laplace_full <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cfs, x, x) + jitter*diag(n)
  gp$x <- x
  gp$fit <- approx_laplace_iterated(gp, K, y, trials=trials)
  gp$log_evidence <- gp$fit$log_evidence # TODO: this if a bit awkward..
  return(gp)
}


gp_laplace_linearized <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  num_inputs <- NCOL(x)
  featuremap <- get_featuremap(gp, num_inputs)
  gp$num_basis <- check_num_basis(gp$cfs, gp$num_basis, NCOL(x))
  z <- featuremap(x)
  gp$fit <- approx_laplace_linearized_iterated(gp, z, y, trials=trials)
  gp$log_evidence <- gp$fit$log_evidence # TODO: this is UGLY
  return(gp)
}








