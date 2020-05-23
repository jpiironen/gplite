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
  gp <- gp_laplace(gp, x, y, trials=trials, jitter=jitter, ...)
  gp$fitted <- TRUE
  return(gp)
}

gp_laplace <- function(object, ...) {
  UseMethod("gp_laplace", object)
}

gp_laplace.gp <- function(object, ...) {
  gp_laplace(object$approx, object, ...)
}

gp_laplace.approx_full <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cfs, x, x) + jitter*diag(n)
  gp$x <- x
  gp$fit <- laplace(object, gp, K, y, trials=trials)
  return(gp)
}

gp_laplace.approx_fitc <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  set.seed(gp$approx$seed)
  # inducing points via k-means
  z <- get_inducing(gp, x)
  Kz <- eval_cf(gp$cfs, z, z) + jitter*diag(gp$approx$num_inducing)
  Kxz <- eval_cf(gp$cfs, x, z)
  Kz_chol <- t(chol(Kz))
  D <- rep(0,n)
  for (i in 1:n) {
    # TODO: this is slow
    D[i] <- eval_cf(gp$cfs, x[i,,drop=F], x[i,,drop=F])
  }
  D <- D - colSums(forwardsolve(Kz_chol, t(Kxz))^2)
  gp$x <- x
  gp$x_inducing <- z
  gp$fit <- tryCatch({
    laplace(object, gp, Kz, Kz_chol, Kxz, D, y, trials=trials)
  },
  error = function(err) {
    print(err)
    list(log_evidence=-Inf)
  })
  return(gp)
}

gp_laplace.approx_rf <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  num_inputs <- NCOL(x)
  featuremap <- get_featuremap(gp, num_inputs)
  gp$approx$num_basis <- check_num_basis(gp$cfs, gp$approx$num_basis, NCOL(x))
  z <- featuremap(x)
  gp$fit <- laplace(object, gp, z, y, trials=trials)
  return(gp)
}

gp_laplace.approx_rbf <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  gp$x <- x
  num_inputs <- NCOL(x)
  featuremap <- get_featuremap(gp, num_inputs)
  gp$approx$num_basis <- check_num_basis(gp$cfs, gp$approx$num_basis, NCOL(x))
  z <- featuremap(x)
  gp$fit <- laplace(object, gp, z, y, trials=trials)
  return(gp)
}


get_inducing <- function(gp, x) {
  if (!is.null(gp$x_inducing))
    return(gp$x_inducing)
  xscaled <- scale(x)
  cl <- kmeans(xscaled, gp$approx$num_inducing)
  zscaled <- cl$centers
  z <- t( t(zscaled)*attr(xscaled, 'scaled:scale') + attr(xscaled, 'scaled:center') )
  return(z)
}











