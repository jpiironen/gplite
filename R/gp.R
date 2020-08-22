




#' Initialize a GP model
#'
#' Initializes a GP model with given covariance function(s) and likelihood. The model can then be fitted using \code{\link{gp_fit}} or \code{\link{gp_mcmc}}. For hyperparameter optimization, see \code{\link{gp_optim}}
#' 
#' @param cfs The covariance function(s). Either a single covariance function or a list of them. See \code{\link{cf}}.
#' @param lik Likelihood (observation model). See \code{\link{lik}}.
#' @param method Method for approximating the covariance function, can be one of 
#' \code{'full'} (full exact GP),
#' \code{'rf'} (random features), or 
#' \code{'fitc'} (fully independent training conditional, FITC). See below for details.
#' @param num_basis Number of basis functions in the covariance approximation for 'rf' and other
#' basis function methods.
#' @param num_inducing Number of inducing points for FITC approximation.
#' @param seed Seed for reproducible results.
#' 
#' @return A GP model object that can be passed to other functions, for example when optimizing the hyperparameters or making predictions.
#' 
#' @details The argument \code{method} defines the method for approximating the covariance
#' function calculation. The choices are:
#' \itemize{
#'  \item{\code{'full'}:} {
#'  Full exact covariance function is used, meaning
#' that the inference will be for the \code{n} latent
#' function values (fitting time scales cubicly in \code{n}).
#' }
#'  \item{\code{'rf'}:} {
#'  Uses random features 
#' (or basis functions) for approximating the covariance function, which means the inference
#' time scales cubicly in the number of approximating basis functions \code{num_basis}.
#' For stationary covariance functions random Fourier features (Rahimi and Recht, 2007) is used,
#' and for non-stationary kernels using case specific method when possible (for example, drawing
#' the hidden layer parameters randomly for \code{cf_nn}).
#' }
#'  \item{\code{'fitc'}:} {
#'  Uses the fully independent training conditional, FITC, approximation 
#'  (see Quiñonero-Candela and Rasmussen, 2005; Snelson and Ghahramani, 2006). 
#'  The fitting time scales O(n*m^2), where n is the number of data points and 
#'  m the number of inducing points \code{num_inducing}.
#'  The inducing point locations are chosen using the k-means algorithm.
#'  }
#'  }
#' 
#' 
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. 
#' MIT Press.
#' 
#' Rahimi, A. and Recht, B. (2008). Random features for large-scale kernel machines. 
#' In Advances in Neural Information Processing Systems 20.
#' 
#' Quiñonero-Candela, J. and Rasmussen, C. E (2005). A unifying view of sparse approximate 
#' Gaussian process regression. Journal of Machine Learning Research 6:1939-1959.
#' 
#' Snelson, E. and Ghahramani, Z. (2006). Sparse Gaussian processes using pseudo-inputs. 
#' In Advances in Neural Information Processing Systems 18.
#'
#' @examples
#' \donttest{
#' 
#' # Generate some toy data
#' set.seed(1242)
#' n <- 500
#' x <- matrix(rnorm(n*3), nrow=n)
#' f <- sin(x[,1]) + 0.5*x[,2]^2 + x[,3]
#' y <- f + 0.5*rnorm(n)
#' x <- data.frame(x1=x[,1], x2=x[,2], x3=x[,3])
#' 
#' # Full GP
#' gp <- gp_init(cf_sexp())
#' gp <- gp_optim(gp, x, y)
#' 
#' # Approximate solution using random features
#' gp <- gp_init(cf_sexp(), method='rf', num_basis=300)
#' gp <- gp_optim(gp, x, y)
#' 
#' # Approximate solution using FITC
#' gp <- gp_init(cf_sexp(), method='fitc', num_inducing=100)
#' gp <- gp_optim(gp, x, y)
#' 
#' }
#'

#' @export
gp_init <- function(cfs=cf_sexp(), lik=lik_gaussian(), method='full', 
                    num_basis=100, num_inducing=100, seed=12345) {
  gp <- list()
  if (!('list' %in% class(cfs)))
    cfs <- list(cfs)
  gp$cfs <- cfs
  gp$lik <- lik
  gp$approx <- get_approx(method, seed=seed, num_basis=check_num_basis(cfs, num_basis), 
                          num_inducing=num_inducing)
  gp$fitted <- FALSE
  class(gp) <- 'gp'
  gp
}


#' Energy of a GP model
#'
#' Returns the energy (negative log marginal likelihood) of a fitted GP model with the current hyperparameters. The result is exact for the Gaussian likelihood and based on Laplace approximation for other cases.
#' 
#' @param gp The fitted GP model.
#' @param include_prior Whether to add log density of the prior to the result (in which case
#' the result is -(log marginal likelihood + log prior))
#' 
#' @return The energy value (negative log marginal likelihood).
#' 
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' 
#' # Generate some toy data
#' set.seed(1242)
#' n <- 500
#' x <- matrix(rnorm(n*3), nrow=n)
#' f <- sin(x[,1]) + 0.5*x[,2]^2 + x[,3]
#' y <- f + 0.5*rnorm(n)
#' x <- data.frame(x1=x[,1], x2=x[,2], x3=x[,3])
#' 
#' # Basic usage
#' gp <- gp_init(cf_sexp(), lik_gaussian())
#' gp <- gp_fit(gp, x, y)
#' e <- gp_energy(gp)
#' 
#' }
#'

#' @export
gp_energy <- function(gp, include_prior=T) {
  if (!is_fitted(gp, type='analytic'))
    stop('The GP must be fitted. Call gp_fit first.')
  energy <- -gp$fit$log_evidence
  if (include_prior)
    energy <- energy - lpdf_prior(gp)
  energy
}

#' @export
get_param.gp <- function(object, ...) {
  # pack the covariance function and likelihood parameters
  param <- get_param(object$cfs)
  param <- c(param, get_param(object$lik))
  param
}

#' @export
set_param.gp <- function(object, param, ...) {
  # unpack the covariance function and likelihood parameters
  object$cfs <- set_param(object$cfs, param)
  np_lik <- length(get_param(object$lik))
  object$lik <- set_param(object$lik, utils::tail(param,np_lik))
  object$fitted <- FALSE
  object
}

lpdf_prior.gp <- function(object, ...) {
  lpdf_prior(object$cfs) + lpdf_prior(object$lik)
}

get_featuremap.gp <- function(object, num_inputs, ...) {
  if ('approx_rf' %in% class(object$approx)) {
    featmap <- rf_featmap(object$cfs, object$approx$num_basis, 
                          num_inputs=num_inputs, seed=object$approx$seed)
    return(featmap)
  } else if ('approx_rbf' %in% class(object$approx)) {
    x <- object$x
    if (is.null(x))
      stop('Cannot compute RBF feature map when the x matrix is not known.')
    featmap <- rbf_featmap(object$cfs, object$approx$num_basis, 
                           num_inputs=num_inputs, x=x, seed=object$approx$seed)
    return(featmap)
  } else
    stop('No feature map implementation for method: ', object$approx$name)
}

is_fitted.gp <- function(object, type, ...) {
  if (type=='analytic')
    fit_found <- is.list(object$fit)
  else if (type=='sampling')
    fit_found <- ifelse(is.null(object$fsample) && is.null(object$wsample), F, T)
  if (fit_found && object$fitted==FALSE)
    stop('The GP object seems to contain a posterior fit, but is not refitted after setting new hyperparameter values. Please refit using gp_fit or gp_mcmc after calling set_param.')
  return(fit_found)
}



# below are some functions for handling the linearized gp

check_num_basis <- function(cfs, num_basis, num_inputs=NA) {
  if (is.null(num_basis))
    return(NULL)
  if (length(num_basis) == 1)
    num_basis <- rep(num_basis, length(cfs))
  if (length(num_basis) != length(cfs))
    stop('The length of num_basis must match to the number of covariance functions.')
  for (k in seq_along(cfs))
    if (get_name(cfs[[k]]) == 'cf_const')
      num_basis[k] <- 1
    else if (get_name(cfs[[k]]) == 'cf_lin')
      num_basis[k] <- num_inputs
  return(num_basis)
}

get_weight_inds <- function(gp, cfind=NULL) {
  if (is.null(cfind))
    cfind <- seq_along(gp$cfs)
  end_points <- c(0, cumsum(gp$approx$num_basis))
  inds <- c()
  for (k in cfind)
    inds <- c(inds, (end_points[k]+1):end_points[k+1])
  return(inds)
}

get_w_mean <- function(gp, cfind=NULL) {
  if (is.null(cfind))
    cfind <- seq_along(gp$cfs)
  inds <- get_weight_inds(gp, cfind)
  w <- gp$fit$wmean[inds]
  return(w)
}

get_w_cov <- function(gp, cfind=NULL) {
  if (is.null(cfind))
    cfind <- seq_along(gp$cfs)
  inds <- get_weight_inds(gp, cfind)
  return(gp$fit$wcov[inds,inds])
}


get_w_sample <- function(gp, cfind=NULL) {
  if (is.null(cfind))
    cfind <- seq_along(gp$cfs)
  inds <- get_weight_inds(gp, cfind)
  return(gp$wsample[inds,,drop=F])
}


# function for determining the default amount of jitter on the covarince diagonal
# for different likelihoods
get_jitter <- function(gp, jitter) {
  if (!is.null(jitter))
    return(jitter)
  return(1e-6)
}












