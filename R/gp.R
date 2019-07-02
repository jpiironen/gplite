
# main gp functions





#' Initialize a GP model
#'
#' Initializes a GP model with given covariance function(s) and likelihood. The model can then be fitted using \code{\link{gp_fit}} or \code{\link{gp_sample}}. For hyperparameter optimization, see \code{\link{gp_optim}}
#' 
#' @param cfs The covariance function(s). Either a single covariance function or a list of them. See \code{\link{cf}}.
#' @param lik Likelihood (observation model). See \code{\link{lik}}.
#' @param method Method for approximating the covariance function, can be one of 'full' or 'rff'.
#' 'full' means that exact covariance function is used, and the inference will be on the n latent
#' function values (inference time scales cubicly in n). 'rff' uses random Fourier features for 
#' approximating the covariance function, which means the inference scales cubicly in the number of
#' approximating basis functions \code{num_basis}.
#' @param num_basis Number of basis functions in the covariance approximation for 'rff' and other
#' basis function methods.
#' @param rff_seed Seed for random Fourier features for reproducible results.
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
gp_init <- function(cfs=gpcf_sexp(), lik=lik_gaussian(), method='full', num_basis=100,
                    rff_seed=12345) {
  gp <- list()
  if (class(cfs) != 'list')
    cfs <- list(cfs)
  gp$cfs <- cfs
  gp$lik <- lik
  gp$method <- method
  gp$num_basis <- num_basis
  gp$rff_seed <- rff_seed
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
    ifelse(is.null(object$fmean) && is.null(object$wmean), F, T)
  else if (type=='sampling')
    ifelse(is.null(object$fsample) && is.null(object$wsample), F, T)
}


















