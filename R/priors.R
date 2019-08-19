
# implementations for the priors

#' Initialize prior for hyperparameter
#'
#' Functions for initializing hyperparameter priors which can then be passed 
#' to \code{\link{gp_init}}. See section Details for the prior explanations.
#' 
#' The supported priors are:
#' \describe{
#'  \item{\code{prior_fixed}}{ The hyperparameter is fixed to its initial value, and is not optimized by 
#'                             \code{gp_optim}. }
#'  \item{\code{prior_logunif}}{ Improper uniform prior on the log of the parameter. }
#' }
#' 
#'
#' @name priors
#'
#' @return The hyperprior object.
#' 
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#' 
#' @examples
#' \donttest{
#' # Fix one of the hyperparameters
#' cf <- cf_nn(sigma0 = 0, prior_sigma0 = prior_fixed())
#' lik <- lik_gaussian()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x, y)
#' 
#' 
#' }
#'
NULL




#' @rdname priors
#' @export
prior_fixed <- function() {
  prior <- list()
  class(prior) <- 'prior_fixed'
  return(prior)
}


#' @rdname priors
#' @export
prior_logunif <- function() {
  prior <- list()
  class(prior) <- 'prior_logunif'
  return(prior)
}



