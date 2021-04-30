
# implementations for the priors

#' Initialize prior for hyperparameter
#'
#' Functions for initializing hyperparameter priors which can then be passed
#' to \code{\link{gp_init}}. See section Details for the prior explanations.
#'
#' The supported priors are:
#' \describe{
#'  \item{\code{prior_fixed}}{ The hyperparameter is fixed to its initial value,
#'  and is not optimized by \code{gp_optim}. }
#'  \item{\code{prior_logunif}}{ Improper uniform prior on the log of the parameter. }
#'  \item{\code{prior_lognormal}}{ Log-normal prior (Gaussian prior on the logarithm of the parameter). }
#'  \item{\code{prior_half_t}}{ Half Student-t prior for a positive parameter. }
#' }
#'
#' @param df Degrees of freedom
#' @param loc Location parameter of the distribution
#' @param scale Scale parameter of the distribution
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
#'
#' # Quasi-periodic covariance function, with fixed period
#' cf1 <- cf_periodic(
#'   period = 5,
#'   prior_period = prior_fixed(),
#'   cf_base = cf_sexp(lscale = 2)
#' )
#' cf2 <- cf_sexp(lscale = 40)
#' cf <- cf1 * cf2
#' gp <- gp_init(cf)
#'
#' # draw from the prior
#' set.seed(104930)
#' xt <- seq(-10, 10, len = 500)
#' plot(xt, gp_draw(gp, xt), type = "l")
#' 
#'
NULL




#' @rdname priors
#' @export
prior_fixed <- function() {
  prior <- list()
  class(prior) <- c("prior_fixed", "prior")
  return(prior)
}


#' @rdname priors
#' @export
prior_logunif <- function() {
  prior <- list()
  class(prior) <- c("prior_logunif", "prior")
  return(prior)
}

#' @rdname priors
#' @export
prior_lognormal <- function(loc = 0, scale = 1) {
  prior <- list()
  prior$loc <- loc
  prior$scale <- scale
  class(prior) <- c("prior_lognormal", "prior")
  return(prior)
}


#' @rdname priors
#' @export
prior_half_t <- function(df = 1, scale = 1) {
  prior <- list()
  prior$df <- df
  prior$scale <- scale
  class(prior) <- c("prior_half_t", "prior")
  return(prior)
}



# lpdf_prior functions

lpdf_prior.prior_fixed <- function(object, param) {
  0
}

lpdf_prior.prior_logunif <- function(object, param) {
  0
}

lpdf_prior.prior_half_t <- function(object, param) {
  theta <- exp(param) # actual parameter, positively constrained
  logdet_jacobian <- param
  log(2) - log(object$scale) + stats::dt(theta / object$scale, df = object$df, log = TRUE) + logdet_jacobian
}

lpdf_prior.prior_lognormal <- function(object, param) {
  stats::dnorm(param, mean = object$loc, sd = object$scale, log = TRUE)
}
