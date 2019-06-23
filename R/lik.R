
# implementations for the observation models (likelihoods)

# TODO: ADD DOCUMENTATION
# TODO: add beta binomial likelihood

#' Initialize likelihood
#'
#' Functions for initializing the likelihood (observations model) which can then be passed to \code{\link{gp_init}}.
#'
#' @name lik
#'
#' @param link Link function if the likelihood supports non-identity links.
#' @param sigma Initial value for the noise standard deviation.
#' 
#' @details Different likelihoods have the following possible link functions:
#' \describe{
#'  \item{\code{lik_gaussian}}{No links (uses identity link).}
#'  \item{\code{lik_binomial}}{'logit' or 'probit'}
#' }
#'
#' @return The likelihood object.
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
NULL

# constructors

#' @rdname lik
#' @export
lik_gaussian <- function(sigma=1.0) {
  lik <- list()
  lik$sigma <- sigma
  lik$stanmodel <- stanmodels$gaussian
  class(lik) <- 'lik_gaussian'
  lik
}

#' @rdname lik
#' @export
lik_binomial <- function(link='logit') {
  lik <- list()
  lik$link <- link
  if (link == 'probit')
    lik$stanmodel <- stanmodels$binomial_probit
  else
    lik$stanmodel <- stanmodels$binomial_logit
  class(lik) <- 'lik_binomial'
  lik
}



# get_param functions

get_param.lik_gaussian <- function(object, ...) {
  log(object$sigma)
}

get_param.lik_binomial <- function(object, ...) {
  c()
}



# set_param functions

set_param.lik_gaussian <- function(object, param, ...) {
  object$sigma <- exp(param[1])
  object
}

set_param.lik_binomial <- function(object, param, ...) {
  object
}




# get_standata functions

get_standata.lik_gaussian <- function(object, ...) {
  list(sigma=object$sigma)
}

get_standata.lik_binomial <- function(object, ...) {
  args <- list(...)
  if (is.null(args$trials))
    stop('trials must be provided for the binomial likelihood.')
  list(trials=args$trials)
}





# get_response functions

get_response.gp <- function(object, ...) {
  get_response(object$lik, ...)
}

get_response.lik_gaussian <- function(object, f, ...) {
  f
}

get_response.lik_binomial <- function(object, f, ...) {
  if (object$link == 'probit')
    return(pnorm(f))
  else if (object$link == 'logit')
    return(1/(1+exp(-f)))
  else
    stop('Unknown link: ', object$links)
}
