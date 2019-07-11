
# implementations for the observation models (likelihoods)



#' Initialize likelihood
#'
#' Functions for initializing the likelihood (observation model) which can then be passed to \code{\link{gp_init}}.
#' 
#' The supported likelihoods are:
#' \describe{
#'  \item{\code{lik_gaussian}}{Gaussian likelihood. Has no links (uses identity link).}
#'  \item{\code{lik_binomial}}{Binomial likelihood. Possible links: 'logit' or 'probit'.}
#'  \item{\code{lik_betabinom}}{Beta binomial likelihood. Possible links: 'logit' or 'probit'.}
#' }
#' 
#'
#' @name lik
#'
#' @param link Link function if the likelihood supports non-identity links. See Details for 
#' information about possible links for each likelihood.
#' @param sigma Initial value for the noise standard deviation.
#' @param phi The over dispersion parameter for beta binomial likelihood.
#' 
#'
#' @return The likelihood object.
#' 
#' @examples
#' \donttest{
#' # Basic usage 
#' cf <- cf_sexp()
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
lik_gaussian <- function(sigma=0.5) {
  lik <- list()
  lik$sigma <- sigma
  class(lik) <- 'lik_gaussian'
  lik
}

#' @rdname lik
#' @export
lik_binomial <- function(link='logit') {
  lik <- list()
  lik$link <- link
  class(lik) <- 'lik_binomial'
  lik
}

#' @rdname lik
#' @export
lik_betabinom <- function(link='logit', phi=0.1) {
  lik <- list()
  lik$phi <- phi
  lik$link <- link
  class(lik) <- 'lik_betabinom'
  lik
}



# get_param functions

get_param.lik_gaussian <- function(object, ...) {
  param <- log(object$sigma)
  names(param) <- c('lik_gaussian.sigma')
  param
}

get_param.lik_binomial <- function(object, ...) {
  c()
}

get_param.lik_betabinom <- function(object, ...) {
  param <- log(object$phi)
  names(param) <- c('lik_betabinom.phi')
  param
}



# set_param functions

set_param.lik_gaussian <- function(object, param, ...) {
  object$sigma <- exp(param[1])
  object
}

set_param.lik_binomial <- function(object, param, ...) {
  object
}

set_param.lik_betabinom <- function(object, param, ...) {
  object$phi <- exp(param[1])
  object
}



# get_stanmodel functions

get_stanmodel.lik_gaussian <- function(object, method, ...) {
  if (method == 'full') 
    return(stanmodels$gp_gaussian)
  else if (method == 'rf') 
    return(stanmodels$gpa_gaussian)
  else
    stop('Got an unknown method: ', method)
}

get_stanmodel.lik_binomial <- function(object, method, ...) {
  if (method == 'full') {
    return(stanmodels$gp_binomial)
  } else if (method == 'rf') {
    return(stanmodels$gpa_binomial)
  } else
    stop('Got an unknown method: ', method)
}

get_stanmodel.lik_betabinom <- function(object, method, ...) {
  if (method == 'full') {
    return(stanmodels$gp_betabinom)
  } else if (method == 'rf') {
    return(stanmodels$gpa_betabinom)
  } else
    stop('Got an unknown method: ', method)
}



# get_standata functions

get_standata.lik_gaussian <- function(object, ...) {
  list(sigma=object$sigma)
}

get_standata.lik_binomial <- function(object, ...) {
  args <- list(...)
  if (is.null(args$trials))
    stop('trials must be provided for the binomial likelihood.')
  if (object$link == 'logit')
    link <- 0
  else if (object$link == 'probit')
    link <- 1
  list(trials=args$trials, link=link)
}

get_standata.lik_betabinom <- function(object, ...) {
  args <- list(...)
  if (is.null(args$trials))
    stop('trials must be provided for the beta binomial likelihood.')
  if (object$link == 'logit')
    link <- 0
  else if (object$link == 'probit')
    link <- 1
  list(trials=args$trials, link=link, phi=object$phi)
}




# get_response functions

get_response.gp <- function(object, f, ...) {
  # TODO: CHANGE SO THAT THIS WILL BE USED
  get_response(object$lik, f, ...)
}

get_response.lik_gaussian <- function(object, f, ...) {
  f
}

get_response.lik_binomial <- function(object, f, ...) {
  if (object$link == 'probit')
    return(stats::pnorm(f))
  else if (object$link == 'logit')
    return(1/(1+exp(-f)))
  else
    stop('Unknown link: ', object$link)
}

get_response.lik_betabinom <- function(object, f, ...) {
  if (object$link == 'probit')
    return(stats::pnorm(f))
  else if (object$link == 'logit')
    return(1/(1+exp(-f)))
  else
    stop('Unknown link: ', object$link)
}



# function for determining the default amount of jitter on the covarince diagonal
# for different likelihoods
get_jitter <- function(gp, jitter) {
  if (!is.null(jitter))
    return(jitter)
  return(1e-4)
}
