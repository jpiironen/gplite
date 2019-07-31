
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
lik_gaussian <- function(sigma=0.5, prior_sigma=prior_logunif()) {
  lik <- list()
  lik$sigma <- sigma
  lik$priors <- list(sigma=prior_sigma)
  class(lik) <- c('lik_gaussian', 'lik')
  lik
}

#' @rdname lik
#' @export
lik_binomial <- function(link='logit') {
  lik <- list()
  lik$link <- link
  class(lik) <- c('lik_binomial', 'lik')
  lik
}

#' @rdname lik
#' @export
lik_betabinom <- function(link='logit', phi=0.1, prior_phi=prior_logunif()) {
  lik <- list() 
  lik$phi <- phi
  lik$link <- link
  lik$priors <- list(phi=prior_phi)
  class(lik) <- c('lik_betabinom', 'lik')
  lik
}


get_name.lik <- function(object, ...) {
  class(object)[1]
}


# get_param_names functions

get_param_names.lik_gaussian <- function(object) {
  c('sigma')
}

get_param_names.lik_binomial <- function(object) {
  c()
}

get_param_names.lik_betabinom <- function(object) {
  c('phi')
}



# get_param functions

get_param.lik <- function(object, ...) {
  param_names <- filter_fixed(object, get_param_names(object))
  if (length(param_names) == 0)
    return(NULL)
  param <- unlist(object[param_names])
  names(param) <- add_obj_name(object, names(param))
  param <- log(param)
  param
}



# set_param functions

set_param.lik <- function(object, param, ...) {
  param_names <- filter_fixed(object, names(param))
  param_names <- rm_obj_name(object, param_names)
  for (j in seq_along(param_names))
    object[[param_names[j]]] <- unname(exp(param[j]))
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


# generate_target functions

generate_target.gp <- function(object, f, ...) {
  generate_target(object$lik, f, ...)
}

generate_target.lik_gaussian <- function(object, f, ...) {
  stats::rnorm(length(f))*object$sigma + f
}

generate_target.lik_binomial <- function(object, f, trials, ...) {
  mu <- get_response(object, f)
  stats::rbinom(length(f), trials, prob = mu)
}

generate_target.lik_betabinom <- function(object, f, trials, ...) {
  mu <- get_response(object, f)
  a <- mu/object$phi
  b <- (1-mu)/object$phi
  pr <- stats::rbeta(length(f), a, b)
  stats::rbinom(length(f), trials, prob = pr)
}


# function for determining the default amount of jitter on the covarince diagonal
# for different likelihoods
get_jitter <- function(gp, jitter) {
  if (!is.null(jitter))
    return(jitter)
  return(1e-6)
}
