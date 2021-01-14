
# implementations for the observation models (likelihoods)



#' Initialize likelihood
#'
#' Functions for initializing the likelihood (observation model) which can then be passed to \code{\link{gp_init}}.
#'
#' The supported likelihoods are:
#' \describe{
#'  \item{\code{lik_gaussian}}{Gaussian likelihood. Has no links (uses identity link).}
#'  \item{\code{lik_bernoulli}}{Bernoulli likelihood. Possible links: 'logit' or 'probit'.}
#'  \item{\code{lik_binomial}}{Binomial likelihood. Possible links: 'logit' or 'probit'.}
#'  \item{\code{lik_betabinom}}{Beta binomial likelihood. Possible links: 'logit' or 'probit'.}
#'  \item{\code{lik_poisson}}{Poisson likelihood. Possible links: 'log'.}
#' }
#'
#'
#' @name lik
#'
#' @param link Link function if the likelihood supports non-identity links. See Details for
#' information about possible links for each likelihood.
#' @param sigma Initial value for the noise standard deviation.
#' @param phi The over dispersion parameter for beta binomial likelihood.
#' @param prior_sigma Prior for hyperparameter \code{sigma}. See \code{\link{priors}}.
#' @param prior_phi Prior for hyperparameter \code{phi}. See \code{\link{priors}}.
#'
#'
#' @return The likelihood object.
#'
#' @examples
#'
#' # Basic usage
#' cf <- cf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' 
#'
NULL

# constructors

#' @rdname lik
#' @export
lik_gaussian <- function(sigma = 0.5, prior_sigma = prior_logunif()) {
  lik <- list()
  lik$sigma <- sigma
  lik$priors <- list(sigma = prior_sigma)
  class(lik) <- c("lik_gaussian", "lik")
  lik
}

#' @rdname lik
#' @export
lik_bernoulli <- function(link = "logit") {
  lik <- list()
  lik$link <- link
  class(lik) <- c("lik_bernoulli", "lik")
  lik
}

#' @rdname lik
#' @export
lik_binomial <- function(link = "logit") {
  lik <- list()
  lik$link <- link
  class(lik) <- c("lik_binomial", "lik")
  lik
}

#' @rdname lik
#' @export
lik_betabinom <- function(link = "logit", phi = 1.0, prior_phi = prior_logunif()) {
  lik <- list()
  lik$phi <- phi
  lik$link <- link
  lik$priors <- list(phi = prior_phi)
  class(lik) <- c("lik_betabinom", "lik")
  lik
}

#' @rdname lik
#' @export
lik_poisson <- function(link = "log") {
  lik <- list()
  lik$link <- link
  class(lik) <- c("lik_poisson", "lik")
  lik
}



#' @export
print.lik <- function(x, quiet = FALSE, ...) {
  object <- x
  param_names <- get_param_names(object)
  param <- unlist(object[param_names])
  digits <- 3
  description <- paste0(get_name(object), "(")
  for (i in seq_along(param_names)) {
    description <- paste0(description, param_names[i], " = ", round(param[i], digits))
    if (i < length(param_names)) {
      description <- paste0(description, "; ")
    }
  }
  description <- paste0(description, ")\n")
  if (!quiet) {
    cat(description)
  }
  invisible(description)
}


get_name.lik <- function(object, ...) {
  class(object)[1]
}


# get_param_names functions

get_param_names.lik_gaussian <- function(object) {
  c("sigma")
}

get_param_names.lik_bernoulli <- function(object) {
  c()
}

get_param_names.lik_binomial <- function(object) {
  c()
}

get_param_names.lik_betabinom <- function(object) {
  c("phi")
}

get_param_names.lik_poisson <- function(object) {
  c()
}



# get_param functions

get_param.lik <- function(object, ...) {
  param_names <- filter_fixed(object, get_param_names(object))
  if (length(param_names) == 0) {
    return(NULL)
  }
  param <- unlist(object[param_names])
  names(param) <- add_obj_name(object, names(param))
  param <- log(param)
  param
}



# set_param functions

set_param.lik <- function(object, param, ...) {
  param_names <- filter_fixed(object, names(param))
  param_names <- rm_obj_name(object, param_names)
  for (j in seq_along(param_names)) {
    object[[param_names[j]]] <- unname(exp(param[j]))
  }
  object
}


# lpdf_prior function(s)

lpdf_prior.lik <- function(object, ...) {
  param <- get_param(object)
  param_names <- rm_obj_name(object, names(param))
  param_names <- filter_fixed(object, param_names)
  lp <- 0
  for (j in seq_along(param_names)) {
    lp <- lp + lpdf_prior(object$priors[[param_names[j]]], unname(param[j]))
  }
  lp
}



# get_pseudodata_la functions

get_pseudodata_la.lik_gaussian <- function(object, f, y, ...) {
  n <- length(y)
  loglik <- sum(stats::dnorm(y, mean = f, sd = object$sigma, log = TRUE))
  list(z = y, var = object$sigma^2 * rep(1, n), loglik = loglik)
}

get_pseudodata_la.lik <- function(object, f, y, min_abs_curvature = 1e-12, ...) {
  f <- as.vector(f)
  out <- get_loglik_d2(object, f, y, ...)
  grad <- out$grad
  grad2 <- out$grad2
  grad2 <- sign(grad2) * pmax(abs(grad2), min_abs_curvature)
  list(z = f - grad / grad2, var = -1 / grad2)
}



# get_pseudodata_ep functions

get_pseudodata_ep.lik <- function(object, post_mean, post_prec,
                                  z_old, P_old, y, quad_order = 11, damping = 0.5, ...) {
  tilted_moments <- get_tilted_moments(object, post_mean, post_prec,
    z_old, P_old, y,
    quad_order = quad_order, ...
  )
  mean_tilted <- tilted_moments$mean_tilted
  var_tilted <- tilted_moments$var_tilted
  mean_cavity <- tilted_moments$mean_cavity
  var_cavity <- tilted_moments$var_cavity
  log_C <- tilted_moments$log_normalizing_const
  quad_ok <- tilted_moments$quad_ok

  P_new <- 1 / var_tilted - 1 / var_cavity
  nu_new <- (mean_tilted / var_tilted - mean_cavity / var_cavity)

  # damping
  nu_old <- z_old * P_old
  nu_new <- damping * nu_new + (1 - damping) * nu_old
  P_new <- damping * P_new + (1 - damping) * P_old

  return(list(
    z = (1 / P_new) * nu_new, var = 1 / P_new, log_normalizing_const = log_C,
    mean_cavity = mean_cavity, var_cavity = var_cavity, quad_ok = quad_ok
  ))
}


# get_tilted_moments functions

get_tilted_moments.lik <- function(object, post_mean, post_prec,
                                   z_old, P_old, y, quad_order = 11, ...) {
  cavity_prec <- post_prec - P_old
  cavity_mean <- 1 / cavity_prec * (post_mean * post_prec - z_old * P_old)

  # compute moments of the tilted distributions using quadrature
  gh <- gauss_hermite_points_scaled(cavity_mean, sqrt(1 / cavity_prec), order = quad_order)
  fgrid <- gh$x
  weights <- gh$weights

  loglik <- get_loglik(object, fgrid, y, sum = FALSE, ...)
  lik <- exp(loglik) # this might be unstable, but unavoidable

  if (any(is.na(lik))) {
    # set NaNs to zero, but return a flag indicating that the result might be unreliable
    lik[is.na(lik)] <- 0
    loglik[is.na(loglik)] <- -Inf
    quad_ok <- FALSE
  } else {
    quad_ok <- TRUE
  }

  log_C <- apply(loglik, 1, logsumexp, weights = weights)
  C <- exp(log_C)
  mean_tilted <- (1 / C) * (fgrid * lik) %*% weights
  mean_tilted <- as.vector(mean_tilted)
  s2_tilted <- (1 / C) * ((fgrid - mean_tilted)^2 * lik) %*% weights
  s2_tilted <- as.vector(s2_tilted)

  list(
    mean_tilted = mean_tilted, var_tilted = s2_tilted, log_normalizing_const = log_C,
    mean_cavity = cavity_mean, var_cavity = 1 / cavity_prec, quad_ok = quad_ok
  )
}




#

required_extra_args.lik <- function(object, ...) {
  c()
}

required_extra_args.lik_binomial <- function(object, ...) {
  c("trials")
}

required_extra_args.lik_betabinom <- function(object, ...) {
  c("trials")
}

required_extra_args.lik_poisson <- function(object, ...) {
  c()
}



# get_loglik functions

get_loglik.lik <- function(object, f, y, ...) {
  stop(paste("No implementation for class", class(object)[1], "yet."))
}

get_loglik.lik_gaussian <- function(object, f, y, sum = TRUE, ...) {
  f <- add_offset(f, ...)
  loglik <- stats::dnorm(y, mean = f, sd = object$sigma, log = TRUE)
  if (sum) {
    return(sum(loglik))
  }
  return(loglik)
}

get_loglik.lik_bernoulli <- function(object, f, y, sum = TRUE, ...) {
  get_loglik.lik_binomial(object, f, y, sum = sum, trials = rep(1, length(y)), ...)
}

get_loglik.lik_binomial <- function(object, f, y, sum = TRUE, ...) {
  args <- list(...)
  if (is.null(args$trials)) {
    stop("trials must be provided for the binomial likelihood.")
  }

  f <- add_offset(f, ...)
  mu <- get_response(object, f)
  successes <- y
  trials <- args$trials
  loglik <- stats::dbinom(y, trials, mu, log = TRUE)

  if (sum) {
    return(sum(loglik))
  }
  return(loglik)
}

get_loglik.lik_betabinom <- function(object, f, y, sum = TRUE, ...) {
  args <- list(...)
  if (is.null(args$trials)) {
    stop("trials must be provided for the beta binomial likelihood.")
  }

  f <- add_offset(f, ...)
  mu <- get_response(object, f)
  a <- mu / object$phi
  b <- (1 - mu) / object$phi
  successes <- y
  trials <- args$trials

  term1 <- lgamma(trials + 1) - lgamma(successes + 1) - lgamma(trials - successes + 1)
  term2 <- lgamma(successes + a) + lgamma(trials - successes + b) - lgamma(trials + a + b)
  term3 <- lgamma(a + b) - lgamma(a) - lgamma(b)
  loglik <- term1 + term2 + term3

  if (sum) {
    return(sum(loglik))
  }
  return(loglik)
}

get_loglik.lik_poisson <- function(object, f, y, sum = TRUE, ...) {
  
  f <- add_offset(f, ...)
  mu <- get_response(object, f)
  loglik <- stats::dpois(y, mu, log=T)
  
  if (sum) {
    return(sum(loglik))
  }
  return(loglik)
}



# get_loglik_d functions (derivative of the log likelihood w.r.t f_i)

get_loglik_d.lik_bernoulli <- function(object, f, y, ...) {
  get_loglik_d.lik_binomial(object, f, y, trials = rep(1, length(y)), ...)
}

get_loglik_d.lik_binomial <- function(object, f, y, ...) {
  args <- list(...)
  if (is.null(args$trials)) {
    stop("trials must be provided for the binomial likelihood.")
  }

  f <- add_offset(f, ...)
  mu <- get_response(object, f)
  trials <- args$trials

  if (object$link == "probit") {
    dmu <- stats::dnorm(f)
  } else if (object$link == "logit") {
    dmu <- mu * (1 - mu)
  } else {
    stop("Unknown link: ", object$link)
  }

  # loglik = y*log(mu) + (trials-y)*log(1-mu)
  grad <- dmu * (y / mu - (trials - y) / (1 - mu))
  return(grad)
}


get_loglik_d.lik_betabinom <- function(object, f, y, ...) {
  args <- list(...)
  if (is.null(args$trials)) {
    stop("trials must be provided for the beta binomial likelihood.")
  }

  f <- add_offset(f, ...)
  mu <- get_response(object, f)
  phi <- object$phi
  a <- mu / phi
  b <- (1 - mu) / phi
  successes <- y
  trials <- args$trials

  if (object$link == "probit") {
    dmu <- stats::dnorm(f)
  } else if (object$link == "logit") {
    dmu <- mu * (1 - mu)
  } else {
    stop("Unknown link: ", object$link)
  }

  term2 <- (dmu / phi) * (digamma(successes + a) - digamma(trials - successes + b))
  term3 <- (dmu / phi) * (-digamma(a) + digamma(b))
  return(term2 + term3)
}

get_loglik_d.lik_poisson <- function(object, f, y, ...) {
  
  f <- add_offset(f, ...)
  mu <- get_response(object, f)
  
  if (object$link == "log") {
    dmu_df <- mu
  } else {
    stop("Unknown link: ", object$link)
  }
  
  # loglik = y*log(mu) - mu - log(y!)
  # mu = exp(f+offset)
  # dmu_df = mu
  grad <- (y/mu - 1)*dmu_df
  return(grad)
}



# get_loglik_d2 functions (second derivative of the log likelihood w.r.t f_i)

get_loglik_d2.lik <- function(object, f, y, eps = 1e-6, ...) {
  grad <- get_loglik_d(object, f, y, ...)
  grad_eps <- get_loglik_d(object, f + eps, y, ...)
  grad2 <- (grad_eps - grad) / eps
  list(grad = grad, grad2 = grad2)
}



# get_response functions

get_response.gp <- function(object, f, ...) {
  get_response(object$lik, f, ...)
}

get_response.lik_gaussian <- function(object, f, ...) {
  f
}

get_response.lik_bernoulli <- function(object, f, ...) {
  get_response.lik_binomial(object, f, ...)
}

get_response.lik_binomial <- function(object, f, ...) {
  if (object$link == "probit") {
    return(stats::pnorm(f))
  } else if (object$link == "logit") {
    return(1 / (1 + exp(-f)))
  } else {
    stop("Unknown link: ", object$link)
  }
}

get_response.lik_betabinom <- function(object, f, ...) {
  if (object$link == "probit") {
    return(stats::pnorm(f))
  } else if (object$link == "logit") {
    return(1 / (1 + exp(-f)))
  } else {
    stop("Unknown link: ", object$link)
  }
}

get_response.lik_poisson <- function(object, f, ...) {
  if (object$link == "log") {
    return(exp(f))
  } else {
    stop("Unknown link: ", object$link)
  }
}



# generate_target functions

generate_target.gp <- function(object, f, ...) {
  if (NCOL(f) > 1) {
    out <- apply(f, 2, function(f_i) generate_target(object$lik, f_i, ...))
  } else {
    out <- generate_target(object$lik, f, ...)
  }
  return(out)
}

generate_target.lik_gaussian <- function(object, f, ...) {
  stats::rnorm(length(f)) * object$sigma + f
}

generate_target.lik_bernoulli <- function(object, f, ...) {
  generate_target.lik_binomial(object, f, trials = rep(1, length(f)))
}

generate_target.lik_binomial <- function(object, f, trials, ...) {
  mu <- get_response(object, f)
  stats::rbinom(length(f), trials, prob = mu)
}

generate_target.lik_betabinom <- function(object, f, trials, ...) {
  mu <- get_response(object, f)
  a <- mu / object$phi
  b <- (1 - mu) / object$phi
  pr <- stats::rbeta(length(f), a, b)
  stats::rbinom(length(f), trials, prob = pr)
}

generate_target.lik_poisson <- function(object, f, ...) {
  mu <- get_response(object, f)
  stats::rpois(length(f), mu)
}



