




#' Initialize a GP model
#'
#' Initializes a GP model with given covariance function(s) and likelihood. 
#' The model can then be fitted using \code{\link{gp_fit}}. For hyperparameter 
#' optimization, see \code{\link{gp_optim}}
#'
#' @param cfs The covariance function(s). Either a single covariance function 
#' or a list of them. See \code{\link{cf}}.
#' @param lik Likelihood (observation model). See \code{\link{lik}}.
#' @param method Method for approximating the covariance function.
#' See \code{\link{method}}.
#' @param approx Approximate inference method for Gaussian approximation
#' for the posterior of the latent values. See \code{\link{approx}}.
#'
#' @return A GP model object that can be passed to other functions, for example 
#' when optimizing the hyperparameters or making predictions.
#'
#'
#'
#'
#' @section References:
#'
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning.
#' MIT Press.
#'
#'
#' @examples
#'
#' # Full exact GP with Gaussian likelihood
#' gp <- gp_init(
#'   cfs = cf_sexp(),
#'   lik = lik_gaussian(),
#'   method = method_full()
#' )
#' 
#' # Binary classification model with EP approximation for the latent values
#' # and FITC sparse approximation to facilitate large datasets
#' gp <- gp_init(
#'   cfs = cf_sexp(),
#'   lik = lik_bernoulli(),
#'   approx = approx_ep(),
#'   method = method_fitc(num_inducing = 100)
#' )
#' 
#'
#' @export
gp_init <- function(cfs = cf_sexp(), lik = lik_gaussian(), method = method_full(),
                    approx = approx_laplace()) {
  gp <- list()
  if (!("list" %in% class(cfs))) {
    cfs <- list(cfs)
  }
  gp$cfs <- cfs
  gp$lik <- lik
  gp$method <- method
  if (!is.null(gp$method$num_basis)) {
    gp$method$num_basis <- check_num_basis(cfs, gp$method$num_basis)
  }
  gp$approx <- approx
  gp$fitted <- FALSE
  class(gp) <- "gp"
  gp
}


#' @export
print.gp <- function(x, ...) {
  object <- x
  indent <- "  "
  indent_str <- function(s) paste(indent, strsplit(s, "\n")[[1]], "\n", sep = " ", collapse = "")

  # likelihood
  str_lik <- indent_str(print(object$lik, quiet = TRUE))
  str <- paste0("Likelihood:\n", str_lik)

  # cfs
  str <- paste0(str, "Covariance functions:\n")
  for (i in seq_along(object$cfs)) {
    str_cf <- indent_str(print(object$cfs[[i]], quiet = TRUE))
    str <- paste0(str, str_cf)
  }

  # method
  str <- paste0(str, "Method:\n")
  str <- paste0(str, indent_str(print(object$method, quiet = TRUE)))

  # approx
  str <- paste0(str, "Latent approximation:\n")
  str <- paste0(str, indent_str(print(object$approx, quiet = TRUE)))

  cat(str)
  invisible(str)
}


#' Energy of a GP model
#'
#' Returns the energy (negative log marginal likelihood) of a fitted GP model with the
#' current hyperparameters. The result is exact for the Gaussian likelihood and
#' dependent on the \code{\link{approx}} for other cases.
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
#' x <- matrix(rnorm(n * 3), nrow = n)
#' f <- sin(x[, 1]) + 0.5 * x[, 2]^2 + x[, 3]
#' y <- f + 0.5 * rnorm(n)
#' x <- data.frame(x1 = x[, 1], x2 = x[, 2], x3 = x[, 3])
#'
#' # Basic usage
#' gp <- gp_init(cf_sexp(), lik_gaussian())
#' gp <- gp_fit(gp, x, y)
#' e <- gp_energy(gp)
#' }
#'
#' @export
gp_energy <- function(gp, include_prior = TRUE) {
  if (!is_fitted(gp, type = "analytic")) {
    stop("The GP must be fitted. Call gp_fit first.")
  }
  energy <- -gp$fit$log_evidence
  if (include_prior) {
    energy <- energy - lpdf_prior(gp)
  }
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
  object$lik <- set_param(object$lik, utils::tail(param, np_lik))
  object$fitted <- FALSE
  object
}

lpdf_prior.gp <- function(object, ...) {
  lpdf_prior(object$cfs) + lpdf_prior(object$lik)
}

get_featuremap.gp <- function(object, num_inputs, ...) {
  if ("method_rf" %in% class(object$method)) {
    featmap <- rf_featmap(object$cfs, object$method$num_basis,
      num_inputs = num_inputs, seed = object$method$seed
    )
    return(featmap)
  } else if ("method_rbf" %in% class(object$method)) {
    x <- object$x
    if (is.null(x)) {
      stop("Cannot compute RBF feature map when the x matrix is not known.")
    }
    featmap <- rbf_featmap(object$cfs, object$method$num_basis,
      num_inputs = num_inputs, x = x, seed = object$method$seed
    )
    return(featmap)
  } else {
    stop("No feature map implementation for method: ", object$method$name)
  }
}

is_fitted.gp <- function(object, type, ...) {
  if (is.list(object$fit)) {
    if (type == "analytic") {
      fit_found <- object$fit$type == "analytic"
    } else if (type == "sampling") {
      fit_found <- object$fit$type == "mcmc"
    }
  } else {
    fit_found <- FALSE
  }
  if (fit_found && object$fitted == FALSE) {
    stop("The GP object seems to contain a posterior fit, but is not refitted after setting new hyperparameter values. Please refit using gp_fit after calling set_param.")
  }
  return(fit_found)
}



# below are some functions for handling the linearized gp

check_num_basis <- function(cfs, num_basis, num_inputs = NA) {
  if (is.null(num_basis)) {
    return(NULL)
  }
  if (length(num_basis) == 1) {
    num_basis <- rep(num_basis, length(cfs))
  }
  if (length(num_basis) != length(cfs)) {
    stop("The length of num_basis must match to the number of covariance functions.")
  }
  for (k in seq_along(cfs)) {
    if (get_name(cfs[[k]]) == "cf_const") {
      num_basis[k] <- 1
    } else if (get_name(cfs[[k]]) == "cf_lin") {
      if (!is.null(cfs[[k]]$vars)) {
        num_basis[k] <- length(cfs[[k]]$vars)
      } else {
        num_basis[k] <- num_inputs
      }
    }
  }
  return(num_basis)
}

get_weight_inds <- function(gp, cfind = NULL) {
  if (is.null(cfind)) {
    cfind <- seq_along(gp$cfs)
  }
  end_points <- c(0, cumsum(gp$method$num_basis))
  inds <- c()
  for (k in cfind) {
    inds <- c(inds, (end_points[k] + 1):end_points[k + 1])
  }
  return(inds)
}

get_w_mean <- function(gp, cfind = NULL) {
  if (is.null(cfind)) {
    cfind <- seq_along(gp$cfs)
  }
  inds <- get_weight_inds(gp, cfind)
  w <- gp$fit$wmean[inds]
  return(w)
}

get_w_cov <- function(gp, cfind = NULL) {
  if (is.null(cfind)) {
    cfind <- seq_along(gp$cfs)
  }
  inds <- get_weight_inds(gp, cfind)
  return(gp$fit$wcov[inds, inds])
}


get_w_sample <- function(gp, cfind = NULL) {
  if (is.null(cfind)) {
    cfind <- seq_along(gp$cfs)
  }
  inds <- get_weight_inds(gp, cfind)
  return(gp$fit$wsample[inds, , drop = FALSE])
}


# function for determining the default amount of jitter on the covarince diagonal
# for different likelihoods
get_jitter <- function(gp, jitter) {
  if (!is.null(jitter)) {
    return(jitter)
  }
  return(1e-6)
}
