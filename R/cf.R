
# implementations for the covariance functions




#' Initialize covariance function
#'
#' Functions for initializing the covariance functions which can then be passed
#' to \code{\link{gp_init}}. See section Details for explanation of what covariance
#' function is what.
#'
#' The supported covariance functions are (see Rasmussen and Williams, 2006):
#' \describe{
#'  \item{\code{cf_const}}{ Constant covariance function. Can be used to model the intercept. }
#'  \item{\code{cf_lin}}{ Linear covariance function. Produces linear functions. }
#'  \item{\code{cf_sexp}}{ Squared exponential (or exponentiated quadratic, or Gaussian) covariance function.}
#'  \item{\code{cf_matern32}}{ Matern nu=3/2 covariance function. }
#'  \item{\code{cf_matern52}}{ Matern nu=5/2 covariance function. }
#'  \item{\code{cf_nn}}{ Neural network covariance function. }
#'  \item{\code{cf_periodic}}{ Periodic covariance function. The periodicity is achieved by mapping the
#'  original inputs through sine and cosine functions, and then applying the base kernel in this new space.}
#'  \item{\code{cf_prod}}{ Product of two or more covariance functions. }
#' }
#'
#' @name cf
#'
#' @param vars Indices of the inputs which are taken into account when calculating this
#'  covariance. If the input matrix has named columns, can also be a vector of column names.
#'  Default is all the inputs.
#' @param normalize Whether to automatically scale and center the inputs for the given
#' covariance function. Can be useful for inputs with mean and variance far from 0 and 1, respectively.
#' @param lscale Initial value for the length-scale hyperparameter.
#' @param magn Initial value for the magnitude hyperparameter (depicts the magnitude of
#' the variation captured by the given covariance function).
#' @param sigma0 Prior std for the bias in the neural network covariance function.
#' @param sigma Prior std for the weights in the hidden layers of the neural network
#' covariance function.
#' @param period Period length for the periodic covariance function.
#' @param cf_base Base covariance function that is used to model the variability within each period
#' in periodic covariance function.
#' @param prior_magn Prior for hypeparameter \code{magn}. See \code{\link{priors}}.
#' @param prior_lscale Prior for hyperparameter \code{lscale}. See \code{\link{priors}}.
#' @param prior_sigma0 Prior for hyperparameter \code{sigma0}. See \code{\link{priors}}.
#' @param prior_sigma Prior for hyperparameter \code{sigma}. See \code{\link{priors}}.
#' @param prior_period Prior for hyperparameter \code{period}. See \code{\link{priors}}.
#' @param ... Meaning depends on context. For \code{cf_prod} pass in the covariance functions in the product.
#' @param cf1 Instance of a covariance function.
#' @param cf2 Instance of a covariance function.
#'
#' @return The covariance function object.
#'
#' @section References:
#'
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#'
#' # Generate some toy data
#' set.seed(1242)
#' n <- 50
#' x <- matrix(rnorm(n * 3), nrow = n)
#' f <- sin(x[, 1]) + 0.5 * x[, 2]^2 + x[, 3]
#' y <- f + 0.5 * rnorm(n)
#' x <- data.frame(x1 = x[, 1], x2 = x[, 2], x3 = x[, 3])
#' 
#' # Basic usage (single covariance function)
#' cf <- cf_sexp()
#' lik <- lik_gaussian()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x, y)
#' plot(gp_pred(gp, x)$mean, y)
#' 
#' # More than one covariance function; one for x1 and x2, and another one for x3
#' cf1 <- cf_sexp(c("x1", "x2"))
#' cf2 <- cf_lin("x3")
#' cfs <- list(cf1, cf2)
#' lik <- lik_gaussian()
#' gp <- gp_init(cfs, lik)
#' gp <- gp_optim(gp, x, y, maxiter = 500)
#' plot(gp_pred(gp, x)$mean, y)
#' plot(x[, 3], gp_pred(gp, x, cfind = 2)$mean) # plot effect w.r.t x3 only
#' 
#'
NULL


# constructors

#' @rdname cf
#' @export
cf_const <- function(magn = 1.0, prior_magn = prior_logunif()) {
  cf <- list()
  cf$magn <- magn
  cf$priors <- list(magn = prior_magn)
  class(cf) <- c("cf_const", "cf")
  cf
}

#' @rdname cf
#' @export
cf_lin <- function(vars = NULL, magn = 1.0, prior_magn = prior_logunif(), normalize = FALSE) {
  cf <- list()
  cf$vars <- vars
  cf$magn <- magn
  cf$priors <- list(magn = prior_magn)
  cf$normalize <- normalize
  class(cf) <- c("cf_lin", "cf")
  cf
}

#' @rdname cf
#' @export
cf_sexp <- function(vars = NULL, lscale = 0.3, magn = 1.0,
                    prior_lscale = prior_logunif(), prior_magn = prior_logunif(),
                    normalize = FALSE) {
  cf <- list()
  cf$vars <- vars
  cf$lscale <- lscale
  cf$magn <- magn
  cf$priors <- list(lscale = prior_lscale, magn = prior_magn)
  cf$normalize <- normalize
  class(cf) <- c("cf_sexp", "cf")
  cf
}

#' @rdname cf
#' @export
cf_matern32 <- function(vars = NULL, lscale = 0.3, magn = 1.0,
                        prior_lscale = prior_logunif(), prior_magn = prior_logunif(),
                        normalize = FALSE) {
  cf <- list()
  cf$vars <- vars
  cf$lscale <- lscale
  cf$magn <- magn
  cf$priors <- list(lscale = prior_lscale, magn = prior_magn)
  cf$normalize <- normalize
  class(cf) <- c("cf_matern32", "cf")
  cf
}

#' @rdname cf
#' @export
cf_matern52 <- function(vars = NULL, lscale = 0.3, magn = 1.0,
                        prior_lscale = prior_logunif(), prior_magn = prior_logunif(),
                        normalize = FALSE) {
  cf <- list()
  cf$vars <- vars
  cf$lscale <- lscale
  cf$magn <- magn
  cf$priors <- list(lscale = prior_lscale, magn = prior_magn)
  cf$normalize <- normalize
  class(cf) <- c("cf_matern52", "cf")
  cf
}

#' @rdname cf
#' @export
cf_nn <- function(vars = NULL, sigma0 = 1.0, sigma = 3.0, magn = 1.0,
                  prior_sigma0 = prior_half_t(), prior_sigma = prior_half_t(),
                  prior_magn = prior_logunif(), normalize = TRUE) {
  cf <- list()
  cf$vars <- vars
  cf$sigma0 <- sigma0
  cf$sigma <- sigma
  cf$magn <- magn
  cf$priors <- list(sigma0 = prior_sigma0, sigma = prior_sigma, magn = prior_magn)
  cf$normalize <- normalize
  class(cf) <- c("cf_nn", "cf")
  cf
}

#' @rdname cf
#' @export
cf_periodic <- function(vars = NULL, period = 1, cf_base = cf_sexp(), prior_period = prior_logunif()) {
  cf <- list()
  cf$vars <- vars
  cf$period <- period
  cf$base <- cf_base
  cf$base$normalize <- FALSE # ensure no normalization for the base kernel
  cf$base$vars <- NULL # ensure base kernel uses both of the transformed features
  cf$priors <- list(period = prior_period)
  class(cf) <- c("cf_periodic", "cf")
  cf
}

#' @rdname cf
#' @export
cf_prod <- function(...) {
  cf <- list()
  cf$cfs <- list(...)
  class(cf) <- c("cf_prod", "cf")
  cf
}

#' @rdname cf
#' @export
"*.cf" <- function(cf1, cf2) {
  if ("cf_prod" %in% class(cf1)) {
    if ("cf_prod" %in% class(cf2)) {
      cf1$cfs <- c(cf1$cfs, cf2$cfs)
    } else {
      cf1$cfs <- c(cf1$cfs, list(cf2))
    }
    return(cf1)
  }
  if ("cf_prod" %in% class(cf2)) {
    return(cf2 * cf1)
  }
  return(cf_prod(cf1, cf2))
}


#' @export
print.cf_prod <- function(x, quiet = FALSE, ...) {
  object <- x
  indent <- "  "
  str <- paste0(get_name(object), ":\n")
  for (i in seq_along(object$cfs)) {
    str <- paste0(str, indent, "cf", i, ": ", print(object$cfs[[i]], quiet = TRUE))
  }
  if (!quiet) {
    cat(str)
  }
  invisible(str)
}

#' @export
print.cf_periodic <- function(x, quiet = FALSE, ...) {
  object <- x
  str <- print.cf(object, quiet = TRUE)
  str <- strsplit(str, ")")[[1]][1] # remove ')' and newline in the end

  str_base <- print(object$base, quiet = TRUE)
  str_base <- strsplit(str_base, ")")[[1]][1] # remove ')' and newline in the end

  str <- paste(str, "; cf_base = ", str_base, "))\n", sep = "")

  if (!quiet) {
    cat(str)
  }
  invisible(str)
}

#' @export
print.cf <- function(x, quiet = FALSE, ...) {
  object <- x
  param_names <- get_param_names(object)
  param <- unlist(object[param_names])
  digits <- 3
  description <- paste0(get_name(object), "(")
  if (!is.null(object$vars)) {
    description <- paste0(description, "vars = ", paste0(object$vars, collapse = ","), "; ")
  }
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




# for figuring out the name of the cf conveniently

get_name.cf <- function(object, ...) {
  class(object)[1]
}



# get_param_names functions

get_param_names.cf_const <- function(object) {
  c("magn")
}

get_param_names.cf_lin <- function(object) {
  c("magn")
}

get_param_names.cf_sexp <- function(object) {
  c("lscale", "magn")
}

get_param_names.cf_matern32 <- function(object) {
  c("lscale", "magn")
}

get_param_names.cf_matern52 <- function(object) {
  c("lscale", "magn")
}

get_param_names.cf_nn <- function(object) {
  c("sigma0", "sigma", "magn")
}

get_param_names.cf_periodic <- function(object) {
  c("period")
}

get_param_names.cf_prod <- function(object) {
  sapply(object$cfs, function(cf) get_param_names(cf))
}



# get_param functions

get_param.list <- function(object, ...) {
  param <- c()
  for (k in seq_along(object)) {
    param <- c(param, get_param(object[[k]]))
  }
  param
}

get_param.cf <- function(object, ...) {
  param_names <- filter_fixed(object, get_param_names(object))
  if (length(param_names) == 0) {
    return(NULL)
  }
  param <- unlist(object[param_names])
  names(param) <- add_obj_name(object, names(param))
  param <- log(param)
  param
}

get_param.cf_periodic <- function(object, ...) {
  param <- get_param(object$base)
  if (!is_fixed(object, "period")) {
    param <- c(log(object$period), param)
    names(param)[1] <- "cf_periodic.period"
  }
  # overwrite the parameter names of the base kernel
  names(param) <- add_obj_name(object, rm_obj_name(object$base, names(param)))
  param
}

get_param.cf_prod <- function(object, ...) {
  get_param(object$cfs)
}



# set_param functions

set_param.list <- function(object, param, ...) {
  j <- 1
  for (k in seq_along(object)) {
    np <- length(get_param(object[[k]]))
    if (np > 0) {
      object[[k]] <- set_param(object[[k]], param[j:(j + np - 1)])
    }
    j <- j + np
  }
  object
}

set_param.cf <- function(object, param, ...) {
  if (is.null(names(param))) {
    stop(paste0(
      "Caught unnamed vector of parameter values; please provide a named vector ",
      "(similar to get_param(gp))."
    ))
  }
  param_names <- rm_obj_name(object, names(param))
  param_names <- filter_fixed(object, param_names)
  for (j in seq_along(param_names)) {
    object[[param_names[j]]] <- unname(exp(param[j]))
  }
  object
}

set_param.cf_periodic <- function(object, param, ...) {
  fixed_period <- is_fixed(object, "period")
  if (!fixed_period) {
    object$period <- exp(param[1])
  }
  object$base <- set_param(object$base, utils::tail(param, length(param) - !fixed_period))
  object
}

set_param.cf_prod <- function(object, param, ...) {
  object$cfs <- set_param(object$cfs, param)
  object
}



# eval_cf functions

eval_cf.list <- function(object, x1, x2, cfind = NULL, diag_only = FALSE, ...) {
  if (is.null(cfind)) {
    cfind <- seq_along(object)
  }
  K <- 0
  for (k in cfind) {
    K <- K + eval_cf(object[[k]], x1, x2, diag_only = diag_only, ...)
  }
  K
}

eval_cf.cf_const <- function(object, x1, x2, diag_only = FALSE, ...) {
  x1 <- prepare_inputmat(object, x1)
  x2 <- prepare_inputmat(object, x2)
  n1 <- NROW(x1)
  n2 <- NROW(x2)
  if (diag_only) {
    K <- rep(object$magn^2, min(n1, n2))
  } else {
    K <- matrix(object$magn^2, nrow = n1, ncol = n2)
  }
  return(K)
}

eval_cf.cf_lin <- function(object, x1, x2, diag_only = FALSE, ...) {
  x1 <- prepare_inputmat(object, x1)
  x2 <- prepare_inputmat(object, x2)
  K <- object$magn^2 * x1 %*% t(x2)
  if (diag_only) {
    # TODO: this not the most efficient thing to do
    K <- as.vector(diag(K))
  }
  return(K)
}

eval_cf.cf_sexp <- function(object, x1, x2, diag_only = FALSE, ...) {
  x1 <- prepare_inputmat(object, x1)
  x2 <- prepare_inputmat(object, x2)
  K <- cf_sexp_c(x1, x2, object$lscale, object$magn, diag_only = diag_only)
  if (diag_only) {
    K <- as.vector(K)
  }
  return(K)
}

eval_cf.cf_matern32 <- function(object, x1, x2, diag_only = FALSE, ...) {
  x1 <- prepare_inputmat(object, x1)
  x2 <- prepare_inputmat(object, x2)
  K <- cf_matern32_c(x1, x2, object$lscale, object$magn, diag_only = diag_only)
  if (diag_only) {
    K <- as.vector(K)
  }
  return(K)
}

eval_cf.cf_matern52 <- function(object, x1, x2, diag_only = FALSE, ...) {
  x1 <- prepare_inputmat(object, x1)
  x2 <- prepare_inputmat(object, x2)
  K <- cf_matern52_c(x1, x2, object$lscale, object$magn, diag_only = diag_only)
  if (diag_only) {
    K <- as.vector(K)
  }
  return(K)
}

eval_cf.cf_nn <- function(object, x1, x2, diag_only = FALSE, ...) {
  d <- NCOL(x1)
  x1 <- prepare_inputmat(object, x1)
  x2 <- prepare_inputmat(object, x2)
  K <- cf_nn_c(x1, x2, object$sigma0, object$sigma, object$magn, diag_only = diag_only)
  if (diag_only) {
    K <- as.vector(K)
  }
  return(K)
}

eval_cf.cf_periodic <- function(object, x1, x2, diag_only = FALSE, ...) {
  x1 <- prepare_inputmat(object, x1)
  x2 <- prepare_inputmat(object, x2)
  period <- object$period
  x1_transf <- cbind(sin(2 * pi / period * x1), cos(2 * pi / period * x1))
  x2_transf <- cbind(sin(2 * pi / period * x2), cos(2 * pi / period * x2))
  K <- eval_cf(object$base, x1_transf, x2_transf, diag_only = diag_only)
  if (diag_only) {
    K <- as.vector(K)
  }
  return(K)
}

eval_cf.cf_prod <- function(object, x1, x2, diag_only = FALSE, ...) {
  K <- 1
  for (k in seq_along(object$cfs)) {
    K <- K * eval_cf(object$cfs[[k]], x1, x2, diag_only = diag_only, ...)
  }
  K
}



# lpdf_prior functions

lpdf_prior.list <- function(object, ...) {
  lp <- 0
  for (k in seq_along(object)) {
    lp <- lp + lpdf_prior(object[[k]])
  }
  lp
}

lpdf_prior.cf <- function(object, ...) {
  param <- get_param(object)
  param_names <- rm_obj_name(object, names(param))
  param_names <- filter_fixed(object, param_names)
  lp <- 0
  for (j in seq_along(param_names)) {
    lp <- lp + lpdf_prior(object$priors[[param_names[j]]], unname(param[j]))
  }
  lp
}

lpdf_prior.cf_prod <- function(object, ...) {
  lpdf_prior(object$cfs)
}

lpdf_prior.cf_periodic <- function(object, ...) {
  lpdf_prior(object$priors$period) + lpdf_prior(object$base)
}



# rf_featmap functions

rf_featmap.list <- function(object, num_feat, ...) {
  fmaps <- list()
  for (k in seq_along(object)) {
    fmaps[[k]] <- rf_featmap(object[[k]], num_feat[k], ...)
  }

  featuremap <- function(x, cfind = NULL) {
    if (is.null(cfind)) {
      cfind <- seq_along(fmaps)
    }
    z <- c()
    for (k in cfind) {
      z <- cbind(z, fmaps[[k]](x))
    }
    return(z)
  }
  return(featuremap)
}

rf_featmap.cf_const <- function(object, ...) {
  featuremap <- function(x) {
    n <- NROW(x)
    object$magn * rep(1, n)
  }
  return(featuremap)
}

rf_featmap.cf_lin <- function(object, ...) {
  # for linear kernel, the linearization feature mapping is simply the identity
  # (with the features scaled by the magnitude)
  featuremap <- function(x) {
    x <- prepare_inputmat(object, x)
    object$magn * x
  }
  return(featuremap)
}

rf_featmap.cf_sexp <- function(object, num_feat, num_inputs, seed = NULL, ...) {
  #
  # spectral density of sexp kernel is given by:
  #     C*N(0,s^2), where
  # s = 1/(2*pi*lscale) and C = (2*pi)^((d-1)/2) * lscale^(d-1)
  #

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  if (is.null(object$vars)) {
    object$vars <- c(1:num_inputs)
  } else {
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  }
  if (num_feat %% 2 == 1) {
    stop("number of features must be an even number.")
  }

  m <- num_feat / 2

  # this is the tricky part; the equations commented out should be correct, but the
  # length-scale and magnitude do not match to the full GP. the simpler equations
  # just seem to work correctly instead..
  scale <- 1 / object$lscale # 1/(2*pi*object$lscale) # scale of the spectral density
  C <- 1 # (2*pi)^((num_inputs-1)/2) * object$lscale^(num_inputs-1) # normalization constant
  w <- matrix(stats::rnorm(m * num_inputs), nrow = num_inputs, ncol = m) * scale # frequences

  featuremap <- function(x) {
    x <- prepare_inputmat(object, x)
    h <- x %*% w
    object$magn * sqrt(C / m) * cbind(cos(h), sin(h))
  }
  return(featuremap)
}

rf_featmap.cf_matern32 <- function(object, num_feat, num_inputs, seed = NULL, ...) {
  # TODO: implement this
  stop("Random Fourier features for Matern kernels not implemented yet.")
}

rf_featmap.cf_matern52 <- function(object, num_feat, num_inputs, seed = NULL, ...) {
  # TODO: implement this
  stop("Random Fourier features for Matern kernels not implemented yet.")
}

rf_featmap.cf_nn <- function(object, num_feat, num_inputs, seed = NULL, ...) {
  #
  # neural network kernel does not have a spectral density (because it's non-stationary),
  # but we can draw the random features by drawing the hidden layer weights from the prior,
  # and then using the probit activations as the features
  #

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  if (is.null(object$vars)) {
    object$vars <- c(1:num_inputs)
  } else {
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  }

  # draw the hidden layer weights randomly
  m <- num_feat
  w <- matrix(stats::rnorm(m * num_inputs), nrow = num_inputs, ncol = m) * object$sigma
  w0 <- matrix(stats::rnorm(m), nrow = 1, ncol = m) * object$sigma0 * object$sigma
  w_aug <- rbind(w0, w)
  erf <- function(t) stats::pnorm(t, sd = sqrt(0.5)) - stats::pnorm(-t, sd = sqrt(0.5))

  featuremap <- function(x) {
    x <- prepare_inputmat(object, x)
    x_aug <- cbind(1, x)
    h <- erf(x_aug %*% w_aug) # hidden layer activations
    object$magn / sqrt(m) * h
  }
  return(featuremap)
}

rf_featmap.cf_periodic <- function(object, num_feat, num_inputs, seed = NULL, ...) {
  if (is.null(object$vars)) {
    object$vars <- c(1:num_inputs)
  } else {
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  }

  featuremap_base <- rf_featmap(object$base, num_feat, num_inputs = 2 * num_inputs, seed = seed)
  featuremap <- function(x) {
    x <- prepare_inputmat(object, x)
    x_transf <- cbind(sin(2 * pi / object$period * x), cos(2 * pi / object$period * x))
    featuremap_base(x_transf)
  }
  return(featuremap)
}

rf_featmap.cf_prod <- function(object, num_feat, num_inputs, seed = NULL, ...) {
  cf_types <- sapply(object$cfs, class)
  if ("cf_lin" %in% cf_types || "cf_const" %in% cf_types) {
    stop("Random features for product kernel containing constant or linear kernel not implemented yet.")
  }
  fmaps <- list()
  for (k in seq_along(object$cfs)) {
    fmaps[[k]] <- rf_featmap(object$cfs[[k]], num_feat, num_inputs, seed, ...)
  }
  featuremap <- function(x) {
    # the random features are obtained by taking the product of the random features
    # of the kernels in the product
    z <- 1
    for (k in seq_along(object$cfs)) {
      z <- z * fmaps[[k]](x)
    }
    z
  }
  return(featuremap)
}







rbf_featmap.list <- function(object, num_feat, ...) {
  fmaps <- list()
  for (k in seq_along(object)) {
    fmaps[[k]] <- rbf_featmap(object[[k]], num_feat[k], ...)
  }

  featuremap <- function(x, cfind = NULL) {
    if (is.null(cfind)) {
      cfind <- seq_along(fmaps)
    }
    z <- c()
    for (k in cfind) {
      z <- cbind(z, fmaps[[k]](x))
    }
    return(z)
  }
  return(featuremap)
}

rbf_featmap.cf_const <- function(object, ...) {
  featuremap <- function(x) {
    n <- NROW(x)
    object$magn * rep(1, n)
  }
  return(featuremap)
}

rbf_featmap.cf_lin <- function(object, ...) {
  # for linear kernel, the linearization feature mapping is simply the identity
  # (with the features scaled by the magnitude)
  featuremap <- function(x) {
    x <- prepare_inputmat(object, x)
    object$magn * x
  }
  return(featuremap)
}

rbf_featmap.cf_sexp <- function(object, num_feat, num_inputs, x = NULL, seed = NULL, ...) {
  if (is.null(x)) {
    stop("Cannot get RBF featuremap if inputs are not given.")
  }

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  if (is.null(object$vars)) {
    object$vars <- c(1:num_inputs)
  } else {
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  }

  # center locations
  x <- as.matrix(x)
  x <- x[, object$vars, drop = FALSE]
  n <- NROW(x)
  ind <- sample(n, num_feat)
  centers <- x[ind, , drop = FALSE]

  featuremap <- function(x) {
    x <- prepare_inputmat(object, x)
    z <- matrix(nrow = NROW(x), ncol = num_feat)
    for (j in 1:num_feat) {
      z[, j] <- exp(-colSums((t(x) - centers[j, ])^2) / object$lscale^2)
    }
    object$magn * sqrt(1 / num_feat) * z
  }
  return(featuremap)
}


rbf_featmap.cf_periodic <- function(object, num_feat, num_inputs, x = NULL, seed = NULL, ...) {
  if (is.null(object$vars)) {
    object$vars <- c(1:num_inputs)
  } else {
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  }

  x_transf <- cbind(sin(2 * pi / object$period * x), cos(2 * pi / object$period * x))
  featuremap_base <- rbf_featmap(object$base, num_feat, num_inputs = 2 * num_inputs, x = x_transf, seed = seed)

  featuremap <- function(x) {
    x <- prepare_inputmat(object, x)
    x_transf <- cbind(sin(2 * pi / object$period * x), cos(2 * pi / object$period * x))
    featuremap_base(x_transf)
  }
  return(featuremap)
}


rbf_featmap.cf_prod <- function(object, num_feat, num_inputs, seed = NULL, ...) {
  cf_types <- sapply(object$cfs, class)
  if ("cf_lin" %in% cf_types || "cf_const" %in% cf_types) {
    stop("RBF features for product kernel containing constant or linear kernel not implemented yet.")
  }
  fmaps <- list()
  for (k in seq_along(object$cfs)) {
    fmaps[[k]] <- rbf_featmap(object$cfs[[k]], num_feat, num_inputs, seed, ...)
  }
  featuremap <- function(x) {
    # the random features are obtained by taking the product of the random features
    # of the kernels in the product
    z <- 1
    for (k in seq_along(object$cfs)) {
      z <- z * fmaps[[k]](x)
    }
    z
  }
  return(featuremap)
}




# learn_scales functions

learn_scales.list <- function(object, x, ...) {
  for (k in seq_along(object)) {
    object[[k]] <- learn_scales(object[[k]], x, ...)
  }
  object
}

learn_scales.gp <- function(object, x, ...) {
  object$cfs <- learn_scales(object$cfs, x, ...)
  object
}

learn_scales.cf <- function(object, x, ...) {
  if (is.null(object$vars)) {
    x <- as.matrix(x)
  } else {
    x <- as.matrix(x)[, object$vars, drop = FALSE]
  }
  object$means <- colMeans(x)
  object$scales <- apply(x, 2, stats::sd)
  object
}

learn_scales.cf_prod <- function(object, x, ...) {
  object$cfs <- learn_scales(object$cfs, x, ...)
  object
}

learn_scales.cf_periodic <- function(object, x, ...) {
  # periodic cf does not implement normalization on purpose
  object
}
