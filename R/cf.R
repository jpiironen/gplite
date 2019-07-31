
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
#'
#' @name cf
#'
#' @param vars Indices of the inputs which are taken into account when calculating this
#'  covariance. If the input matrix has named columns, can also be a vector of column names.
#'  Default is all the inputs.
#' @param lscale Initial value for the length-scale hyperparameter.
#' @param magn Initial value for the magnitude hyperparameter (depicts the magnitude of 
#' the variation captured by the given covariance function).
#' @param sigma0 Prior std for the bias in the neural network covariance function.
#' @param sigma Prior std for the weights in the hidden layers of the neural network 
#' covariance function. 
#' @param period Period length for the periodic covariance function.
#' @param cf_base Base covariance function that is used to model the variability within each period 
#' in periodic covariance function.
#' @param ... Meaning depends on context. For \code{cf_prod} pass in the covariance functions in the product. 
#'
#' @return The covariance function object.
#' 
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#' 
#' @examples
#' \donttest{
#' # Basic usage (single covariance function)
#' cf <- cf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x ,y, trials)
#' 
#' 
#' # More than one covariance function, and the function is additive 
#' # in the first input, but has an interaction between inputs 2 and 3
#' # (the functional form is f(x) = f(x_1) + f(x_2,x_3))
#' cfs <- list(cf_const(), cf_sexp(1), cf_sexp(c(2,3)))
#' lik <- lik_gaussian()
#' gp <- gp_init(cfs, lik)
#' gp <- gp_optim(gp, x, y)
#' 
#' }
#'
NULL


# constructors

#' @rdname cf
#' @export
cf_const <- function(magn=1.0, prior_magn=prior_logunif()) {
  cf <- list()
  cf$magn <- magn
  cf$priors <- list(magn=prior_magn)
  class(cf) <- c('cf_const', 'cf')
  cf
}

#' @rdname cf
#' @export
cf_lin <- function(vars=NULL, magn=1.0, prior_magn=prior_logunif()) {
  cf <- list()
  cf$vars <- vars
  cf$magn <- magn
  cf$priors <- list(magn=prior_magn)
  class(cf) <- c('cf_lin', 'cf')
  cf
}

#' @rdname cf
#' @export
cf_sexp <- function(vars=NULL, lscale=0.3, magn=1.0,
                    prior_lscale=prior_logunif(), prior_magn=prior_logunif()) {
  cf <- list()
  cf$vars <- vars
  cf$lscale <- lscale
  cf$magn <- magn
  cf$priors <- list(lscale=prior_lscale, magn=prior_magn)
  class(cf) <- c('cf_sexp', 'cf')
  cf
}

#' @rdname cf
#' @export
cf_matern32 <- function(vars=NULL, lscale=0.3, magn=1.0,
                        prior_lscale=prior_logunif(), prior_magn=prior_logunif()) {
  cf <- list()
  cf$vars <- vars
  cf$lscale <- lscale
  cf$magn <- magn
  cf$priors <- list(lscale=prior_lscale, magn=prior_magn)
  class(cf) <- c('cf_matern32', 'cf')
  cf
}

#' @rdname cf
#' @export
cf_matern52 <- function(vars=NULL, lscale=0.3, magn=1.0,
                        prior_lscale=prior_logunif(), prior_magn=prior_logunif()) {
  cf <- list()
  cf$vars <- vars
  cf$lscale <- lscale
  cf$magn <- magn
  cf$priors <- list(lscale=prior_lscale, magn=prior_magn)
  class(cf) <- c('cf_matern52', 'cf')
  cf
}

#' @rdname cf
#' @export
cf_nn <- function(vars=NULL, sigma0=1.0, sigma=1.0, magn=1.0,
                  prior_sigma0=prior_logunif(), prior_sigma=prior_logunif(), prior_magn=prior_logunif()) {
  cf <- list()
  cf$vars <- vars
  cf$sigma0 <- sigma0
  cf$sigma <- sigma
  cf$magn <- magn
  cf$priors <- list(sigma0=prior_sigma0, sigma=prior_sigma, magn=prior_magn)
  class(cf) <- c('cf_nn', 'cf')
  cf
}

#' @rdname cf
#' @export 
cf_periodic <- function(vars=NULL, period=1, cf_base=cf_sexp(), prior_period=prior_logunif()) {
  cf <- list()
  cf$vars <- vars
  cf$period <- period
  cf$base <- cf_base
  cf$priors <- list(period=prior_period)
  class(cf) <- c('cf_periodic', 'cf')
  cf
}

#' @rdname cf
#' @export
cf_prod <- function(...) {
  cf <- list()
  cf$cfs <- list(...)
  class(cf) <- c('cf_prod', 'cf')
  cf
}




get_name.cf <- function(object, ...) {
  class(object)[1]
}



# get_param_names functions

get_param_names.cf_const <- function(object) {
  c('magn')
}

get_param_names.cf_lin <- function(object) {
  c('magn')
}

get_param_names.cf_sexp <- function(object) {
  c('lscale', 'magn')
}

get_param_names.cf_matern32 <- function(object) {
  c('lscale', 'magn')
}

get_param_names.cf_matern52 <- function(object) {
  c('lscale', 'magn')
}

get_param_names.cf_nn <- function(object) {
  c('sigma0', 'sigma', 'magn')
}

get_param_names.cf_periodic <- function(object) {
  c('period')
}



# get_param functions

get_param.list <- function(object, ...) {
  param <- c()
  for (k in seq_along(object))
    param <- c(param, get_param(object[[k]]))
  param
}

get_param.cf <- function(object, ...) {
  param_names <- filter_fixed(object, get_param_names(object))
  if (length(param_names) == 0)
    return(NULL)
  param <- unlist(object[param_names])
  names(param) <- add_obj_name(object, names(param))
  param <- log(param)
  param
}

get_param.cf_periodic <- function(object, ...) {
  param <- get_param(object$base)
  if (!is_fixed(object, 'period')) {
    param <- c(log(object$period), param)
    names(param)[1] <- 'cf_periodic.period'
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
    if (np > 0)
      object[[k]] <- set_param(object[[k]], param[j:(j+np-1)])
    j <- j + np
  }
  object
}

set_param.cf <- function(object, param, ...) {
  param_names <- filter_fixed(object, names(param))
  param_names <- rm_obj_name(object, param_names)
  for (j in seq_along(param_names))
    object[[param_names[j]]] <- unname(exp(param[j]))
  object
}

set_param.cf_periodic <- function(object, param, ...) {
  fixed_period <- is_fixed(object, 'period')
  if (!fixed_period)
    object$period <- exp(param[1])
  object$base <- set_param(object$base, tail(param, length(param)-fixed_period))
  object
}

set_param.cf_prod <- function(object, param, ...) {
  object$cfs <- set_param(object$cfs, param)
  object
}



# eval_cf functions

eval_cf.list <- function(object, x1, x2, cfind=NULL, ...) {
  if (is.null(cfind))
    cfind <- seq_along(object)
  K <- 0
  for (k in cfind) {
    K <- K + eval_cf(object[[k]], x1, x2, ...)
  }
  K
}

eval_cf.cf_const <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1)
  x2 <- prepare_inputmat(x2)
  n1 <- NROW(x1)
  n2 <- NROW(x2)
  K <- matrix(object$magn^2, nrow=n1,ncol=n2)
  return(K)
}

eval_cf.cf_lin <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$vars)
  x2 <- prepare_inputmat(x2, object$vars)
  K <- object$magn^2* x1 %*% t(x2)
  return(K)
}

eval_cf.cf_sexp <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$vars)
  x2 <- prepare_inputmat(x2, object$vars)
  K <- cf_sexp_c(x1, x2, object$lscale, object$magn)
  return(K)
}

eval_cf.cf_matern32 <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$vars)
  x2 <- prepare_inputmat(x2, object$vars)
  K <- cf_matern32_c(x1, x2, object$lscale, object$magn)
  return(K)
}

eval_cf.cf_matern52 <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$vars)
  x2 <- prepare_inputmat(x2, object$vars)
  K <- cf_matern52_c(x1, x2, object$lscale, object$magn)
  return(K)
}

eval_cf.cf_nn <- function(object, x1, x2, ...) {
  d <- NCOL(x1)
  x1 <- cbind(1, prepare_inputmat(x1, object$vars))
  x2 <- cbind(1, prepare_inputmat(x2, object$vars))
  K <- cf_nn_c(x1, x2, object$sigma0, object$sigma, object$magn)
  return(K)
}

eval_cf.cf_periodic <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$vars)
  x2 <- prepare_inputmat(x2, object$vars)
  period <- object$period
  x1_transf <- cbind(sin(2*pi/period*x1), cos(2*pi/period*x1))
  x2_transf <- cbind(sin(2*pi/period*x2), cos(2*pi/period*x2))
  K <- eval_cf(object$base, x1_transf, x2_transf)
  return(K)
}

eval_cf.cf_prod <- function(object, x1, x2, ...) {
  K <- 1
  for (k in seq_along(object$cfs))
    K <- K*eval_cf(object$cfs[[k]], x1, x2, ...)
  K
}





# rf_featmap functions

rf_featmap.list <- function(object, num_feat, ...) {
  fmaps <- list()
  for (k in seq_along(object))
    fmaps[[k]] <- rf_featmap(object[[k]], num_feat[k], ...)
  
  featuremap <- function(x, cfind=NULL) {
    if (is.null(cfind))
      cfind <- seq_along(fmaps)
    z <- c()
    for (k in cfind)
      z <- cbind(z,fmaps[[k]](x))
    return(z)
  }
  return(featuremap)
}

rf_featmap.cf_const <- function(object, ...) {
  featuremap <- function(x) {
    n <- NROW(x)
    object$magn*rep(1,n)
  }
  return(featuremap)
}

rf_featmap.cf_lin <- function(object, ...) {
  # for linear kernel, the linearization feature mapping is simply the identity
  # (with the features scaled down by the magnitude)
  featuremap <- function(x) {
    object$magn*x
  }
  return(featuremap)
}

rf_featmap.cf_sexp <- function(object, num_feat, num_inputs, seed=NULL, ...) {
  #
  # spectral density of sexp kernel is given by:
  #     C*N(0,s^2), where
  # s = 1/(2*pi*lscale) and C = (2*pi)^((d-1)/2) * lscale^(d-1)
  #
  
  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(seed)
  
  if (is.null(object$vars))
    object$vars <- c(1:num_inputs)
  else
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  if (num_feat %% 2 == 1)
    stop('number of features must be an even number.')
  
  m <- num_feat/2
  scale <- 1/(2*pi*object$lscale) # scale of the spectral density
  w <- matrix(stats::rnorm(m*num_inputs), nrow=num_inputs, ncol=m)*scale # frequences
  C <- (2*pi)^((num_inputs-1)/2) * object$lscale^(num_inputs-1)
  
  featuremap <- function(x) {
    x <- as.matrix(x)
    h <- x[,object$vars,drop=F] %*% w
    object$magn*sqrt(C/m)*cbind(cos(h),sin(h))
  }
  return(featuremap)
}

rf_featmap.cf_matern32 <- function(object, num_feat, num_inputs, seed=NULL, ...) {
  # TODO: implement this
  stop('Random Fourier features for Matern kernels not implemented yet.')
}

rf_featmap.cf_matern52 <- function(object, num_feat, num_inputs, seed=NULL, ...) {
  # TODO: implement this
  stop('Random Fourier features for Matern kernels not implemented yet.')
}

rf_featmap.cf_nn <- function(object, num_feat, num_inputs, seed=NULL, ...) {
  #
  # neural network kernel does not have a spectral density (because it's non-stationary),
  # but we can draw the random features by drawing the hidden layer weights from the prior, 
  # and then using the probit activations as the features
  #
  
  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(seed)
  
  if (is.null(object$vars))
    object$vars <- c(1:num_inputs)
  else
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  
  # draw the hidden layer weights randomly
  m <- num_feat
  w <- matrix(stats::rnorm(m*num_inputs), nrow=num_inputs, ncol=m)*object$sigma
  w0 <- matrix(stats::rnorm(m), nrow=1, ncol=m)*object$sigma0
  w_aug <- rbind(w0, w)
  
  featuremap <- function(x) {
    x_aug <- cbind(1, as.matrix(x)[,object$vars,drop=F])
    h <- stats::pnorm(x_aug %*% w_aug) # hidden layer activations
    object$magn/sqrt(m)*h
  }
  return(featuremap)
}

rf_featmap.cf_periodic <- function(object, num_feat, num_inputs, seed=NULL, ...) {
  
  if (is.null(object$vars))
    object$vars <- c(1:num_inputs)
  else
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$vars)
  
  featuremap_base <- rf_featmap(object$base, num_feat, num_inputs=2*length(object$vars), seed=seed)
  featuremap <- function(x) {
    x_transf <- cbind(sin(2*pi/object$period*x), cos(2*pi/object$period*x))
    featuremap_base(x_transf)
  }
  return(featuremap)
}

rf_featmap.cf_prod <- function(object, num_feat, num_inputs, seed=NULL, ...) {
  cf_types <- sapply(object$cfs, class)
  if ('cf_lin' %in% cf_types || 'cf_const' %in% cf_types)
    stop('Random features for product kernel containing constant or linear kernel not implemented yet.')
  fmaps <- list()
  for (k in seq_along(object$cfs))
    fmaps[[k]] <- rf_featmap(object$cfs[[k]], num_feat, num_inputs, seed, ...)
  featuremap <- function(x) {
    # the random features are obtained by taking the product of the random features 
    # of the kernels in the product
    z <- 1
    for (k in seq_along(object$cfs))
      z <- z*fmaps[[k]](x)
    z
  }
  return(featuremap)
}









