
# implementations for the covariance functions




#' Initialize covariance function
#'
#' Functions for initializing the covariance functions which can then be passed to \code{\link{gp_init}}.
#'
#' @name cf
#'
#' @param ind Indices of the inputs which are taken into account when calculating this covariance. Default is all the inputs.
#' @param lscale Initial value for the length-scale hyperparameter.
#' @param magn Initial value for the magnitude hyperparameter (depicts the magnitude of the variation captured by the given covariance function).
#'
#' @return The covariance function object.
#' 
#' @examples
#' \donttest{
#' # Basic usage (single covariance function)
#' cf <- gpcf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x ,y, trials)
#' 
#' 
#' # More than one covariance function, and the function is additive 
#' # in the first input, but has an interaction between inputs 2 and 3
#' # (the functional form is f(x) = f(x_1) + f(x_2,x_3))
#' cfs <- list(gpcf_const(), gpcf_sexp(1), gpcf_sexp(c(2,3)))
#' lik <- lik_gaussian()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x, y)
#' 
#' }
#'
NULL


# constructors

#' @rdname cf
#' @export
gpcf_const <- function(magn=1.0) {
  cf <- list()
  cf$magn <- magn
  class(cf) <- 'cf_const'
  cf
}

#' @rdname cf
#' @export
gpcf_lin <- function(ind=NULL, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$magn <- magn
  class(cf) <- 'cf_lin'
  cf
}

#' @rdname cf
#' @export
gpcf_sexp <- function(ind=NULL, lscale=0.1, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_sexp'
  cf
}

#' @rdname cf
#' @export
gpcf_matern32 <- function(ind=NULL, lscale=0.1, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_matern32'
  cf
}

#' @rdname cf
#' @export
gpcf_matern52 <- function(ind=NULL, lscale=0.1, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_matern52'
  cf
}

#' @rdname cf
#' @export
gpcf_nn <- function(ind=NULL, sigma0=1.0, sigma=1.0, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$sigma0 <- sigma0
  cf$sigma <- sigma
  cf$magn <- magn
  class(cf) <- 'cf_nn'
  cf
}



# get_param functions

get_param.cf_const <- function(object, ...) {
  param <- log(object$magn)
  names(param) <- 'cf_const.magn'
  param
}

get_param.cf_lin <- function(object, ...) {
  param <- log(object$magn)
  names(param) <- 'cf_lin.magn'
  param
}

get_param.cf_sexp <- function(object, ...) {
  param <- log(c(object$lscale, object$magn))
  names(param) <- c('cf_sexp.lscale','cf_sexp.magn')
  param
}

get_param.cf_matern32 <- function(object, ...) {
  param <- log(c(object$lscale, object$magn))
  names(param) <- c('cf_matern32.lscale','cf_matern32.magn')
  param
}

get_param.cf_matern52 <- function(object, ...) {
  param <- log(c(object$lscale, object$magn))
  names(param) <- c('cf_matern52.lscale','cf_matern52.magn')
  param
}

get_param.cf_nn <- function(object, ...) {
  param <- log(c(object$sigma0, object$sigma, object$magn))
  names(param) <- c('cf_nn.sigma0', 'cf_nn.sigma', 'cf_nn.magn')
  param
}



# set_param functions

set_param.cf_const <- function(object, param, ...) {
  object$magn <- exp(param[1])
  object
}

set_param.cf_lin <- function(object, param, ...) {
  object$magn <- exp(param[1])
  object
}

set_param.cf_sexp <- function(object, param, ...) {
  object$lscale <- exp(param[1])
  object$magn <- exp(param[2])
  object
}

set_param.cf_matern32 <- function(object, param, ...) {
  object$lscale <- exp(param[1])
  object$magn <- exp(param[2])
  object
}

set_param.cf_matern52 <- function(object, param, ...) {
  object$lscale <- exp(param[1])
  object$magn <- exp(param[2])
  object
}

set_param.cf_nn <- function(object, param, ...) {
  object$sigma0 <- exp(param[1])
  object$sigma <- exp(param[2])
  object$magn <- exp(param[3])
  object
}



# eval_cf functions

eval_cf.list <- function(object, x1, x2, ...) {
  K <- 0
  for (k in seq_along(object)) {
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
  x1 <- prepare_inputmat(x1, object$ind)
  x2 <- prepare_inputmat(x2, object$ind)
  K <- object$magn^2* x1 %*% t(x2)
  return(K)
}

eval_cf.cf_sexp <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$ind)
  x2 <- prepare_inputmat(x2, object$ind)
  K <- cf_sexp_c(x1, x2, object$lscale, object$magn)
  return(K)
}

eval_cf.cf_matern32 <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$ind)
  x2 <- prepare_inputmat(x2, object$ind)
  K <- cf_matern32_c(x1, x2, object$lscale, object$magn)
  return(K)
}

eval_cf.cf_matern52 <- function(object, x1, x2, ...) {
  x1 <- prepare_inputmat(x1, object$ind)
  x2 <- prepare_inputmat(x2, object$ind)
  K <- cf_matern52_c(x1, x2, object$lscale, object$magn)
  return(K)
}

eval_cf.cf_nn <- function(object, x1, x2, ...) {
  d <- NCOL(x1)
  x1 <- cbind(1, prepare_inputmat(x1, object$ind))
  x2 <- cbind(1, prepare_inputmat(x2, object$ind))
  K <- cf_nn_c(x1, x2, object$sigma0, object$sigma, object$magn)
  return(K)
}


prepare_inputmat <- function(x, ind=NULL) {
  if (is.null(ind))
    return(as.matrix(x))
  else
    return(as.matrix(x)[,ind,drop=F])
}



# rf_featmap functions

rf_featmap.list <- function(object, ...) {
  fmaps <- list()
  for (k in seq_along(object))
    fmaps[[k]] <- rf_featmap(object[[k]], ...)
  
  featuremap <- function(x) {
    z <- c()
    for (k in seq_along(fmaps))
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

rf_featmap.cf_sexp <- function(object, num_inputs, num_feat, seed=NULL, ...) {
  #
  # spectral density of sexp kernel is given by:
  #     C*N(0,scale^2), where
  # scale = 1/(2*pi*lscale) and C = (2*pi)^(d-1)/2*lscale^(d-1)
  #
  
  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(seed)
  
  if (is.null(object$ind))
    object$ind <- c(1:num_inputs)
  else
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$ind)
  if (num_feat %% 2 == 1)
    stop('number of features must be an even number.')
  
  m <- num_feat/2
  scale <- 1/(2*pi*object$lscale) # scale of the spectral density
  w <- matrix(stats::rnorm(m*num_inputs), nrow=num_inputs, ncol=m)*scale # frequences
  C <- (2*pi)^(num_inputs-1)/2*object$lscale^(num_inputs-1)
  
  featuremap <- function(x) {
    x <- as.matrix(x)
    h <- x[,object$ind,drop=F] %*% w
    object$magn*sqrt(C/m)*cbind(cos(h),sin(h))
  }
  return(featuremap)
}

rf_featmap.cf_matern32 <- function(object, num_inputs, num_feat, seed=NULL, ...) {
  # TODO: implement this
  stop('Random Fourier feature for Matern kernels not implemented yet.')
}

rf_featmap.cf_matern52 <- function(object, num_inputs, num_feat, seed=NULL, ...) {
  # TODO: implement this
  stop('Random Fourier feature for Matern kernels not implemented yet.')
}

rf_featmap.cf_nn <- function(object, num_inputs, num_feat, seed=NULL, ...) {
  #
  # neural network kernel does not have a spectral density (because it's non-stationary),
  # but we can draw the random features by drawing the hidden layer weights from the prior, 
  # and then using the probit activations as the features
  #
  
  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(seed)
  
  if (is.null(object$ind))
    object$ind <- c(1:num_inputs)
  else
    # override the number of inputs, because using only a subset of inputs
    num_inputs <- length(object$ind)
  if (num_feat %% 2 == 1)
    stop('number of features must be an even number.')
  
  # draw the hidden layer weights randomly
  m <- num_feat
  w <- matrix(stats::rnorm(m*num_inputs), nrow=num_inputs, ncol=m)*object$sigma
  w0 <- matrix(stats::rnorm(m), nrow=1, ncol=m)*object$sigma0
  w_aug <- rbind(w0, w)
  
  featuremap <- function(x) {
    x_aug <- cbind(1, as.matrix(x)[,object$ind,drop=F])
    h <- stats::pnorm(x_aug %*% w_aug) # hidden layers activations
    object$magn/sqrt(m)*h
  }
  return(featuremap)
  
}












