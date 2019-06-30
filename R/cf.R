
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
gpcf_sexp <- function(ind=NULL, lscale=0.5, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_sexp'
  cf
}

#' @rdname cf
#' @export
gpcf_matern32 <- function(ind=NULL, lscale=0.5, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_matern32'
  cf
}

#' @rdname cf
#' @export
gpcf_matern52 <- function(ind=NULL, lscale=0.5, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_matern52'
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


prepare_inputmat <- function(x, ind=NULL) {
  if (is.null(ind))
    return(as.matrix(x))
  else
    return(as.matrix(x)[,ind,drop=F])
}


