
# implementations for the covariance functions

# TODO: add matern kernel



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
gpcf_sexp <- function(ind=NULL, lscale=0.5, magn=1.0) {
  cf <- list()
  cf$ind <- ind
  cf$lscale <- lscale
  cf$magn <- magn
  class(cf) <- 'cf_sexp'
  cf
}




# get_param functions

get_param.cf_const <- function(object, ...) {
  param <- log(object$magn)
  names(param) <- 'cf_const.magn'
  param
}

get_param.cf_sexp <- function(object, ...) {
  param <- log(c(object$lscale, object$magn))
  names(param) <- c('cf_sexp.lscale','cf_sexp.magn')
  param
}




# set_param functions

set_param.cf_const <- function(object, param, ...) {
  object$magn <- exp(param[1])
  object
}

set_param.cf_sexp <- function(object, param, ...) {
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

eval_cf.cf_sexp <- function(object, x1, x2, ...) {
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  n1 <- NROW(x1)
  n2 <- NROW(x2)
  K <- matrix(nrow=n1,ncol=n2)
  lscale <- object$lscale
  magn <- object$magn
  ind <- object$ind
  if (is.null(ind))
    ind <- c(1:NCOL(x1))
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      K[i,j] <- magn^2*exp(-sum((x1[i,]-x2[j,])^2/lscale^2))
    }
  }
  return(K)
}

eval_cf.cf_const <- function(object, x1, x2, ...) {
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  n1 <- NROW(x1)
  n2 <- NROW(x2)
  K <- matrix(object$magn^2, nrow=n1,ncol=n2)
  return(K)
}
