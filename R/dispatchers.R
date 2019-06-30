
#' Get or set GP model parameters
#'
#' \code{get_param} returns the current hyperparameters of the GP model in a vector. \code{set_param} can be used to set the parameters.
#' 
#' @name param
#' 
#' @param object The model object.
#' @param param The parameters to be set. Call first \code{get_param} to see the order in which the parameters should be given for a particular model. Notice that all positive parameters should be given in a log-scale.
#' @param ... Ignored currently.
#' 
#' @return \code{get_param} returns the current hyperparameters and \code{set_param} the GP model structure with the new parameter values.
#' 
#'
#' @examples
#' \donttest{
#' # Basic usage 
#' gp <- gp_init(cf=gpcf_sexp(), lik=lik_gaussian())
#' print(get_param(gp)) # print out to see the parameter ordering
#' gp <- set_param(gp, log(c(0.1,0.8,0.3))) # set some new values
#' print(get_param(gp)) # check the result
#' 
#' }
#'
NULL

#' @rdname param
#' @export
get_param <- function (object, ...) {
  UseMethod("get_param", object)
}

#' @rdname param
#' @export
set_param <- function (object, param, ...) {
  UseMethod("set_param", object)
}

get_stanmodel <- function(object, ...) {
  UseMethod("get_stanmodel", object)
}

get_standata <- function(object, ...) {
  UseMethod("get_standata", object)
}

get_response <- function(object, ...) {
  UseMethod("get_response", object)
}

eval_cf <- function (object, ...) {
  UseMethod("eval_cf", object)
}

rff_featmap <- function (object, ...) {
  UseMethod("rff_featmap", object)
}

is_fitted <- function(object, type, ...) {
  UseMethod("is_fitted", object)
}



