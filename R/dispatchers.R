
#' Get or set GP model parameters
#'
#' \code{get_param} returns the current hyperparameters of the GP model in a vector.
#' \code{set_param} can be used to set the parameters. Note that these functions
#' are intended mainly for internal usage, and there is typically
#' no need to use these functions directly but instead create a new GP model using
#' \code{gp_init}.
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
#'
#' # Set up some model
#' gp <- gp_init(cf = cf_sexp(), lik = lik_gaussian())
#'
#' # print out to see the parameter ordering
#' param <- get_param(gp)
#' print(param)
#'
#' # set some new values
#' param_new <- log(c(0.1, 0.8, 0.3))
#' names(param_new) <- names(param)
#' gp <- set_param(gp, param_new)
#'
#' # check the result
#' print(get_param(gp))
#' 
#'
NULL

#' @rdname param
#' @export
get_param <- function(object, ...) {
  UseMethod("get_param", object)
}

#' @rdname param
#' @export
set_param <- function(object, param, ...) {
  UseMethod("set_param", object)
}

get_name <- function(object, ...) {
  UseMethod("get_name")
}

get_param_names <- function(object, ...) {
  UseMethod("get_param_names", object)
}

get_featuremap <- function(object, ...) {
  UseMethod("get_featuremap", object)
}

get_pseudodata_la <- function(object, ...) {
  UseMethod("get_pseudodata_la", object)
}

get_pseudodata_ep <- function(object, ...) {
  UseMethod("get_pseudodata_ep", object)
}

get_tilted_moments <- function(object, ...) {
  UseMethod("get_tilted_moments", object)
}

get_loglik <- function(object, ...) {
  UseMethod("get_loglik", object)
}

get_loglik_d <- function(object, ...) {
  UseMethod("get_loglik_d", object)
}

get_loglik_d2 <- function(object, ...) {
  UseMethod("get_loglik_d2", object)
}

get_response <- function(object, ...) {
  UseMethod("get_response", object)
}

generate_target <- function(object, ...) {
  UseMethod("generate_target", object)
}

eval_cf <- function(object, ...) {
  UseMethod("eval_cf", object)
}

rf_featmap <- function(object, ...) {
  UseMethod("rf_featmap", object)
}

rbf_featmap <- function(object, ...) {
  UseMethod("rbf_featmap", object)
}

is_fitted <- function(object, type, ...) {
  UseMethod("is_fitted", object)
}

lpdf_prior <- function(object, ...) {
  UseMethod("lpdf_prior", object)
}

learn_scales <- function(object, ...) {
  UseMethod("learn_scales", object)
}

required_extra_args <- function(object, ...) {
  UseMethod("required_extra_args", object)
}
