

get_param <- function (object, ...) {
  UseMethod("get_param", object)
}

set_param <- function (object, ...) {
  UseMethod("set_param", object)
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