.onAttach <- function(...) {
  ver <- utils::packageVersion("gplite")
  packageStartupMessage("This is gplite version ", ver)
}

add_offset <- function(f, ...) {
  # add offset to the vector of latent values if offset was given
  # in the additional arguments
  args <- list(...)
  offset <- args$offset
  if (is.null(offset)) {
    offset <- 0
  }
  return(f + offset)
}

prepare_inputmat <- function(cf, x) {
  # ensure x is a matrix, and pick the correct columns
  vars <- cf$vars
  if (is.null(vars)) {
    x <- as.matrix(x)
  } else {
    x <- as.matrix(x)[, vars, drop = FALSE]
  }
  if (!is.null(cf$normalize) && cf$normalize) {
    if (!is.null(cf$means) && !is.null(cf$scales)) {
      x <- t((t(x) - cf$means) / cf$scales)
    }
  }
  return(x)
}

is_fixed <- function(object, param_names) {
  # identify fixed hyperparameter(s) of a given object (cf or lik)
  fixed <- sapply(param_names, function(name) {
    "prior_fixed" %in% class(object$priors[[name]])
  })
  return(fixed)
}

filter_fixed <- function(object, param_names) {
  # take only parameter names that are not fixed
  param_names[!is_fixed(object, param_names)]
}

add_obj_name <- function(object, param_names) {
  # attach the object name (cf of lik) to the parameter names
  paste0(get_name(object), ".", param_names)
}

rm_obj_name <- function(object, param_names) {
  # remove the object name from parameter names
  sapply(param_names, function(name) {
    unlist(strsplit(name, ".", fixed = TRUE))[2]
  })
}

check_if_overparametrized_magnitude <- function(cf) {
  if (!("cf_prod" %in% class(cf))) {
    stop(paste(
      "Internal error encountered: checking for overparametrized magnitude for a non-product kernel.",
      "This is a bug, please report to the developers!"
    ))
  }
  # check that there are no multiple free magnitude parameters
  param <- get_param(cf)
  param_names <- names(param)
  num_free_magn_params <- 0
  for (i in seq_along(param_names)) {
    if (grepl("magn", param_names[i], fixed = TRUE)) {
      num_free_magn_params <- num_free_magn_params + 1
    }
  }
  if (num_free_magn_params > 1)
    warning(paste(
      "Detected a product of covariance functions with more than one non-fixed magnitude parameters.",
      "This means the magnitude of the product kernel is unnecessarily overparametrized, which",
      "is likely to lead to to problems in hyperparameter optimization. Consider setting",
      "prior_magn = prior_fixed() for all except one of the kernels in the product."
    ))
}
