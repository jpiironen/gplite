

get_method <- function(name, ...) {
  if (name == 'full') {
    method <- method_full(...)
  } else if (name == 'rf') {
    method <- method_rf(...)
  } else if (name == 'rbf') {
    method <- method_rbf(...)
  } else if (name == 'fitc') {
    method <- method_fitc(...)
  } else {
    stop('Unknown method: ', name)
  }
  method
}

method_full <- function(...) {
  method <- list(name='full')
  class(method) <- c('method_full', 'method')
  method
}

method_fitc <- function(seed=NULL, num_inducing=NULL, bin_inducing=NULL, ...) {
  method <- list(name='fitc', seed=seed, num_inducing=num_inducing, bin_inducing=bin_inducing)
  class(method) <- c('method_fitc', 'method')
  method
}

method_rf <- function(seed=NULL, num_basis=NULL, ...) {
  method <- list(name='rf', seed=seed, num_basis=num_basis)
  class(method) <- c('method_rf', 'method')
  method
}

method_rbf <- function(seed=NULL, num_basis=NULL, ...) {
  method <- list(name='rbf', seed=seed, num_basis=num_basis)
  class(method) <- c('method_rbf', 'method')
  method
}

