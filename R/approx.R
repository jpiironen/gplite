

get_approx <- function(name, ...) {
  if (name == 'full') {
    approx <- approx_full(...)
  } else if (name == 'rf') {
    approx <- approx_rf(...)
  } else if (name == 'rbf') {
    approx <- approx_rbf(...)
  } else if (name == 'fitc') {
    approx <- approx_fitc(...)
  } else {
    stop('Unknown approximation: ', name)
  }
  approx
}

approx_full <- function(...) {
  approx <- list(name='full')
  class(approx) <- c('approx_full', 'approx')
  approx
}

approx_fitc <- function(seed=NULL, num_inducing=NULL, ...) {
  approx <- list(name='fitc', seed=seed, num_inducing=num_inducing)
  class(approx) <- c('approx_fitc', 'approx')
  approx
}

approx_rf <- function(seed=NULL, num_basis=NULL, ...) {
  approx <- list(name='rf', seed=seed, num_basis=num_basis)
  class(approx) <- c('approx_rf', 'approx')
  approx
}

approx_rbf <- function(seed=NULL, num_basis=NULL, ...) {
  approx <- list(name='rbf', seed=seed, num_basis=num_basis)
  class(approx) <- c('approx_rbf', 'approx')
  approx
}

