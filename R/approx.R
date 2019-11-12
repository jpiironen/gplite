

get_approx <- function(name) {
  if (name == 'full') {
    approx <- approx_full()
  } else if (name == 'rf') {
    approx <- approx_rf()
  } else if (name == 'fitc') {
    approx <- approx_fitc()
  } else {
    stop('Unknown approximation: ', name)
  }
  approx
}

approx_full <- function() {
  approx <- list(name='full')
  class(approx) <- c('approx_full', 'approx')
  approx
}

approx_fitc <- function() {
  approx <- list(name='fitc')
  class(approx) <- c('approx_fitc', 'approx')
  approx
}

approx_rf <- function() {
  approx <- list(name='rf')
  class(approx) <- c('approx_rf', 'approx')
  approx
}

