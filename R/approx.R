

get_approx_method <- function(name, ...) {
  if (name == 'laplace') {
    method <- approx_laplace(...)
  } else if (name == 'ep') {
    method <- approx_ep(...)
  } else {
    stop('Unknown approximation: ', name)
  }
  method
}


approx_laplace <- function(seed=NULL,  ...) {
  method <- list(name='laplace')
  class(method) <- c('approx_laplace', 'approx')
  method
}

approx_ep <- function(seed=NULL, ep_damping=NULL, ep_quad_order=NULL, ...) {
  method <- list(name='ep', damping=ep_damping, quad_order=ep_quad_order)
  class(method) <- c('approx_ep', 'approx')
  method
}