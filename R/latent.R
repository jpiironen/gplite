

get_latent_method <- function(name, ...) {
  if (name == 'laplace') {
    method <- latent_laplace(...)
  } else if (name == 'ep') {
    method <- latent_ep(...)
  } else {
    stop('Unknown approximation: ', name)
  }
  method
}


latent_laplace <- function(seed=NULL,  ...) {
  method <- list(name='laplace')
  class(method) <- c('latent_laplace', 'latent')
  method
}

latent_ep <- function(seed=NULL, ep_damping=NULL, ep_quad_order=NULL, ...) {
  method <- list(name='ep', damping=ep_damping, quad_order=ep_quad_order)
  class(method) <- c('latent_ep', 'latent')
  method
}