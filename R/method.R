

# implementations for the methods/types of GPs



#' Initialize method or type of the model
#'
#' Functions for initializing the method or type of the model, which can then be passed to
#' \code{\link{gp_init}}.
#' The supported methods are:
#' \describe{
#'  \item{\code{method_full}}{ Full GP, so full exact covariance function is used, meaning
#' that the inference will be for the \code{n} latent
#' function values (fitting time scales cubicly in \code{n}).}
#'  \item{\code{method_fitc}}{ Fully independent training (and test) conditional, 
#'  or FITC, approximation (see Quiñonero-Candela and Rasmussen, 2005; 
#'  Snelson and Ghahramani, 2006). 
#'  The fitting time scales \code{O(n*m^2)}, where n is the number of data points and 
#'  m the number of inducing points \code{num_inducing}.
#'  The inducing point locations are chosen using the k-means algorithm.
#'  }
#'  \item{\code{method_rf}}{ Random features, that is, linearized GP.
#'  Uses random features (or basis functions) for approximating the covariance function,
#'  which means the inference
#' time scales cubicly in the number of approximating basis functions \code{num_basis}.
#' For stationary covariance functions random Fourier features (Rahimi and Recht, 2007)
#' is used, and for non-stationary kernels using case specific method when possible 
#' (for example, drawing the hidden layer parameters randomly for \code{cf_nn}). For
#' \code{cf_const} and \code{cf_lin} this means using standard linear model, and the
#' inference is performed on the weight space (not in the function space). Thus if
#' the model is linear (only \code{cf_const} and \code{cf_lin} are used), this will give
#' a potentially huge speed-up if the number of features is considerably smaller than
#' the number of data points.}
#' }
#' 
#' 
#'
#' @name method
#'
#' @param num_inducing Number of inducing points for the approximation.
#' @param num_basis Number of basis functions for the approximation.
#' @param bin_inducing TODO: fill this in.
#' @param seed Random seed for reproducible results.
#' 
#'
#' @return The method object.
#' 
#' 
#' @section References:
#' 
#' Rahimi, A. and Recht, B. (2008). Random features for large-scale kernel machines. 
#' In Advances in Neural Information Processing Systems 20.
#' 
#' Quiñonero-Candela, J. and Rasmussen, C. E (2005). A unifying view of sparse approximate 
#' Gaussian process regression. Journal of Machine Learning Research 6:1939-1959.
#' 
#' Snelson, E. and Ghahramani, Z. (2006). Sparse Gaussian processes using pseudo-inputs. 
#' In Advances in Neural Information Processing Systems 18.
#' 
#' @examples
#' \donttest{
#' 
#' # Basic usage 
#' gp <- gp_init(
#'   cfs=cf_sexp(), 
#'   lik=lik_gaussian(), 
#'   method=method_fitc(num_inducing=100)
#' )
#' 
#' }
#'
NULL

# constructors

#' @rdname method
#' @export
method_full <- function() {
  method <- list(name='full')
  class(method) <- c('method_full', 'method')
  method
}

#' @rdname method
#' @export
method_fitc <- function(num_inducing=100, bin_inducing=NULL, seed=12345) {
  method <- list(name='fitc', seed=seed, num_inducing=num_inducing, bin_inducing=bin_inducing)
  class(method) <- c('method_fitc', 'method')
  method
}

#' @rdname method
#' @export
method_rf <- function(num_basis=400, seed=12345) {
  method <- list(name='rf', seed=seed, num_basis=num_basis)
  class(method) <- c('method_rf', 'method')
  method
}

method_rbf <- function(num_basis=400, seed=12345) {
  method <- list(name='rbf', seed=seed, num_basis=num_basis)
  class(method) <- c('method_rbf', 'method')
  method
}

