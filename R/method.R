

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
#' @param inducing Inducing points to use. If not given, then \code{num_inducing}
#' points will be placed in the input space using a clustering algorithm.
#' @param num_inducing Number of inducing points for the approximation. Will be ignored
#' if the inducing points are given by the user.
#' @param bin_along Either an index or a name of the input variable along which to bin the
#' values before placing the inducing inputs. For example, if \code{bin_along=3}, then the
#' input data is divided into \code{bin_count} bins along 3rd input variable, and each bin
#' will have the same number inducing points (or as close as possible). This can sometimes
#' be useful to ensure that inducing points are spaced evenly with respect to some particular
#' variable, for example time in spatio-temporal models.
#' @param bin_count The number of bins to use if \code{bin_along} given. Has effect only if
#' \code{bin_along} is given.
#' @param num_basis Number of basis functions for the approximation.
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
#'
#' 
#' #' # Generate some toy data
#' # NOTE: this is so small dataset that in reality there would be no point
#' # use sparse approximation here; we use this small dataset only to make this
#' # example run fast
#' set.seed(1242)
#' n <- 50
#' x <- matrix(rnorm(n * 3), nrow = n)
#' f <- sin(x[, 1]) + 0.5 * x[, 2]^2 + x[, 3]
#' y <- f + 0.5 * rnorm(n)
#' x <- data.frame(x1 = x[, 1], x2 = x[, 2], x3 = x[, 3])
#'
#' # Full exact GP with Gaussian likelihood
#' gp <- gp_init(cf_sexp())
#' gp <- gp_optim(gp, x, y)
#'
#' # Approximate solution using random features (here we use a very small 
#' # number of random features only to make this example run fast)
#' gp <- gp_init(cf_sexp(), method = method_rf(num_basis = 30))
#' gp <- gp_optim(gp, x, y)
#'
#' # Approximate solution using FITC (here we use a very small 
#' # number of incuding points only to make this example run fast)
#' gp <- gp_init(cf_sexp(), method = method_fitc(num_inducing = 10))
#' gp <- gp_optim(gp, x, y)
#'
NULL

# constructors

#' @rdname method
#' @export
method_full <- function() {
  method <- list(name = "full")
  class(method) <- c("method_full", "method")
  method
}

#' @rdname method
#' @export
method_fitc <- function(inducing = NULL, num_inducing = 100,
                        bin_along = NULL, bin_count = 10, seed = 12345) {
  if (!is.null(inducing)) {
    num_inducing <- NROW(inducing)
  }
  bin_count <- min(bin_count, num_inducing)
  method <- list(
    name = "fitc", seed = seed, inducing = inducing,
    num_inducing = num_inducing, bin_along = bin_along, bin_count = bin_count
  )
  class(method) <- c("method_fitc", "method")
  method
}

#' @rdname method
#' @export
method_rf <- function(num_basis = 400, seed = 12345) {
  method <- list(name = "rf", seed = seed, num_basis = num_basis)
  class(method) <- c("method_rf", "method")
  method
}

method_rbf <- function(num_basis = 400, seed = 12345) {
  method <- list(name = "rbf", seed = seed, num_basis = num_basis)
  class(method) <- c("method_rbf", "method")
  method
}


#' @export
print.method <- function(x, quiet = FALSE, ...) {
  object <- x
  str <- class(object)[1]
  if (!quiet) {
    cat(str)
  }
  invisible(str)
}
