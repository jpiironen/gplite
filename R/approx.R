


# implementations for the methods/types of GPs



#' Approximations to the posterior of the latent values
#'
#' Functions for initializing the approximation for the latent values, which can
#' then be passed to \code{\link{gp_init}}.
#' The supported methods are:
#' \describe{
#'  \item{\code{approx_laplace}}{ Laplace's method, that is, based on local second
#'   order approximation to the log likelihood. For Gaussian likelihood, this means exact inference
#'    (no approximation). }
#'  \item{\code{approx_ep}}{ Expectation propagation, EP. Approximates the likelihood by
#'   introducing Gaussian pseudo-data so that the posterior marginals match to the so called
#'    tilted distributions (leave-one-out posterior times the true likelihood factor) as
#'     closely as possible.  Typically more accurate than
#'   Laplace, but slower. }
#' }
#'
#' @name approx
#'
#' @param tol Convergence tolerance.
#' @param maxiter Maximum number of iterations in the Laplace/EP iteration.
#' @param damping Damping factor for EP. Should be between 0 and 1. Smaller values
#' typically lead to more stable iterations, but also increase the number of iterations,
#' and thus make the algorithm slower.
#' @param quad_order Order of the Gauss-Hermite quadrature used to evaluate the required
#'  tilted moments in EP.
#'
#'
#' @return The approximation object.
#'
#'
#' @section References:
#'
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning.
#' MIT Press.
#'
#'
#' @examples
#'
#' # Basic usage
#' gp <- gp_init(
#'   cfs = cf_sexp(),
#'   lik = lik_bernoulli(),
#'   method = method_fitc(num_inducing = 100),
#'   approx = approx_ep()
#' )
#'
NULL

# constructors

#' @rdname approx
#' @export
approx_laplace <- function(maxiter = 30, tol = 1e-4) {
  approx <- list(name = "laplace", maxiter = maxiter, tol = tol)
  class(approx) <- c("approx_laplace", "approx")
  approx
}

#' @rdname approx
#' @export
approx_ep <- function(damping = 0.9, quad_order = 11, maxiter = 100) {
  if (damping <= 0 || damping >= 1) {
    stop("Damping should be between 0 and 1 (exclusive)")
  }
  approx <- list(name = "ep", damping = damping, quad_order = quad_order, maxiter = maxiter)
  class(approx) <- c("approx_ep", "approx")
  approx
}

#' @export
print.approx <- function(x, quiet = FALSE, ...) {
  object <- x
  str <- class(object)[1]
  if (!quiet) {
    cat(str)
  }
  invisible(str)
}
