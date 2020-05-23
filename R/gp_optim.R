#' Optimize hyperparameters of a GP model
#' 
#' This function can be used to optimize the hyperparameters of the model to the maximum
#'  marginal
#' likelihood solution (type-II maximum likelihood) based on the Laplace approximation.
#' 
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input
#'  dimension). 
#' Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param method Optimization method that will be passed to \code{\link{optim}} function.
#' @param tol Relative change in the objective function value below the optimization is
#'  terminated. 
#' @param maxiter Maximum number of iterations.
#' @param verbose If TRUE, then some information about the progress of the optimization is
#'  printed to the console.
#' @param warnings Whether to print out some potential warnigns (such as maximum number of 
#' iterations reached) during the optimization.
#' @param ... Further arguments to be passed to \code{\link{gp_fit}} that are needed in the fitting
#' process, for example \code{trials} in the case of binomial likelihood.
#'
#'
#' @return An updated GP model object.
#'  
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' # Basic usage (single covariance function)
#' cf <- cf_sexp()
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x, y, trials=trials)
#' 
#' }
#'
#'
#' 
#' @export
gp_optim <- function(gp, x, y, method='Nelder-Mead', tol=1e-4, 
                     maxiter=500, verbose=T, warnings=T, ...) {

  iter <- 0
  energy <- function(param) {
    gp <- set_param(gp, param)
    gp <- gp_fit(gp,x,y, ...)
    optim_iter_message(gp, iter, verbose)
    iter <<- iter + 1
    gp_energy(gp)
  }
  
  optim_start_message(gp, verbose)
  param0 <- get_param(gp)
  control <- list(reltol = tol, fnscale = length(y), maxit = maxiter, 
                  warn.1d.NelderMead = warnings)
  res <- stats::optim(param0, energy, method = method, control = control)
  if (res$convergence == 1 && warnings)
    warning('Maximum number of iterations reached, the optimization may not have converged.')
  param <- res$par
  gp <- set_param(gp, param)
  gp <- gp_fit(gp,x,y, ...)
  gp
}

optim_start_message <- function(gp, verbose=T) {
  if (!verbose)
    return()
  nam <- names(get_param(gp))
  items <- sapply(seq_along(nam), function(i) paste0(sprintf('p%d: log ',i), nam[i]) )
  cat('Optimizing parameters\n')
  cat(paste(unname(items), collapse="\n"))
  cat('\n\n')
  
  symbols <- sapply(seq_along(nam), function(i) sprintf('p%d',i))
  row_items <- c(sprintf('%8s', symbols), sprintf('%10s','Energy'), sprintf('%9s\n', 'Iteration'))
  cat(paste0(row_items, collapse=' '))
}

optim_iter_message <- function(gp, iter, verbose=T) {
  if (!verbose)
    return()
  row_items <- c(sprintf('%8.2f', get_param(gp)), sprintf('%10.2f', gp_energy(gp)), sprintf('%9d', iter))
  cat(paste0(row_items, collapse=' '), '\n')
}

