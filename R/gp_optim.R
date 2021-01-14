#' Optimize hyperparameters of a GP model
#'
#' This function can be used to optimize the hyperparameters of the model to the maximum
#' marginal likelihood (or maximum marginal posterior if priors are used), using Nelder-Mead
#' algorithm.
#'
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input
#'  dimension). Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param tol Relative change in the objective function value after which the optimization is
#'  terminated.
#' @param maxiter Maximum number of iterations.
#' @param verbose If TRUE, then some information about the progress of the optimization is
#'  printed to the console.
#' @param warnings Whether to print out some potential warnings (such as maximum number of
#' iterations reached) during the optimization.
#' @param ... Further arguments to be passed to \code{\link{gp_fit}} that are needed
#' in the fitting process, for example \code{trials} in the case of binomial likelihood.
#'
#'
#' @return An updated GP model object.
#'
#' @section References:
#'
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#'
#' # Generate some toy data
#' set.seed(1242)
#' n <- 50
#' x <- matrix(rnorm(n * 3), nrow = n)
#' f <- sin(x[, 1]) + 0.5 * x[, 2]^2 + x[, 3]
#' y <- f + 0.5 * rnorm(n)
#' x <- data.frame(x1 = x[, 1], x2 = x[, 2], x3 = x[, 3])
#'
#' # Basic usage
#' cf <- cf_sexp()
#' lik <- lik_gaussian()
#' gp <- gp_init(cf, lik)
#' gp <- gp_optim(gp, x, y)
#' 
#'
#' @export
gp_optim <- function(gp, x, y, tol = 1e-4, maxiter = 500, verbose = TRUE, warnings = TRUE, ...) {
  iter <- 0
  energy <- function(param) {
    gp <- set_param(gp, param)
    gp <- gp_fit(gp, x, y, ...)
    optim_iter_message(gp, iter, verbose)
    iter <<- iter + 1
    gp_energy(gp)
  }

  optim_start_message(gp, verbose)
  param0 <- get_param(gp)
  n <- length(y)
  control <- list(
    reltol = tol, fnscale = n, maxit = maxiter,
    warn.1d.NelderMead = warnings
  )
  res <- stats::optim(param0, energy, method = "Nelder-Mead", control = control)
  if (res$convergence == 1 && warnings) {
    warning("Maximum number of iterations reached, optimization may not have converged.")
  }
  param_opt <- res$par
  e_opt <- res$value

  if (warnings) {
    check_convergence(energy, e_opt, param_opt, verbose = verbose)
  }

  gp <- set_param(gp, param_opt)
  gp <- gp_fit(gp, x, y, ...)
  gp
}

check_convergence <- function(lossfun, loss_opt, param_opt, verbose = FALSE) {

  # check convergence
  delta <- 1e-1
  tol <- 1e-2
  if (verbose) {
    cat("Assessing convergence...\n")
  }
  for (k in seq_along(param_opt)) {
    dparam <- rep(0, length(param_opt))
    dparam[k] <- delta
    loss1 <- lossfun(param_opt + dparam)
    loss2 <- lossfun(param_opt - dparam)
    loss_diff <- loss_opt - min(loss1, loss2)
    if (loss_diff > tol) {
      warning(sprintf("Not all hyperparameters have reached convergence within tolerance %.2f. Try reoptimizing starting from the current hyperparameter values, or reduce tol.", delta))
      break
    }
  }
}

optim_start_message <- function(gp, verbose = TRUE) {
  if (!verbose) {
    return()
  }
  nam <- names(get_param(gp))
  items <- sapply(seq_along(nam), function(i) paste0(sprintf("p%d: log ", i), nam[i]))
  cat("Optimizing parameters\n")
  cat(paste(unname(items), collapse = "\n"))
  cat("\n\n")

  symbols <- sapply(seq_along(nam), function(i) sprintf("p%d", i))
  row_items <- c(sprintf("%8s", symbols), sprintf("%10s", "Energy"), sprintf("%9s\n", "Iteration"))
  cat(paste0(row_items, collapse = " "))
}

optim_iter_message <- function(gp, iter, verbose = TRUE) {
  if (!verbose) {
    return()
  }
  row_items <- c(sprintf("%8.2f", get_param(gp)), sprintf("%10.2f", gp_energy(gp)), sprintf("%9d", iter))
  cat(paste0(row_items, collapse = " "), "\n")
}
