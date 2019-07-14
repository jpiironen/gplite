#' Fit a GP model
#'
#' Function \code{gp_fit} fits a GP model with the current hyperparameters. 
#' Notice that this function does not optimize the hyperparameters in any way, 
#' but only finds the Laplace approximation (or the analytical 
#' true posterior in the case of Gaussian likelihood) to the latent values. 
#' Function \code{gp_mcmc} draws from the posterior of the latent values 
#' given the current hyperparameter estimates using MCMC. For optimizing the hyperparameter
#' values, see \code{gp_optim}.
#' 
#' @name gp_fit
#' 
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input dimension). 
#' Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param trials Vecton of length n giving the number of trials for each observation in binomial 
#' (and beta binomial) model.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability.
#'  Default is 1e-6.
#' @param ... Further arguments to be passed to \link{rstan}'s function 
#' \code{\link[rstan]{optimizing}} (if \code{gp_fit} was called) or 
#' \code{\link[rstan]{sampling}} (if \code{gp_mcmc} was called).
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
#' # Analytic approximation
#' cf <- cf_sexp(lscale=0.3, magn=0.5)
#' lik <- lik_binomial()
#' gp <- gp_init(cf, lik)
#' gp <- gp_fit(gp, x, y)
#' 
#' # MCMC solution
#' gpmc <- gp_mcmc(gp, x, y, trials=trials, chains=2, iter=1000)
#' }
#'
NULL

#' @rdname gp_fit
#' @export
gp_fit <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  if (gp$method == 'full') {
    gp <- gp_laplace_full(gp, x, y, trials=trials, jitter=jitter, ...)
  } else if (gp$method == 'rf') {
    num_inputs <- NCOL(x)
    featuremap <- get_featuremap(gp, num_inputs)
    gp <- gp_laplace_linearized(gp, x, y, featuremap, trials=trials, jitter=jitter, ...)
  } else
    stop('Unknown method: ', gp$method)
  gp$fitted <- TRUE
  return(gp)
}


gp_laplace_full <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cf, x, x) + jitter*diag(n)
  gp$x <- x
  gp$K <- K
  gp$K_chol <- t(chol(K)) # lower triangular
  data <- c(list(n=n,L=gp$K_chol,y=y), get_standata(gp$lik, trials=trials))
  model <- get_stanmodel(gp$lik, gp$method)
  gp$fit <- rstan::optimizing(model, data=data, hessian=T, as_vector=F, init=0, ...)
  gp$fmean <- as.vector(gp$fit$par$f) # posterior mean for f
  fw_prec_chol <- t(chol(-as.matrix(gp$fit$hessian))) # cholesky of precision for the whitened f
  aux <- solve(t(gp$K_chol),fw_prec_chol)
  gp$fprec_chol <- t(chol(aux %*% t(aux))) # cholesky of precision for f
  gp$log_evidence <- gp$fit$value + 0.5*n*log(2*pi) - sum(log(diag(fw_prec_chol)))
  return(gp)
}

gp_laplace_linearized <- function(gp, x, y, featuremap, trials=NULL, jitter=NULL, ...) {
  gp$featuremap <- featuremap
  z <- gp$featuremap(x)
  n <- length(y)
  data <- c(list(n=n,m=ncol(z),Z=z,y=y), get_standata(gp$lik, trials=trials))
  model <- get_stanmodel(gp$lik, gp$method)
  gp$fit <- rstan::optimizing(model, data=data, hessian=T, as_vector=F, init=0, ...)
  gp$wmean <- as.vector(gp$fit$par$w)
  gp$wprec_chol <- t(chol(-as.matrix(gp$fit$hessian))) # cholesky of precision
  gp$log_evidence <- gp$fit$value + 0.5*n*log(2*pi) - sum(log(diag(gp$wprec_chol)))
  return(gp)
}




#' Optimize hyperparameters of a GP model
#' 
#' This function can be used to optimize the hyperparameters of the model to the maximum marginal
#' likelihood solution (type-II maximum likelihood) based on the Laplace approximation.
#' 
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input dimension). 
#' Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param method Optimization method that will be passed to \link{optim} function.
#' @param tol Relative change in the objective function value below the optimization is terminated. 
#' @param verbose If TRUE, then some information about the progress of the optimization is printed to the console.
#' @param ... Further arguments to be passed to \link{gp_fit} that are needed in the fitting process, for
#' example \code{trials} in the case of binomial likelihood.
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
gp_optim <- function(gp, x, y, method='Nelder-Mead', tol=1e-4, verbose=T, ...) {
  energy <- function(param) {
    gp <- set_param(gp, param)
    gp <- gp_fit(gp,x,y, ...)
    gp_optim_iter_message(gp, verbose)
    gp_energy(gp)
  }
  
  gp_optim_start_message(gp, verbose)
  param0 <- get_param(gp)
  res <- stats::optim(param0, energy, method=method, control = list(reltol=tol))
  param <- res$par
  gp <- set_param(gp, param)
  gp <- gp_fit(gp,x,y, ...)
  gp
}

gp_optim_start_message <- function(gp, verbose=T) {
  if (!verbose)
    return()
  nam <- names(get_param(gp))
  items <- sapply(seq_along(nam), function(i) paste0(sprintf('p%d: log ',i), nam[i]) )
  cat('Optimizing parameters\n')
  cat(paste(unname(items), collapse="\n"))
  cat('\n\n')
  
  symbols <- sapply(seq_along(nam), function(i) sprintf('p%d',i))
  row_items <- c(sprintf('%8s', symbols), sprintf('%10s\n','Energy'))
  cat(paste0(row_items, collapse=' '))
}

gp_optim_iter_message <- function(gp, verbose=T) {
  if (!verbose)
    return()
  row_items <- c(sprintf('%8.2f', get_param(gp)), sprintf('%10.2f', gp_energy(gp)))
  cat(paste0(row_items, collapse=' '), '\n')
}


