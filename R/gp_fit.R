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
#' @param trials Vector of length n giving the number of trials for each observation in binomial 
#' (and beta binomial) model.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability.
#'  Default is 1e-6.
#' @param ... Further arguments to be passed to \link{rstan}'s function 
#' \code{\link[rstan]{sampling}} if \code{gp_mcmc} was called.
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
#' 
#' # Generate some toy data
#' set.seed(32004)
#' n <- 150
#' sigma <- 0.1
#' x <- rnorm(n)
#' ycont <- sin(3*x)*exp(-abs(x)) +  rnorm(n)*sigma
#' y <- rep(0,n)
#' y[ycont > 0] <- 1
#' trials <- rep(1,n)
#' 
#' # Analytic approximation (Laplace)
#' cf <- cf_sexp(lscale=0.3, magn=3)
#' gp <- gp_init(cf, lik_binomial())
#' gp <- gp_fit(gp, x, y, trials=trials)
#' 
#' # MCMC solution
#' gpmc <- gp_mcmc(gp, x, y, trials=trials, chains=2, iter=1000)
#' }
#'
NULL

#' @rdname gp_fit
#' @export
gp_fit <- function(gp, x, y, trials=NULL, jitter=NULL, ...) {
  gp <- learn_scales(gp, x)
  gp <- fit_model(gp, x, y, trials=trials, jitter=jitter, ...)
  gp$fit$type <- 'analytic'
  gp$fitted <- TRUE
  return(gp)
}

fit_model <- function(object, ...) {
  UseMethod("fit_model", object)
}

fit_model.gp <- function(object, ...) {
  fit_model(object$latent, object, ...)
}

fit_model.latent_laplace <- function(object, gp, ...) {
  fit_laplace(gp$approx, gp, ...)
}

fit_model.latent_ep <- function(object, gp, ...) {
  fit_ep(gp$approx, gp, ...)
}

fit_laplace <- function(object, ...) {
  UseMethod("fit_laplace", object)
}

fit_ep <- function(object, ...) {
  UseMethod("fit_ep", object)
}

fit_ep.gp <- function(object, ...) {
  fit_ep(object$approx, object, ...)
}


fit_laplace.approx_full <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cfs, x, x) + jitter*diag(n)
  gp$x <- x
  gp$fit <- laplace(object, gp, K, y, trials=trials)
  return(gp)
}

fit_laplace.approx_fitc <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  z <- get_inducing(gp, x)
  Kz <- eval_cf(gp$cfs, z, z) + jitter*diag(gp$approx$num_inducing)
  Kxz <- eval_cf(gp$cfs, x, z)
  Kz_chol <- t(chol(Kz))
  D <- eval_cf(gp$cfs, x, x, diag_only=T) + jitter
  D <- D - colSums(forwardsolve(Kz_chol, t(Kxz))^2)
  gp$x <- x
  gp$x_inducing <- z
  gp$fit <- tryCatch({
    laplace(object, gp, Kz, Kz_chol, Kxz, D, y, trials=trials)
  },
  error = function(err) {
    print(err)
    list(log_evidence=-Inf)
  })
  return(gp)
}

fit_laplace.approx_rf <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  num_inputs <- NCOL(x)
  featuremap <- get_featuremap(gp, num_inputs)
  gp$approx$num_basis <- check_num_basis(gp$cfs, gp$approx$num_basis, NCOL(x))
  z <- featuremap(x)
  gp$fit <- laplace(object, gp, z, y, trials=trials)
  return(gp)
}

fit_laplace.approx_rbf <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  gp$x <- x
  num_inputs <- NCOL(x)
  featuremap <- get_featuremap(gp, num_inputs)
  gp$approx$num_basis <- check_num_basis(gp$cfs, gp$approx$num_basis, NCOL(x))
  z <- featuremap(x)
  gp$fit <- laplace(object, gp, z, y, trials=trials)
  return(gp)
}



fit_ep.approx <- function(object, ...) {
  stop(paste('No EP implementation for', class(object)[1], 'yet.'))
}

fit_ep.approx_full <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  K <- eval_cf(gp$cfs, x, x) + jitter*diag(n)
  gp$x <- x
  gp$fit <- ep(object, gp, K, y, trials=trials, quad_order=gp$latent$quad_order,
               damping=gp$latent$damping, damping_min=0.1, ...)
  return(gp)
}

fit_ep.approx_fitc <- function(object, gp, x, y, trials=NULL, jitter=NULL, ...) {
  x <- as.matrix(x)
  n <- length(y)
  jitter <- get_jitter(gp,jitter)
  z <- get_inducing(gp, x)
  Kz <- eval_cf(gp$cfs, z, z) + jitter*diag(gp$approx$num_inducing)
  Kxz <- eval_cf(gp$cfs, x, z)
  Kz_chol <- t(chol(Kz))
  K_diag <- eval_cf(gp$cfs, x, x, diag_only=T) + jitter
  D <- K_diag - colSums(forwardsolve(Kz_chol, t(Kxz))^2)
  gp$x <- x
  gp$x_inducing <- z
  gp$fit <- tryCatch({
    ep(object, gp, Kz, Kz_chol, Kxz, K_diag, D, y, trials=trials, ...)
  },
  error = function(err) {
    print(err)
    list(log_evidence=-Inf)
  })
  return(gp)
}



get_inducing <- function(gp, x) {
  # set random seed but ensure the old RNG state is restored on exit
  if (exists('.Random.seed')) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(gp$approx$seed)
  
  if (!is.null(gp$x_inducing))
    return(gp$x_inducing)
  xscaled <- scale(x)
  binning <- gp$approx$bin_inducing
  num_inducing <- gp$approx$num_inducing
  
  if (!is.null(binning)) {
    varname <- names(binning)[1] # TODO: currently handles only one variable
    nbins <- binning[[varname]]
    xscaled_binned <- bin(xscaled, nbins=nbins, var=varname)
    bin_sizes <- sapply(xscaled_binned, function(bin) NROW(bin))
    ordering <- order(bin_sizes, decreasing = T)
    xscaled_binned <- xscaled_binned[ordering]
    bin_sizes <- bin_sizes[ordering]
    per_bin <- rep(floor(num_inducing/nbins), nbins)
    num_leftover <- num_inducing - sum(per_bin)
    if (num_leftover > 0)
      per_bin[1:num_leftover] <- per_bin[1:num_leftover] + 1
    zscaled_binned <- lapply(1:nbins, function(i) {
      kmeans(xscaled_binned[[i]], per_bin[i])$centers
    })
    zscaled <- do.call(rbind, zscaled_binned)
    rownames(zscaled) <- NULL
  } else {
    cl <- stats::kmeans(xscaled, num_inducing)
    zscaled <- cl$centers
    rownames(zscaled) <- NULL
  }
  z <- t( t(zscaled)*attr(xscaled, 'scaled:scale') + attr(xscaled, 'scaled:center') )
  return(z)
}


bin <- function(x, nbins=NULL, cutpoints=NULL, var=1) {
  
  x <- as.matrix(x)
  xv_min <- min(x[,var]) - 2*.Machine$double.eps
  xv_max <- max(x[,var]) + 2*.Machine$double.eps
  
  if (is.null(cutpoints))
    if (!is.null(nbins))
      cutpoints <- seq(xv_min, xv_max, len=nbins+1)
  else
    stop('Either bins or cutpoints must be provided')
  else
    cutpoints <- c(xv_min, cutpoints, xv_max)
  
  nbins <- length(cutpoints) - 1
  xbinned <- lapply(1:nbins, function(i) {
    ind <- x[,var] > cutpoints[i] & x[,var] < cutpoints[i+1]
    x[ind,,drop=F]
  })
  xbinned
}










