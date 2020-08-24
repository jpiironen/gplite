#' Model assessment and comparison
#'
#' Function \code{gp_loo} computes the approximate leave-one-out (LOO) 
#' cross-validation statistics for the given GP model with the current 
#' hyperparameters. 
#' Function \code{gp_compare} estimates the difference in the expected 
#' predictive accuracy of two or more GP models given their LOO statistics.
#' 
#' @name gp_loo
#' 
#' @param gp The gp model object to be fitted.
#' @param x n-by-d matrix of input values (n is the number of observations and d the input dimension). 
#' Can also be a vector of length n if the model has only a single input.
#' @param y Vector of n output (target) values.
#' @param trials Vector of length n giving the number of trials for each observation in binomial 
#' (and beta binomial) model.
#' @param draws Number of posterior draws to estimate the required integrals.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability.
#'  Default is 1e-6.
#' @param seed Random seed.
#' @param ... LOO statistics for the models to compare.
#'
#'
#' @return \code{gp_loo} returns a list with LOO statistics. 
#' \code{gp_compare} returns a matrix with comparison statistics.
#'
#' @section References:
#' 
#' Vehtari A., Mononen T., Tolvanen V., Sivula T. and Winther O (2016). 
#' Bayesian Leave-One-Out Cross-Validation Approximations for Gaussian Latent
#' Variable Models. Journal of Machine Learning Research 17(103):1-38.
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
#' # Set up two models
#' gp1 <- gp_init(cf_sexp(), lik_binomial())
#' gp2 <- gp_init(cf_matern32(), lik_binomial())
#' 
#' # Optimize
#' gp1 <- gp_optim(gp1, x, y, trials=trials)
#' gp2 <- gp_optim(gp2, x, y, trials=trials)
#' 
#' # Compare
#' loo1 <- gp_loo(gp1, x, y, trials=trials)
#' loo2 <- gp_loo(gp2, x, y, trials=trials)
#' gp_compare(loo1, loo2)
#' 
#' }
#'
NULL


#' @rdname gp_loo
#' @export
gp_loo <- function(gp, x, y, trials=NULL, draws=4000, jitter=NULL, seed=NULL) {
  
  
  fhat <- as.vector(gp$fit$fmean)
  pobs <- get_pseudodata(gp$lik, fhat, y, trials=trials)
  z <- pobs$z
  V <- pobs$var
  grad <- (z - fhat)/V
  pred_var <- gp_pred(gp, x, var=TRUE)$var
  
  # mean and variance of LOO posteriors
  loo_var <- 1/(1/pred_var - 1/V)
  loo_mean <- fhat - loo_var*grad
  
  # sample from LOO posteriors and evaluate predictive distribution in using Monte Carlo
  n <- length(y)
  fsample <- matrix(stats::rnorm(n*draws, mean = loo_mean, sd = sqrt(loo_var)), nrow=n) 
  loglik <- get_loglik(gp$lik, fsample, y, trials=trials, sum=FALSE)
  
  loos <- apply(loglik, 1, logsumexp) - log(draws)
  
  res <- list(loo=sum(loos), sd=stats::sd(loos), loos=loos)
  class(res) <- 'loores'
  return(res)
  
}


logsumexp <- function(x) {
  c <- max(x)
  c + log(sum(exp(x - c)))
}


#' @rdname gp_loo
#' @export
gp_compare <- function(...) {
  args <- list(...)
  loo <- c()
  loos <- list()
  for (i in seq_along(args)) {
    if (!('loores' %in% class(args[[i]]))) {
      stop('Expected objects of type \'loores\', but found an object with type \'', class(args[[i]]), '\'')
    }
    
    loo[i] <- args[[i]]$loo
    loos[[i]] <- args[[i]]$loos
  }
  
  imax <- which.max(loo)
  
  res <- matrix(nrow=length(loos), ncol=2)
  for (i in seq_along(loos)) {
    d <- loos[[i]] - loos[[imax]]
    res[i,1] <- sum(d)
    res[i,2] <- stats::sd(d)
  }
  rownames(res) <- sapply(1:length(loos), function(i) sprintf('model%i', i))
  colnames(res) <- c('loo-diff', 'se')
  return(res)
}
