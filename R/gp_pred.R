#' Make predictions with a GP model
#' 
#' This function can be used to make predictions with a fitted model.
#' 
#' @param gp A fitted GP model object.
#' @param xt N-by-d matrix of input values (N is the number of test points and d the input dimension). 
#' Can also be a vector of length N if the model has only a single input.
#' @param var Whether to compute the predictive variances along with predictive mean
#' @param draws Number of draws to generate from the predictive distribution for the latent values. 
#' @param transform Whether to transform the draws of latent values to the same scale as the target y.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability. Default is 1e-4.
#'
#'
#' @return If draws is given (sampling from the predictive distribution), then returns an N-by-draws
#' matrix of random draws from the predictive distribution. If draws is not set, then returns either only
#' the predictive means (if \code{var=F}) or a list with both the predictive means and variances (if \code{var=T}).
#'  
#' @section References:
#' 
#' Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian processes for machine learning. MIT Press.
#'
#' @examples
#' \donttest{
#' # fit GP
#' gp <- gp_init()
#' gp <- gp_optim(gp, x, y, trials=trials)
#' 
#' # analytic prediction
#' pred <- gp_pred(gp, xt, var=T)
#' pred$mean 
#' pred$var 
#' 
#' # draw from the predictive distribution based on the analytic posterior approximation
#' draws <- gp_pred(gp, xt, draws=1000) 
#' 
#' # fit using mcmc and predict using that
#' gpmc <- gp_sample(gp, x, y, trials=trials, iter=1000, chains=2)
#' draws <- gp_pred(gpmc, xt, draws=1000)
#' }
#'
#'
#' @export
gp_pred <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=NULL) {
  
  if (is_fitted(gp, 'sampling')) {
    # model fitted using mcmc, so predict using the draws from the posterior
    if (gp$method == 'full')
      pred <- gp_pred_full_mcmc(gp, xt, draws=draws, transform=transform, jitter=jitter)
    else if (gp$method == 'rff')
      pred <- gp_pred_linearized_mcmc(gp, xt, draws=draws, transform=transform, jitter=jitter)
    else
      stop('Unknown method: ', gp$method)
  } else {
    # model fitted using analytical approximation, so predict based on that
    if (!is_fitted(gp, 'analytic')) 
      stop('The given model appears not to be fitted yet.')
    if (gp$method == 'full')
      pred <- gp_pred_full_analytic(gp, xt, var=var, draws=draws, 
                                    transform=transform, jitter=jitter)
    else if (gp$method == 'rff')
      pred <- gp_pred_linearized_analytic(gp, xt, var=var, draws=draws,
                                          transform=transform, jitter=jitter)
    else
      stop('Unknown method: ', gp$method)
  }
  return(pred)
}

gp_pred_full_mcmc <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {
  
  draws <- NCOL(gp$fsample)
  nt <- NROW(xt)
  jitter <- get_jitter(gp,jitter)
  Kt <- eval_cf(gp$cf, xt, gp$x)
  Ktt <- eval_cf(gp$cf, xt, xt)
  K_chol <- gp$K_chol
  aux <- solve(K_chol, t(Kt))
  pred_cov <- Ktt - t(aux) %*% aux + jitter*diag(nt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, gp$fsample))
  sample <- t(chol(pred_cov)) %*% matrix(stats::rnorm(nt*draws), nrow=nt) + pred_mean
  if (transform)
    sample <- get_response(gp$lik, sample)
  return(sample)
}

gp_pred_linearized_mcmc <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {
  
  zt <- gp$featuremap(xt)
  sample <- zt %*% gp$wsample
  if (transform)
    sample <- get_response(gp$lik, sample)
  return(sample)
}



gp_pred_full_analytic <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=NULL) {
  
  # compute the latent mean first
  K_chol <- gp$K_chol
  Kt <- eval_cf(gp$cf, xt, gp$x)
  Ktt <- eval_cf(gp$cf, xt, xt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, gp$fmean))
  pred_mean <- as.vector(pred_mean)
  
  if (var == T || !is.null(draws)) {
    # covariance of the latent function
    nt <- length(pred_mean)
    jitter <- get_jitter(gp,jitter)
    aux1 <- solve(K_chol, t(Kt))
    aux2 <- solve(gp$fprec_chol, solve(t(K_chol), solve(K_chol, t(Kt))))
    pred_cov <- Ktt - t(aux1) %*% aux1 + t(aux2) %*% aux2 + jitter*diag(nt)
    
    if (is.null(draws))
      return(list(mean=pred_mean, var=diag(pred_cov)))
    else {
      # sample from the predictive distribution
      sample <- t(chol(pred_cov)) %*% matrix(stats::rnorm(nt*draws), nrow=nt) + pred_mean
      if (transform)
        sample <- get_response(gp$lik, sample)
      return(sample)
    }
  }
  return(pred_mean)
}

gp_pred_linearized_analytic <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=NULL) {
  
  # compute the latent mean first
  zt <- gp$featuremap(xt)
  pred_mean <- as.vector(zt %*% gp$wmean)
  
  if (var == T || !is.null(draws)) {
    # covariance of the latent function
    nt <- length(pred_mean)
    jitter <- get_jitter(gp,jitter) 
    aux <- solve(gp$wprec_chol, t(zt))
    pred_cov <- t(aux) %*% aux + jitter*diag(nt)
    
    if (is.null(draws))
      return(list(mean=pred_mean, var=diag(pred_cov)))
    else {
      # sample from the predictive distribution
      sample <- t(chol(pred_cov)) %*% matrix(stats::rnorm(nt*draws), nrow=nt) + pred_mean
      if (transform)
        sample <- get_response(gp$lik, sample)
      return(sample)
    }
  }
  return(pred_mean)
}






