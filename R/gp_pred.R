#' Make predictions with a GP model
#' 
#' Function \code{gp_pred} can be used to make analytic predictions for the latent function
#' values at test points, whereas \code{gp_draw}
#' can be used to draw from the predictive distribution (or from the prior if the GP has
#' not been fitted yet.)
#' 
#' @name pred
#' 
#' @param gp A fitted GP model object.
#' @param xnew N-by-d matrix of input values (N is the number of test points and d 
#' the input dimension). 
#' Can also be a vector of length N if the model has only a single input.
#' @param var Whether to compute the predictive variances along with predictive mean.
#' @param draws Number of draws to generate from the predictive distribution for the 
#' latent values. 
#' @param transform Whether to transform the draws of latent values to the same scale
#'  as the target y.
#' @param target If TRUE, draws values for the target variable \code{y} instead of the latent
#'  function values.
#' @param cfind Indices of covariance functions to be used in the prediction. By default uses
#' all covariance functions.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability.
#'  Default is 1e-6.
#' @param seed Random seed for draws.
#'
#'
#' @return \code{gp_pred} returns a vector of predictive mean (one value for each row of
#'  \code{xnew}) for the latent mean, or a list with fields having both the 
#'  mean and variance for each 
#'  observation if \code{var = TRUE}. \code{gp_draw} returns an N-by-draws
#' matrix of random draws from the predictive distribution. 
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
#' pred <- gp_pred(gp, xnew, var=T)
#' pred$mean 
#' pred$var 
#' 
#' # draw from the predictive distribution based on the analytic posterior approximation
#' draws <- gp_pred(gp, xnew, draws=1000) 
#' 
#' # fit using mcmc and draw from the predictive distribution
#' gpmc <- gp_mcmc(gp, x, y, trials=trials, iter=1000, chains=2)
#' draws <- gp_draw(gpmc, xnew)
#' }
#'
NULL

#' @rdname pred
#' @export
gp_pred <- function(gp, xnew, var=F, cfind=NULL, jitter=NULL) {
  
  if (is_fitted(gp, 'sampling')) 
    stop('Only gp_draw currently available for models fitted using gp_mcmc.')

  if (!is_fitted(gp, 'analytic')) {
    # model not fitted, so predict based on the prior
    if (gp$method == 'full')
      pred <- gp_pred_full_prior(gp, xnew, var=var, cfind=cfind, jitter=jitter)
    else if (gp$method == 'rf')
      pred <- gp_pred_linearized_prior(gp, xnew, var=var, cfind=cfind, jitter=jitter)
    else
      stop('Unknown method: ', gp$method)
  } else {
    # model fitted using analytical approximation
    if (gp$method == 'full')
      pred <- gp_pred_full_post(gp, xnew, var=var, cfind=cfind, jitter=jitter)
    else if (gp$method == 'rf')
      pred <- gp_pred_linearized_post(gp, xnew, var=var, cfind=cfind, jitter=jitter)
    else
      stop('Unknown method: ', gp$method)
  }
  return(pred)
}


gp_pred_full_post <- function(gp, xt, var=F, cov=F, cfind=NULL, jitter=NULL) {
  
  # compute the latent mean first
  Kt <- eval_cf(gp$cfs, xt, gp$x, cfind)
  Ktt <- eval_cf(gp$cfs, xt, xt, cfind)
  K_chol <- gp$fit$K_chol
  fmean <- gp$fit$fmean
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, fmean))
  pred_mean <- as.vector(pred_mean)
  
  if (var || cov) {
    # (co)variance of the latent function
    nt <- length(pred_mean)
    jitter <- get_jitter(gp,jitter)
    C_chol <- gp$fit$C_chol
    aux <- solve(C_chol, t(Kt))
    pred_cov <- Ktt - t(aux) %*% aux
    if (cov)
      return(list(mean = pred_mean, cov = pred_cov + jitter*diag(nt)))
    else
      return(list(mean = pred_mean, var = diag(pred_cov)))
  }
  return(pred_mean)
}


gp_pred_full_prior <- function(gp, xt, var=F, cov=F, cfind=NULL, jitter=NULL) {
  
  nt <- NROW(xt)
  pred_mean <- rep(0,nt)
  
  if (var || cov) {
    jitter <- get_jitter(gp, jitter)
    pred_cov <- eval_cf(gp$cfs, xt, xt, cfind)
    if (cov)
      return(list(mean = pred_mean, cov = pred_cov + jitter*diag(nt)))
    else
      return(list(mean = pred_mean, var = diag(pred_cov)))
  }
  return(pred_mean)
}



gp_pred_linearized_post <- function(gp, xt, var=F, cfind=NULL, jitter=NULL) {
  
  # compute the latent mean first
  featuremap <- get_featuremap(gp, num_inputs = NCOL(xt))
  zt <- featuremap(xt, cfind)
  wmean <- get_w_mean(gp, cfind)
  pred_mean <- as.vector(zt %*% wmean)
  
  if (var == T) {
    # covariance of the latent function
    nt <- length(pred_mean)
    wcov <- get_w_cov(gp, cfind)
    pred_cov <- zt %*% (wcov %*% t(zt))
    return(list(mean=pred_mean, var=diag(pred_cov)))
  }
  return(pred_mean)
}

gp_pred_linearized_prior <- function(gp, xt, var=F, cfind=NULL, jitter=NULL) {
  
  # mean is zero
  nt <- NROW(xt)
  pred_mean <- rep(0,nt)
  
  if (var == T) {
    num_inputs <- NCOL(xt)
    featuremap <- get_featuremap(gp, num_inputs)
    zt <- featuremap(xt, cfind)
    pred_cov <- zt %*% t(zt) 
    return(list(mean=pred_mean, var=diag(pred_cov)))
  }
  return(pred_mean)
}





