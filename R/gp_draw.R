




#' @rdname pred
#' @export
gp_draw <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {

  if (is_fitted(gp, 'sampling')) {
    #
    # model fitted using mcmc, so predict using the draws from the posterior
    #
    if (gp$method == 'full')
      pred <- gp_draw_full_mcmc(gp, xt, draws=draws, transform=transform, jitter=jitter)
    else if (gp$method == 'rf')
      pred <- gp_draw_linearized_mcmc(gp, xt, draws=draws, transform=transform)
    else
      stop('Unknown method: ', gp$method)
  } else {
    #
    # model fitted using analytical gaussian approximation (or not fitted at all),
    # so predict based on that
    #
    if (is_fitted(gp, 'analytic')) {
      # draw from the analytical posterior approximation
      if (gp$method == 'full')
        pred <- gp_draw_full_analytic(gp, xt, draws=draws,
                                      transform=transform, jitter=jitter)
      else if (gp$method == 'rf')
        pred <- gp_draw_linearized_analytic(gp, xt, draws=draws,
                                            transform=transform, jitter=jitter)
      else
        stop('Unknown method: ', gp$method)
    } else {
      # draw from the prior
      if (gp$method == 'full')
        pred <- gp_draw_full_prior(gp, xt, draws=draws,
                                   transform=transform, jitter=jitter)
      else if (gp$method == 'rf')
        pred <- gp_draw_linearized_prior(gp, xt, draws=draws,
                                         transform=transform, jitter=jitter)
      else
        stop('Unknown method: ', gp$method)
    }
  }
  return(pred)
}


gp_draw_full_mcmc <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {

  # TODO: this ignores draws, and sets it equal to the number of posterior draws
  
  draws <- NCOL(gp$fsample)
  nt <- NROW(xt)
  jitter <- get_jitter(gp,jitter)
  Kt <- eval_cf(gp$cf, xt, gp$x)
  Ktt <- eval_cf(gp$cf, xt, xt)
  K_chol <- gp$K_chol
  aux <- solve(K_chol, t(Kt))
  pred_cov <- Ktt - t(aux) %*% aux + jitter*diag(nt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, gp$fsample))
  sample <- mvnrnd(draws, pred_mean, chol_cov = t(chol(pred_cov)))
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_full_analytic <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {
  
  pred <- gp_pred_full_post(gp, xt, cov=T, transform=F, jitter=jitter)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_full_prior <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {

  pred <- gp_pred_full_prior(gp, xt, cov=T, transform=F, jitter=jitter)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}


gp_draw_linearized_mcmc <- function(gp, xt, draws=NULL, transform=T) {

  # TODO: this ignores draws, and sets it equal to the number of posterior draws
  
  zt <- gp$featuremap(xt)
  sample <- zt %*% gp$wsample
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_linearized_analytic <- function(gp, xt, draws=NULL, transform=T, jitter=NULL) {

  # draw from the posterior of w
  zt <- gp$featuremap(xt)
  w <- mvnrnd(draws, gp$wmean, chol_prec=gp$wprec_chol)
  sample <- zt %*% w
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_linearized_prior <- function(gp, xt, var=F, draws=NULL, transform=T, jitter=NULL) {

  # draw from the prior of w
  num_inputs <- NCOL(xt)
  featuremap <- get_featuremap(gp, num_inputs)
  zt <- featuremap(xt)
  num_feat <- NCOL(zt)
  w <- matrix(stats::rnorm(num_feat*draws), nrow=num_feat)
  sample <- zt %*% w
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

mvnrnd <- function(draws, mean, chol_cov=NULL, chol_prec=NULL) {
  # Draw from multivariate normal. Lower Cholesky factor of either
  # covariance or precision must be provided.
  d <- NROW(mean)
  if (!is.null(chol_cov))
    r <- as.vector(mean) + chol_cov %*% matrix(stats::rnorm(d*draws), nrow=d)
  else if (!is.null(chol_prec))
    r <- as.vector(mean) + solve(t(chol_prec), matrix(stats::rnorm(d*draws), nrow=d))
  else
    stop('Both cov_chol and prec_chol can\'t be NULL.')
  return(r)
}



