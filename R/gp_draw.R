




#' @rdname pred
#' @export
gp_draw <- function(gp, xnew, draws=NULL, transform=T, cfind=NULL, jitter=NULL, seed=NULL) {
  
  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(seed)

  if (is_fitted(gp, 'sampling')) {
    #
    # model fitted using mcmc, so predict using the draws from the posterior
    #
    if (gp$method == 'full')
      pred <- gp_draw_full_mcmc(gp, xnew, draws=draws, transform=transform, 
                                cfind=cfind, jitter=jitter)
    else if (gp$method == 'rf')
      pred <- gp_draw_linearized_mcmc(gp, xnew, draws=draws, transform=transform)
    else
      stop('Unknown method: ', gp$method)
  } else {
    #
    # model fitted using analytical gaussian approximation (or not fitted at all),
    # so predict based on that
    #
    if (is.null(draws))
      stop('Please provide the number of draws.')
    if (is_fitted(gp, 'analytic')) {
      # draw from the analytical posterior approximation
      if (gp$method == 'full')
        pred <- gp_draw_full_analytic(gp, xnew, draws=draws, transform=transform,
                                      cfind=cfind, jitter=jitter)
      else if (gp$method == 'rf')
        pred <- gp_draw_linearized_analytic(gp, xnew, draws=draws,
                                            transform=transform, jitter=jitter)
      else
        stop('Unknown method: ', gp$method)
    } else {
      # draw from the prior
      if (gp$method == 'full')
        pred <- gp_draw_full_prior(gp, xnew, draws=draws, transform=transform,
                                   cfind=cfind, jitter=jitter)
      else if (gp$method == 'rf')
        pred <- gp_draw_linearized_prior(gp, xnew, draws=draws,
                                         transform=transform, jitter=jitter)
      else
        stop('Unknown method: ', gp$method)
    }
  }
  return(pred)
}


gp_draw_full_mcmc <- function(gp, xt, draws=NULL, transform=T, cfind=NULL, jitter=NULL) {

  fsample <- gp$fsample
  if (is.null(draws))
    draws <- NCOL(fsample)
  else if (draws > NCOL(fsample))
    stop('Can\'t draw more than the provided GP has posterior draws.')
  # permute the posterior draws and pick right number of draws
  fsample <- fsample[,sample(1:NCOL(fsample))]
  fsample <- fsample[,1:draws]
  
  nt <- NROW(xt)
  jitter <- get_jitter(gp,jitter)
  Kt <- eval_cf(gp$cfs, xt, gp$x, cfind)
  Ktt <- eval_cf(gp$cfs, xt, xt, cfind)
  K_chol <- gp$K_chol
  aux <- solve(K_chol, t(Kt))
  pred_cov <- Ktt - t(aux) %*% aux + jitter*diag(nt)
  pred_mean <- Kt %*% solve(t(K_chol), solve(K_chol, fsample))
  sample <- mvnrnd(draws, pred_mean[,1:draws], chol_cov = t(chol(pred_cov)))
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_full_analytic <- function(gp, xt, draws=NULL, transform=T, cfind=NULL, jitter=NULL) {
  
  pred <- gp_pred_full_post(gp, xt, cov=T, cfind=cfind, jitter=jitter)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_full_prior <- function(gp, xt, draws=NULL, transform=T, cfind=NULL, jitter=NULL) {

  pred <- gp_pred_full_prior(gp, xt, cov=T, cfind=cfind, jitter=jitter)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}


gp_draw_linearized_mcmc <- function(gp, xt, draws=NULL, transform=T, cfind=NULL) {
  
  # TODO: take into account cfind
  if (!is.null(cfind))
    stop('cfind not supported for method rf yet.')

  wsample <- gp$wsample
  if (is.null(draws))
    draws <- NCOL(wsample)
  else if (draws > NCOL(wsample))
    stop('Can\'t draw more than the provided GP has posterior draws.')
  # permute the posterior draws and pick right number of draws
  wsample <- wsample[,sample(1:NCOL(wsample)),drop=F]
  wsample <- wsample[,1:draws]
  
  zt <- gp$featuremap(xt)
  sample <- zt %*% wsample
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_linearized_analytic <- function(gp, xt, draws=NULL, transform=T,
                                        cfind=NULL, jitter=NULL) {
  
  # TODO: take into account cfind
  if (!is.null(cfind))
    stop('cfind not supported for method rf yet.')

  # draw from the posterior of w
  zt <- gp$featuremap(xt)
  w <- mvnrnd(draws, gp$wmean, chol_prec=gp$wprec_chol)
  sample <- zt %*% w
  if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_linearized_prior <- function(gp, xt, var=F, draws=NULL, transform=T,
                                     cfind=NULL, jitter=NULL) {
  
  # TODO: take into account cfind
  if (!is.null(cfind))
    stop('cfind not supported for method rf yet.')

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
    r <- mean + chol_cov %*% matrix(stats::rnorm(d*draws), nrow=d)
  else if (!is.null(chol_prec))
    r <- mean + solve(t(chol_prec), matrix(stats::rnorm(d*draws), nrow=d))
  else
    stop('Both cov_chol and prec_chol can\'t be NULL.')
  return(r)
}



