




#' @rdname pred
#' @export
gp_draw <- function(gp, xnew, draws=NULL, transform=T, target=F, marginal=F,
                    cfind=NULL, jitter=NULL, seed=NULL, ...) {
  
  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(seed)

  if (is_fitted(gp, 'sampling')) {
    # model fitted using mcmc, so predict using the draws from the posterior
    pred <- gp_draw_mcmc(gp, xnew, draws=draws, transform=transform,
                         target=target, marginal=marginal, cfind=cfind, 
                         jitter=jitter, ...)
  } else {
    # model fitted using analytical gaussian approximation (or not fitted at all),
    # so predict based on that
    if (is.null(draws))
      draws <- 1
    if (is_fitted(gp, 'analytic')) {
      # draw from the analytical posterior approximation
      pred <- gp_draw_analytic(gp, xnew, draws=draws, transform=transform,
                               target=target, marginal=marginal, cfind=cfind, 
                               jitter=jitter, ...)
    } else {
      # draw from the prior
      pred <- gp_draw_prior(gp, xnew, draws=draws, transform=transform,
                            target=target, cfind=cfind, jitter=jitter, ...)
    }
  }
  return(pred)
}


gp_draw_prior <- function(object, ...) {
  UseMethod('gp_draw_prior', object)
}

gp_draw_analytic <- function(object, ...) {
  UseMethod('gp_draw_analytic', object)
}

gp_draw_mcmc <- function(object, ...) {
  UseMethod('gp_draw_mcmc', object)
}

gp_draw_prior.gp <- function(object, ...) {
  gp_draw_prior(object$approx, object, ...)
}

gp_draw_analytic.gp <- function(object, ...) {
  gp_draw_analytic(object$approx, object, ...)
}

gp_draw_mcmc.gp <- function(object, ...) {
  gp_draw_mcmc(object$approx, object, ...)
}


gp_draw_prior.approx_full <- function(object, gp, xt, draws=NULL, transform=T, target=F,
                                      cfind=NULL, jitter=NULL, ...) {
  
  pred <- gp_pred_prior(object, gp, xt, cov=T, cfind=cfind, jitter=jitter)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_prior.approx_fitc <- function(object, gp, xt, draws=NULL, transform=T, target=F,
                                      cfind=NULL, jitter=NULL, ...) {
  
  pred <- gp_pred_prior(object, gp, xt, cov=T, cfind=cfind, jitter=jitter)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_prior.approx_rf <- function(object, gp, xt, var=F, draws=NULL, transform=T, target=F,
                                    cfind=NULL, ...) {
  
  # draw from the prior of w
  num_inputs <- NCOL(xt)
  featuremap <- get_featuremap(gp, num_inputs)
  zt <- featuremap(xt, cfind)
  num_feat <- NCOL(zt)
  w <- matrix(stats::rnorm(num_feat*draws), nrow=num_feat) # draw from the prior (standard Gaussian)
  sample <- zt %*% w
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}


gp_draw_analytic.approx_full <- function(object, gp, xt, draws=NULL, transform=T, target=F,
                                         marginal=F, cfind=NULL, jitter=NULL, ...) {
  if (marginal) {
    pred <- gp_pred_post(object, gp, xt, cov=F, var=T, cfind=cfind, jitter=jitter)
    sample <- mvnrnd(draws, pred$mean, chol_cov = sqrt(pred$var))
  } else {
    pred <- gp_pred_post(object, gp, xt, cov=T, cfind=cfind, jitter=jitter)
    sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  }
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_analytic.approx_fitc <- function(object, gp, xt, draws=NULL, transform=T, target=F,
                                         marginal=F, cfind=NULL, jitter=NULL, ...) {
  if (marginal) {
    pred <- gp_pred_post(object, gp, xt, cov=F, var=T, cfind=cfind, jitter=jitter)
    sample <- mvnrnd(draws, pred$mean, chol_cov = sqrt(pred$var))
  } else {
    pred <- gp_pred_post(object, gp, xt, cov=T, cfind=cfind, jitter=jitter)
    sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  }
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_analytic.approx_rf <- function(object, gp, xt, draws=NULL, transform=T, target=F,
                                       cfind=NULL, ...) {
  
  # draw from the posterior of w
  num_inputs <- NCOL(xt)
  featuremap <- get_featuremap(gp, num_inputs)
  zt <- featuremap(xt, cfind)
  wmean <- get_w_mean(gp, cfind)
  wcov <- get_w_cov(gp, cfind)
  w <- mvnrnd(draws, wmean, chol_cov = t(chol(wcov)))
  sample <- zt %*% w
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_analytic.approx_rbf <- function(object, gp, xt, draws=NULL, transform=T, target=F,
                                        cfind=NULL, ...) {
  gp_draw_analytic.approx_rf(object, gp, xt, draws=draws, transform=transform, target=target,
                             cfind=cfind, ...)
}




gp_draw_mcmc.approx_full <- function(object, gp, xt, draws=NULL, transform=T, target=T,
                                     marginal=F, cfind=NULL, jitter=NULL, ...) {
  
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
  if (marginal)
    sample <- mvnrnd(draws, pred_mean[,1:draws], chol_cov = sqrt(diag(pred_cov)))
  else
    sample <- mvnrnd(draws, pred_mean[,1:draws], chol_cov = t(chol(pred_cov)))
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}

gp_draw_mcmc.approx_rf <- function(object, gp, xt, draws=NULL, transform=T, target=T, 
                                   cfind=NULL, ...) {
  
  wsample <- get_w_sample(gp, cfind)
  if (is.null(draws))
    draws <- NCOL(wsample)
  else if (draws > NCOL(wsample))
    stop('Can\'t draw more than the provided GP has posterior draws.')
  # permute the posterior draws and pick right number of draws
  wsample <- wsample[,sample(1:NCOL(wsample)),drop=F]
  wsample <- wsample[,1:draws]
  
  num_inputs <- NCOL(xt)
  featuremap <- get_featuremap(gp, num_inputs)
  zt <- featuremap(xt, cfind)
  sample <- zt %*% wsample
  if (target)
    sample <- generate_target(gp, sample, ...)
  else if (transform)
    sample <- get_response(gp, sample)
  return(sample)
}


mvnrnd <- function(draws, mean, chol_cov=NULL, chol_prec=NULL) {
  # Draw from multivariate normal. Lower Cholesky factor of either
  # covariance or precision must be provided.
  d <- NROW(mean)
  if (!is.null(chol_cov)) {
    if (is.vector(chol_cov)) {
      # diagonal covariance matrix
      r <- mean + chol_cov * matrix(stats::rnorm(d*draws), nrow=d)
    } else {
      # non-diagonal covariance matrix
      r <- mean + chol_cov %*% matrix(stats::rnorm(d*draws), nrow=d)
    }
  } else if (!is.null(chol_prec)) {
    if (is.vector(chol_prec)) {
      # diagonal precision
      chol_cov <- 1/chol_prec
      r <- mean + chol_cov * matrix(stats::rnorm(d*draws), nrow=d)
    } else {
      # non-diagonal covariance matrix
      r <- mean + solve(t(chol_prec), matrix(stats::rnorm(d*draws), nrow=d))
    }
  } else
    stop('Both cov_chol and prec_chol can\'t be NULL.')
  return(r)
}



