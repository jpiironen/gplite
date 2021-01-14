




#' @rdname pred
#' @export
gp_draw <- function(gp, xnew, draws = NULL, transform = TRUE, target = FALSE, marginal = FALSE,
                    cfind = NULL, jitter = NULL, seed = NULL, ...) {

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)


  # model fitted using analytical gaussian approximation (or not fitted at all),
  # so predict based on that
  if (is.null(draws)) {
    draws <- 1
  }
  if (is_fitted(gp, "analytic")) {
    # draw from the analytical posterior approximation
    pred <- gp_draw_analytic(gp, xnew,
      draws = draws, transform = transform,
      target = target, marginal = marginal, cfind = cfind,
      jitter = jitter, ...
    )
  } else {
    # draw from the prior
    pred <- gp_draw_prior(gp, xnew,
      draws = draws, transform = transform,
      target = target, cfind = cfind, jitter = jitter, ...
    )
  }

  return(pred)
}


gp_draw_prior <- function(object, ...) {
  UseMethod("gp_draw_prior", object)
}

gp_draw_analytic <- function(object, ...) {
  UseMethod("gp_draw_analytic", object)
}

gp_draw_prior.gp <- function(object, ...) {
  gp_draw_prior(object$method, object, ...)
}

gp_draw_analytic.gp <- function(object, ...) {
  gp_draw_analytic(object$method, object, ...)
}



gp_draw_prior.method_full <- function(object, gp, xt, draws = NULL, transform = TRUE, target = FALSE,
                                      cfind = NULL, jitter = NULL, ...) {
  pred <- gp_pred_prior(object, gp, xt, cov = TRUE, cfind = cfind, jitter = jitter)
  pred$mean <- add_offset(pred$mean, ...)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (target) {
    sample <- generate_target(gp, sample, ...)
  } else if (transform) {
    sample <- get_response(gp, sample)
  }
  return(sample)
}

gp_draw_prior.method_fitc <- function(object, gp, xt, draws = NULL, transform = TRUE, target = FALSE,
                                      cfind = NULL, jitter = NULL, ...) {
  pred <- gp_pred_prior(object, gp, xt, cov = TRUE, cfind = cfind, jitter = jitter)
  pred$mean <- add_offset(pred$mean, ...)
  sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  if (target) {
    sample <- generate_target(gp, sample, ...)
  } else if (transform) {
    sample <- get_response(gp, sample)
  }
  return(sample)
}

gp_draw_prior.method_rf <- function(object, gp, xt, var = FALSE, draws = NULL, transform = TRUE, target = FALSE,
                                    cfind = NULL, ...) {

  # draw from the prior of w
  num_inputs <- NCOL(xt)
  featuremap <- get_featuremap(gp, num_inputs)
  zt <- featuremap(xt, cfind)
  num_feat <- NCOL(zt)
  w <- matrix(stats::rnorm(num_feat * draws), nrow = num_feat) # draw from the prior (standard Gaussian)
  sample <- zt %*% w
  sample <- add_offset(sample, ...)
  if (target) {
    sample <- generate_target(gp, sample, ...)
  } else if (transform) {
    sample <- get_response(gp, sample)
  }
  return(sample)
}


gp_draw_analytic.method_full <- function(object, gp, xt, draws = NULL, transform = TRUE, target = FALSE,
                                         marginal = FALSE, cfind = NULL, jitter = NULL, ...) {
  if (marginal) {
    pred <- gp_pred_post(object, gp, xt, cov = FALSE, var = TRUE, cfind = cfind, jitter = jitter)
    pred$mean <- add_offset(pred$mean, ...)
    sample <- mvnrnd(draws, pred$mean, chol_cov = sqrt(pred$var))
  } else {
    pred <- gp_pred_post(object, gp, xt, cov = TRUE, cfind = cfind, jitter = jitter)
    pred$mean <- add_offset(pred$mean, ...)
    sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  }
  if (target) {
    sample <- generate_target(gp, sample, ...)
  } else if (transform) {
    sample <- get_response(gp, sample)
  }
  return(sample)
}

gp_draw_analytic.method_fitc <- function(object, gp, xt, draws = NULL, transform = TRUE, target = FALSE,
                                         marginal = FALSE, cfind = NULL, jitter = NULL, ...) {
  if (marginal) {
    pred <- gp_pred_post(object, gp, xt, cov = FALSE, var = TRUE, cfind = cfind, jitter = jitter)
    pred$mean <- add_offset(pred$mean, ...)
    sample <- mvnrnd(draws, pred$mean, chol_cov = sqrt(pred$var))
  } else {
    pred <- gp_pred_post(object, gp, xt, cov = TRUE, cfind = cfind, jitter = jitter)
    pred$mean <- add_offset(pred$mean, ...)
    sample <- mvnrnd(draws, pred$mean, chol_cov = t(chol(pred$cov)))
  }
  if (target) {
    sample <- generate_target(gp, sample, ...)
  } else if (transform) {
    sample <- get_response(gp, sample)
  }
  return(sample)
}

gp_draw_analytic.method_rf <- function(object, gp, xt, draws = NULL, transform = TRUE, target = FALSE,
                                       cfind = NULL, ...) {

  # draw from the posterior of w
  num_inputs <- NCOL(xt)
  featuremap <- get_featuremap(gp, num_inputs)
  zt <- featuremap(xt, cfind)
  wmean <- get_w_mean(gp, cfind)
  wcov <- get_w_cov(gp, cfind)
  w <- mvnrnd(draws, wmean, chol_cov = t(chol(wcov)))
  sample <- zt %*% w
  sample <- add_offset(sample, ...)
  if (target) {
    sample <- generate_target(gp, sample, ...)
  } else if (transform) {
    sample <- get_response(gp, sample)
  }
  return(sample)
}

gp_draw_analytic.method_rbf <- function(object, gp, xt, draws = NULL, transform = TRUE, target = FALSE,
                                        cfind = NULL, ...) {
  gp_draw_analytic.method_rf(object, gp, xt,
    draws = draws, transform = transform, target = target,
    cfind = cfind, ...
  )
}










mvnrnd <- function(draws, mean, chol_cov = NULL, chol_prec = NULL) {
  # Draw from multivariate normal. Lower Cholesky factor of either
  # covariance or precision must be provided.
  d <- NROW(mean)
  if (!is.null(chol_cov)) {
    if (is.vector(chol_cov)) {
      # diagonal covariance matrix
      r <- mean + chol_cov * matrix(stats::rnorm(d * draws), nrow = d)
    } else {
      # non-diagonal covariance matrix
      r <- mean + chol_cov %*% matrix(stats::rnorm(d * draws), nrow = d)
    }
  } else if (!is.null(chol_prec)) {
    if (is.vector(chol_prec)) {
      # diagonal precision
      chol_cov <- 1 / chol_prec
      r <- mean + chol_cov * matrix(stats::rnorm(d * draws), nrow = d)
    } else {
      # non-diagonal covariance matrix
      r <- mean + solve(t(chol_prec), matrix(stats::rnorm(d * draws), nrow = d))
    }
  } else {
    stop("Both cov_chol and prec_chol can't be NULL.")
  }
  return(r)
}
