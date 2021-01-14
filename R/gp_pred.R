#' Make predictions with a GP model
#'
#' Function \code{gp_pred} can be used to make analytic predictions for the latent function
#' values at test points, whereas \code{gp_draw}
#' can be used to draw from the predictive distribution (or from the prior if the GP has
#' not been fitted yet.)
#'
#' @name pred
#'
#' @param gp A GP model object.
#' @param xnew N-by-d matrix of input values (N is the number of test points and d
#' the input dimension).
#' Can also be a vector of length N if the model has only a single input.
#' @param var Whether to compute the predictive variances along with predictive mean.
#' @param quantiles Vector of probabilities between 0 and 1 indicating which quantiles are to
#' be predicted.
#' @param draws Number of draws to generate from the predictive distribution for the
#' latent values.
#' @param transform Whether to transform the draws of latent values to the same scale
#'  as the target y, that is, through the response (or inverse-link) function.
#' @param target If TRUE, draws values for the target variable \code{y} instead of the latent
#'  function values.
#' @param marginal If TRUE, then draws for each test point are only marginally correct, but the
#'  covariance structure between test points is not retained. However, this will make the sampling
#'  considerably faster in some cases, and can be useful if one is interested only in looking
#'  at the marginal predictive distributions for a large number of test locations
#'  (for example, in posterior predictive checking).
#' @param cfind Indices of covariance functions to be used in the prediction. By default uses
#' all covariance functions.
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability.
#'  Default is 1e-6.
#' @param quad_order Quadrature order in order to compute the mean and variance on
#' the transformed scale.
#' @param seed Random seed for draws.
#' @param ... Additional parameters that might be needed. For example \code{offset} or 
#' keyword \code{trials} for binomial and beta-binomial likelihoods.
#'
#'
#' @return \code{gp_pred} returns a list with fields giving the predictive mean, variance and
#' quantiles (the last two are computed only if requested). \code{gp_draw} returns an N-by-draws
#' matrix of random draws from the predictive distribution, where N is the number of test points.
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
#' # More than one covariance function; one for x1 and x2, and another one for x3
#' cf1 <- cf_nn(c("x1", "x2"), prior_sigma0 = prior_half_t(df = 4, scale = 2))
#' cf2 <- cf_sexp("x3")
#' cfs <- list(cf1, cf2)
#' lik <- lik_gaussian()
#' gp <- gp_init(cfs, lik)
#' gp <- gp_optim(gp, x, y, maxiter = 500)
#'
#' # plot the predictions with respect to x1, when x2 = x3 = 0
#' xt <- cbind(x1 = seq(-3, 3, len = 50), x2 = 0, x3 = 0)
#' pred <- gp_pred(gp, xt)
#' plot(xt[, "x1"], pred$mean, type = "l")
#'
#' # draw from the predictive distribution
#' xt <- cbind(x1 = seq(-3, 3, len = 50), x2 = 0, x3 = 0)
#' draws <- gp_draw(gp, xt, draws = 100)
#' plot(xt[, "x1"], draws[, 1], type = "l")
#' for (i in 2:50) {
#'   lines(xt[, "x1"], draws[, i])
#' }
#'
#' # plot effect with respect to x3 only
#' xt <- cbind("x3" = seq(-3, 3, len = 50))
#' pred <- gp_pred(gp, xt, cfind = 2)
#' plot(xt, pred$mean, type = "l")
#' 
#'
NULL

#' @rdname pred
#' @export
gp_pred <- function(gp, xnew, var = FALSE, quantiles = NULL, transform = FALSE, cfind = NULL, 
                    jitter = NULL, quad_order = 15, ...) {
  if (!is.null(quantiles) || transform) {
    # we need variances in order to compute the quantiles, or to transform the mean
    var <- TRUE
  }

  if (!is_fitted(gp, "analytic")) {
    # model not fitted, so predict based on the prior
    pred <- gp_pred_prior(gp, xnew, var = var, cfind = cfind, jitter = jitter)
  } else {
    # model fitted using analytical approximation
    pred <- gp_pred_post(gp, xnew, var = var, cfind = cfind, jitter = jitter)
  }
  pred$mean <- add_offset(pred$mean, ...)
  
  if (!is.null(quantiles)) {
    quantiles <- sapply(quantiles, function(q) stats::qnorm(q, mean = pred$mean, sd = sqrt(pred$var)))
    if (transform) {
      quantiles <- get_response(gp, quantiles)
    }
    pred$quantiles <- quantiles
  }
  if (transform) {
    quadrature <- gauss_hermite_points_scaled(pred$mean, sqrt(pred$var), order = quad_order)
    fgrid <- quadrature$x
    weights <- quadrature$weights
    fgrid_transf <- get_response(gp, fgrid)
    mean_transf <- as.vector(fgrid_transf %*% weights)
    var_transf <- as.vector((fgrid_transf - mean_transf)^2 %*% weights)
    pred$mean <- mean_transf
    pred$var <- var_transf
  }
  return(pred)
}


gp_pred_prior <- function(object, ...) {
  UseMethod("gp_pred_prior", object)
}

gp_pred_post <- function(object, ...) {
  UseMethod("gp_pred_post", object)
}

gp_pred_prior.gp <- function(object, ...) {
  gp_pred_prior(object$method, object, ...)
}

gp_pred_post.gp <- function(object, ...) {
  gp_pred_post(object$method, object, ...)
}


gp_pred_prior.method_full <- function(object, gp, xt, var = FALSE, cov = FALSE, cfind = NULL, jitter = NULL) {
  nt <- NROW(xt)
  pred_mean <- rep(0, nt)

  if (var || cov) {
    jitter <- get_jitter(gp, jitter)
    pred_cov <- eval_cf(gp$cfs, xt, xt, cfind)
    if (cov) {
      return(list(mean = pred_mean, cov = pred_cov + jitter * diag(nt)))
    } else {
      return(list(mean = pred_mean, var = diag(pred_cov)))
    }
  }
  return(list(mean = pred_mean))
}

gp_pred_prior.method_fitc <- function(object, gp, xt, var = FALSE, cov = FALSE, cfind = NULL, jitter = NULL) {
  nt <- NROW(xt)
  pred_mean <- rep(0, nt)

  if (var || cov) {
    jitter <- get_jitter(gp, jitter)
    # TODO: the code below is repeating a lot what is already written elsewhere
    if (is.null(gp$method$inducing)) {
      stop("Inducing points not set yet.")
    }
    z <- gp$method$inducing
    Kz <- eval_cf(gp$cfs, z, z) + jitter * diag(gp$method$num_inducing)
    Kxz <- eval_cf(gp$cfs, xt, z)
    Kz_chol <- t(chol(Kz))
    xt <- as.matrix(xt)
    D <- eval_cf(gp$cfs, xt, xt, diag_only = TRUE)
    D <- D - colSums(forwardsolve(Kz_chol, t(Kxz))^2)
    aux <- forwardsolve(Kz_chol, t(Kxz))
    pred_cov <- t(aux) %*% aux + diag(D)
    if (cov) {
      return(list(mean = pred_mean, cov = pred_cov + jitter * diag(nt)))
    } else {
      return(list(mean = pred_mean, var = diag(pred_cov)))
    }
  }
  return(list(mean = pred_mean))
}

gp_pred_prior.method_rf <- function(object, gp, xt, var = FALSE, cfind = NULL, jitter = NULL) {

  # mean is zero
  nt <- NROW(xt)
  pred_mean <- rep(0, nt)

  if (var == TRUE) {
    num_inputs <- NCOL(xt)
    featuremap <- get_featuremap(gp, num_inputs)
    zt <- featuremap(xt, cfind)
    pred_cov <- zt %*% t(zt)
    return(list(mean = pred_mean, var = diag(pred_cov)))
  }
  return(list(mean = pred_mean))
}

gp_pred_post.method_full <- function(object, gp, xt, var = FALSE, cov = FALSE, train = FALSE,
                                     cfind = NULL, jitter = NULL) {

  # compute the latent mean first
  Kt <- eval_cf(gp$cfs, xt, gp$x, cfind)
  Ktt <- eval_cf(gp$cfs, xt, xt, cfind)
  alpha <- gp$fit$alpha
  pred_mean <- Kt %*% alpha
  pred_mean <- as.vector(pred_mean)

  if (var || cov) {
    # (co)variance of the latent function
    nt <- length(pred_mean)
    jitter <- get_jitter(gp, jitter)
    if (!is.null(gp$fit$C_chol)) {
      C_chol <- gp$fit$C_chol
      aux <- forwardsolve(C_chol, t(Kt))
      pred_cov <- Ktt - t(aux) %*% aux
    } else {
      C_lu <- gp$fit$C_lu
      var_reduction <- Kt %*% solve(C_lu$U, solve(C_lu$L, solve(C_lu$P, t(Kt))))
      pred_cov <- Ktt - var_reduction
    }
    if (cov) {
      return(list(mean = pred_mean, cov = pred_cov + jitter * diag(nt)))
    } else {
      return(list(mean = pred_mean, var = diag(pred_cov)))
    }
  }
  return(list(mean = pred_mean))
}

gp_pred_post.method_fitc <- function(object, gp, xt, var = FALSE, cov = FALSE, train = FALSE,
                                     cfind = NULL, jitter = NULL) {

  # compute the latent mean first
  z <- gp$method$inducing
  alpha <- gp$fit$alpha
  Ktz <- eval_cf(gp$cfs, xt, z, cfind)
  Kxz <- gp$fit$Kxz
  Kz_chol <- gp$fit$Kz_chol
  if (var || cov || train) {
    # FITC diagonal correction
    Dt <- eval_cf(gp$cfs, xt, xt, cfind, diag_only = TRUE)
    Dt <- Dt - colSums(forwardsolve(Kz_chol, t(Ktz))^2)
  }
  # with these auxiliary matrices, we have: Kt == t(Kt_aux2) %*% Kt_aux1
  Kt_aux1 <- forwardsolve(Kz_chol, t(Kxz))
  Kt_aux2 <- forwardsolve(Kz_chol, t(Ktz))

  pred_mean <- as.vector(t(Kt_aux2) %*% (Kt_aux1 %*% alpha))
  if (train) {
    # need to take into account the diagonal correction, since training prediction
    pred_mean <- pred_mean + Dt * alpha
  }

  if (var || cov) {
    # (co)variance of the latent function
    nt <- length(pred_mean)
    jitter <- get_jitter(gp, jitter)
    Kz <- gp$fit$Kz
    xt <- as.matrix(xt)

    if (cov) {
      Linv_Ktz <- forwardsolve(Kz_chol, t(Ktz))
      prior_cov <- t(Linv_Ktz) %*% Linv_Ktz + diag(Dt)
      C_inv <- gp$fit$C_inv
      if (train) {
        Kt <- t(Kt_aux2) %*% Kt_aux1 + diag(Dt)
        pred_cov <- prior_cov - Kt %*% inv_lemma_solve(C_inv, t(Kt))
      } else {
        pred_cov <- prior_cov - t(Kt_aux2) %*% (Kt_aux1 %*% inv_lemma_solve(C_inv, t(Kt_aux1))) %*% Kt_aux2
      }
      return(list(mean = pred_mean, cov = pred_cov + jitter * diag(nt)))
    } else {
      C_inv <- gp$fit$C_inv
      Kz_inv_Ktz <- backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Ktz)))
      prior_var <- rowSums(Ktz * t(Kz_inv_Ktz)) + Dt
      if (train) {
        aux <- inv_lemma_solve(C_inv, t(Kt_aux1), Kt_aux1) %*% Kt_aux2 +
          inv_lemma_solve(C_inv, Dt, Kt_aux1, rhs_diag = TRUE)
        diag1 <- colSums(Kt_aux2 * aux)
        diag2 <- rowSums(inv_lemma_solve(C_inv, t(Kt_aux1), Dt, lhs_diag = TRUE) * t(Kt_aux2))
        diag3 <- inv_lemma_solve(C_inv, Dt, Dt, rhs_diag = TRUE, lhs_diag = TRUE, diag_only = TRUE)
        var_reduction <- diag1 + diag2 + diag3
      } else {
        aux <- inv_lemma_solve(C_inv, t(Kt_aux1), Kt_aux1) %*% Kt_aux2
        var_reduction <- colSums(Kt_aux2 * aux)
      }
      pred_var <- prior_var - var_reduction
      return(list(mean = pred_mean, var = as.vector(pred_var)))
    }
  }
  return(list(mean = pred_mean))
}

gp_pred_post.method_rf <- function(object, gp, xt, var = FALSE, train = FALSE, cfind = NULL, jitter = NULL) {

  # compute the latent mean first
  featuremap <- get_featuremap(gp, num_inputs = NCOL(xt))
  zt <- featuremap(xt, cfind)
  wmean <- get_w_mean(gp, cfind)
  pred_mean <- as.vector(zt %*% wmean)

  if (var == TRUE) {
    # covariance of the latent function
    nt <- length(pred_mean)
    wcov <- get_w_cov(gp, cfind)
    pred_cov <- zt %*% (wcov %*% t(zt))
    return(list(mean = pred_mean, var = diag(pred_cov)))
  }
  return(list(mean = pred_mean))
}

gp_pred_post.method_rbf <- function(object, gp, xt, var = FALSE, train = FALSE, cfind = NULL, jitter = NULL) {
  gp_pred_post.method_rf(object, gp, xt, var = var, train = train, cfind = cfind, jitter = jitter)
}
