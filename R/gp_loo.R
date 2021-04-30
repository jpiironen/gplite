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
#' @param quadrature Whether to use deterministic Gauss-Hermite quadrature to estimate the
#'  required integrals. If FALSE, then Monte Carlo estimate is used.
#' @param quad_order Order of the numerical quadrature
#'  (only applicable if \code{quadrature=TRUE}).
#' @param draws Number of posterior draws to estimate the required integrals (only applicable
#'  if \code{quadrature=FALSE}).
#' @param jitter Magnitude of diagonal jitter for covariance matrices for numerical stability.
#'  Default is 1e-6.
#' @param seed Random seed.
#' @param ref Index of the model against which to compare the other models (pairwise
#' comparison for LOO difference). If not given, then the model with the best LOO is
#' used as the reference for comparisons.
#' @param ... For \code{gp_compare}, LOO statistics for the models to compare. For
#' \code{gp_loo}, possible additional data that is required for LOO predictions (for example,
#' argument \code{trials} in case of binomial likelihood).
#'
#'
#' @return \code{gp_loo} returns a list with LOO statistics.
#' \code{gp_compare} returns a matrix with comparison statistics (LOO differences 
#' and stardard errors in the estimates).
#'
#' @section References:
#'
#' Vehtari A., Mononen T., Tolvanen V., Sivula T. and Winther O. (2016).
#' Bayesian Leave-One-Out Cross-Validation Approximations for Gaussian Latent
#' Variable Models. Journal of Machine Learning Research 17(103):1-38.
#'
#' @examples
#'
#' # Generate some toy data
#' set.seed(32004)
#' n <- 50
#' sigma <- 0.1
#' x <- rnorm(n)
#' ycont <- sin(3 * x) * exp(-abs(x)) + rnorm(n) * sigma
#' y <- rep(0, n)
#' y[ycont > 0] <- 1
#' trials <- rep(1, n)
#'
#' # Set up two models
#' gp1 <- gp_init(cf_sexp(), lik_binomial())
#' gp2 <- gp_init(cf_matern32(), lik_binomial())
#'
#' # Optimize
#' gp1 <- gp_optim(gp1, x, y, trials = trials)
#' gp2 <- gp_optim(gp2, x, y, trials = trials)
#'
#' # Compare
#' loo1 <- gp_loo(gp1, x, y, trials = trials)
#' loo2 <- gp_loo(gp2, x, y, trials = trials)
#' gp_compare(loo1, loo2)
#' 
#'
NULL


#' @rdname gp_loo
#' @export
gp_loo <- function(gp, x, y, quadrature = TRUE, quad_order = 11, draws = 4000,
                   jitter = NULL, seed = NULL, ...) {
  if (!is_fitted(gp, "analytic")) {
    stop("The provied GP model is not fitted. Please call gp_fit (or gp_optim) first.")
  }

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  loo_post <- loo_posteriors(gp, x, y, ...)
  loo_mean <- loo_post$mean
  loo_var <- loo_post$var

  if (quadrature) {
    # use Gauss-Hermite quadrature to evaluate predictive distribution
    quadrature <- gauss_hermite_points_scaled(loo_mean, sqrt(loo_var), order = quad_order)
    fgrid <- quadrature$x
    weights <- quadrature$weights
    loglik <- get_loglik(gp$lik, fgrid, y, sum = FALSE, ...)
    loos <- apply(loglik, 1, logsumexp, weights = weights)
  } else {
    # sample from LOO posteriors and evaluate predictive distribution using Monte Carlo
    n <- length(y)
    fsample <- matrix(stats::rnorm(n * draws, mean = loo_mean, sd = sqrt(loo_var)), nrow = n)
    loglik <- get_loglik(gp$lik, fsample, y, sum = FALSE, ...)
    loos <- apply(loglik, 1, logsumexp) - log(draws)
  }

  n <- length(y)
  res <- list(loo = sum(loos), sd = n * (stats::sd(loos) / sqrt(n)), loos = loos)
  class(res) <- "loores"
  return(res)
}

loo_posteriors <- function(object, ...) {
  UseMethod("loo_posteriors", object)
}

loo_posteriors.gp <- function(object, ...) {
  loo_posteriors(object$approx, object, ...)
}

loo_posteriors.approx_laplace <- function(object, gp, x, y, ...) {
  fhat <- as.vector(gp$fit$fmean)
  z <- gp$fit$z
  V <- gp$fit$V
  grad <- (z - fhat) / V
  post_pred <- gp_pred_post(gp$method, gp, x, var = TRUE, train = TRUE)
  post_var <- post_pred$var
  loo_var <- 1 / (1 / post_var - 1 / V)
  loo_mean <- fhat - loo_var * grad
  return(list(mean = loo_mean, var = loo_var))
}

loo_posteriors.approx_ep <- function(object, gp, x, y, ...) {
  return(list(mean = gp$fit$cavity_mean, var = gp$fit$cavity_var))
}


logsumexp <- function(x, weights = NULL) {
  # numerically stable computation of log(sum_i(w_i*exp(x_i))),
  # where w_i are optional weights (will be ones if not given)
  c <- max(x)
  if (is.null(weights)) {
    weights <- rep(1, length(x))
  }
  c + log(sum(weights * exp(x - c)))
}


#' @rdname gp_loo
#' @export
gp_compare <- function(..., ref=NULL) {
  args <- list(...)
  loo <- c()
  loos <- list()
  for (i in seq_along(args)) {
    if (!("loores" %in% class(args[[i]]))) {
      stop("Expected objects of type 'loores', but found an object with type '", class(args[[i]]), "'")
    }

    loo[i] <- args[[i]]$loo
    loos[[i]] <- args[[i]]$loos
  }

  if (is.null(ref))
    ref <- which.max(loo)

  res <- matrix(nrow = length(loos), ncol = 2)
  n <- length(loos[[1]])
  for (i in seq_along(loos)) {
    d <- loos[[i]] - loos[[ref]]
    res[i, 1] <- n * mean(d)
    res[i, 2] <- n * (stats::sd(d) / sqrt(n))
  }
  rownames(res) <- sapply(1:length(loos), function(i) sprintf("model%i", i))
  colnames(res) <- c("loo-diff", "se")
  return(res)
}
