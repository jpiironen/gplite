

laplace_iter <- function(object, ...) {
  UseMethod("laplace_iter", object)
}

laplace <- function(object, ...) {
  UseMethod("laplace", object)
}



laplace_iter.method_full <- function(object, gp, K, y, fhat_old, pobs = NULL, ...) {

  # get the pseudo data first
  if (is.null(pobs)) {
    pobs <- get_pseudodata_la(gp$lik, fhat_old, y, ...)
  }
  z <- pobs$z
  V <- pobs$var
  n <- length(z)

  if (all(V > 0)) {

    # normal case; all pseudo-variances are positive (log-concave likelihood)

    B <- diag(n) + diag_times_dense(sqrt(1 / V), diag_times_dense(K, sqrt(1 / V), diag_right = TRUE))
    C_chol <- diag_times_dense(sqrt(V), t(chol(B)))
    alpha <- backsolve(t(C_chol), forwardsolve(C_chol, z))
    fhat_new <- K %*% alpha

    # compute the log marginal likelihood
    C_logdet <- 2 * sum(log(diag(C_chol)))
    V_logdet <- sum(log(V))
    f_invK_f <- t(fhat_new) %*% alpha
    log_lik <- get_loglik(gp$lik, fhat_new, y, ...)
    log_evidence <- log_lik - 0.5 * f_invK_f - 0.5 * (C_logdet - V_logdet)
    log_evidence <- as.numeric(log_evidence)
    
  } else {

    # non log-concave likelihood; needs special treatment.
    # evaluate the marginal likelihood as p(y) = p(y1,y2) = p(y1)*p(y2 | y1)
    # where y1 are the good data (for which V > 0), and y2 are the bad ones
    good <- which(V > 0)
    bad <- which(V <= 0)
    nbad <- length(bad)

    if (nbad > 500) {
      warning("More than 500 data points with pseudo-variance < 0 (outliers), computation getting slow. Check that the model and hyperparameter values are sensible.")
    }

    args <- c(
      list(
        object = object, gp = gp, K = K[good, good, drop = FALSE], y = y[good], fhat_old = NULL,
        pobs = list(z = z[good], var = V[good])
      ),
      list(...)
    )
    for (argname in required_extra_args(gp$lik)) {
      args[[argname]] <- args[[argname]][good]
    }
    fit1 <- do.call(laplace_iter, args)
    log_evidence1 <- fit1$log_evidence

    # compute predictive mean and covariance for f2, given posterior for y1, and
    # then use these as 'prior' mean and covariance to compute marginal likelihood
    # of y2
    K21 <- K[bad, good, drop = FALSE]
    K22 <- K[bad, bad, drop = FALSE]
    pred_mean <- as.vector(K21 %*% fit1$alpha)
    aux <- forwardsolve(fit1$C_chol, t(K21))
    pred_cov <- K22 - t(aux) %*% aux

    # if the following Cholesky computations fail, then the posterior of the whole
    # data is not proper (the covariance is not positive definite)
    Lk <- t(chol(pred_cov))
    A <- diag(nbad) + t(Lk) %*% diag_times_dense(1 / V[bad], Lk)
    A_chol <- tryCatch(
      {
        t(chol(A))
      },
      error = function(err) {
        NULL
      }
    )
    if (is.null(A_chol)) {
      return(list(fmean = fhat_old, log_evidence = -Inf))
    }
    alpha2 <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% (z[bad] / V[bad]))))

    aux <- backsolve(t(Lk), forwardsolve(Lk, pred_mean))
    aux <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% aux)))
    f2hat <- pred_cov %*% (alpha2 + aux)

    # logdet((inv(pred_cov + V)*V) = logdet((inv(inv(V)*pred_cov + I))
    invC2_V2_logdet <- 2 * sum(log(diag(A_chol)))

    df2_invK2_df2 <- sum(forwardsolve(Lk, f2hat - pred_mean)^2)
    args <- c(list(object = gp$lik, f = f2hat, y = y[bad]), list(...))
    for (argname in required_extra_args(gp$lik)) {
      args[[argname]] <- args[[argname]][bad]
    }
    log_lik2 <- do.call(get_loglik, args)
    log_evidence2 <- log_lik2 - 0.5 * df2_invK2_df2 + 0.5 * invC2_V2_logdet
    log_evidence2 <- as.numeric(log_evidence2)

    log_evidence <- log_evidence1 + log_evidence2

    # finally, fit the model using all the observations.
    # here we use LU-decomposition instead of Cholesky for matrix K+V,
    # as it might not be positive definite
    C_lu <- Matrix::expand(Matrix::lu(K + diag(V)))
    alpha <- solve(C_lu$U, solve(C_lu$L, solve(C_lu$P, z)))
    fhat_new <- K %*% alpha
    return(list(fmean = fhat_new, z = z, V = V, C_lu = C_lu, alpha = alpha, log_evidence = log_evidence))
  }

  list(fmean = fhat_new, z = z, V = V, C_chol = C_chol, alpha = alpha, log_evidence = log_evidence)
}

laplace_iter.method_fitc <- function(object, gp, Kz, Kz_chol, Kxz, D, y, fhat_old,
                                     pobs = NULL, ...) {

  # get the pseudo data first
  if (is.null(pobs)) {
    pobs <- get_pseudodata_la(gp$lik, fhat_old, y, ...)
  }
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  S <- D + V


  if (all(V > 0)) {

    # normal case; all pseudo-variances are positive (log-concave likelihood)

    C_inv <- inv_lemma_get(Kz, Kxz, S)
    alpha <- inv_lemma_solve(C_inv, z)
    fhat_new <- Kxz %*% backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz) %*% alpha)) + D * alpha

    C_logdet <- C_inv$logdet
    V_logdet <- sum(log(V))
    invC_V_logdet <- V_logdet - C_logdet

    f_invK_f <- t(fhat_new) %*% alpha
    log_lik <- get_loglik(gp$lik, fhat_new, y, ...)
    log_evidence <- log_lik - 0.5 * f_invK_f + 0.5 * invC_V_logdet
    log_evidence <- as.numeric(log_evidence)
    log_post_unnorm <- as.numeric(log_lik - 0.5 * f_invK_f)
    
  } else {

    # non log-concave likelihood; needs special treatment.
    # evaluate the marginal likelihood as p(y) = p(y1,y2) = p(y1)*p(y2 | y1)
    # where y1 are the good data (for which V > 0), and y2 are the bad ones
    good <- which(V > 0)
    bad <- which(V <= 0)
    nbad <- length(bad)

    if (nbad > 500) {
      warning("More than 500 data points with pseudo-variance < 0 (outliers), computation getting slow. Check that the model and hyperparameter values are sensible.")
    }

    args <- c(
      list(
        object = object, gp = gp, Kz = Kz, Kz_chol = Kz_chol, Kxz = Kxz[good, , drop = FALSE],
        D = D[good], y = y[good], fhat_old = NULL, pobs = list(z = z[good], var = V[good])
      ),
      list(...)
    )
    for (argname in required_extra_args(gp$lik)) {
      args[[argname]] <- args[[argname]][good]
    }
    fit1 <- do.call(laplace_iter, args)
    log_post_unnorm1 <- fit1$log_post_unnorm
    log_evidence1 <- fit1$log_evidence

    # compute predictive mean and covariance for f2, given posterior for y1, and
    # then use these as 'prior' mean and covariance to compute marginal likelihood
    # of y2
    K21 <- t(forwardsolve(Kz_chol, t(Kxz[bad, , drop = FALSE]))) %*%
      forwardsolve(Kz_chol, t(Kxz[good, , drop = FALSE]))
    Linv_K2z <- forwardsolve(Kz_chol, t(Kxz[bad, , drop = FALSE]))
    prior_cov <- t(Linv_K2z) %*% Linv_K2z + diag(D[bad], nrow = nbad)
    pred_mean <- as.vector(K21 %*% fit1$alpha)
    pred_cov <- prior_cov - K21 %*% inv_lemma_solve(fit1$C_inv, t(K21))

    # if the following Cholesky computations fail, then the posterior of the whole
    # data is not proper (the covariance is not positive definite)
    Lk <- t(chol(pred_cov))
    A <- diag(nbad) + t(Lk) %*% diag_times_dense(1 / V[bad], Lk)
    A_chol <- tryCatch(
      {
        t(chol(A))
      },
      error = function(err) {
        NULL
      }
    )
    if (is.null(A_chol)) {
      return(list(fmean = fhat_old, log_evidence = -Inf))
    }
    alpha2 <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% (z[bad] / V[bad]))))

    aux <- backsolve(t(Lk), forwardsolve(Lk, pred_mean))
    aux <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% aux)))
    f2hat <- pred_cov %*% (alpha2 + aux)

    # logdet((inv(pred_cov + V)*V) = logdet((inv(inv(V)*pred_cov + I))
    invC2_V2_logdet <- 2 * sum(log(diag(A_chol)))

    df2_invK2_df2 <- sum(forwardsolve(Lk, f2hat - pred_mean)^2)
    args <- c(list(object = gp$lik, f = f2hat, y = y[bad]), list(...))
    for (argname in required_extra_args(gp$lik)) {
      args[[argname]] <- args[[argname]][bad]
    }
    log_lik2 <- do.call(get_loglik, args)
    log_post_unnorm2 <- as.numeric(log_lik2 - 0.5 * df2_invK2_df2)
    log_evidence2 <- log_lik2 - 0.5 * df2_invK2_df2 + 0.5 * invC2_V2_logdet
    log_evidence2 <- as.numeric(log_evidence2)
    
    
    log_post_unnorm <- log_post_unnorm1 + log_post_unnorm2
    log_evidence <- log_evidence1 + log_evidence2

    # finally, fit the model using all the observations
    C_inv <- inv_lemma_get(Kz, Kxz, S, logdet = FALSE)
    alpha <- inv_lemma_solve(C_inv, z)
    fhat_new <- Kxz %*% backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz) %*% alpha)) + D * alpha
  }

  list(
    fmean = fhat_new, z = z, V = V, Kz = Kz, Kxz = Kxz, Kz_chol = Kz_chol, C_inv = C_inv,
    diag = D, alpha = alpha, log_evidence = log_evidence, log_post_unnorm = log_post_unnorm
  )
}


laplace_iter.method_rf <- function(object, gp, Z, y, fhat_old, ...) {

  # calculate first the new estimate for posterior mean for f
  pobs <- get_pseudodata_la(gp$lik, fhat_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  d <- NCOL(Z)
  wprec_prior <- diag(d)
  H <- t(Z) %*% (Z / V)
  L <- t(chol(H + wprec_prior))
  what <- backsolve(t(L), forwardsolve(L, t(Z) %*% (z / V)))
  wcov <- backsolve(t(L), forwardsolve(L, diag(d)))
  wcov_logdet <- -2 * sum(log(diag(L)))
  fhat_new <- Z %*% what

  # compute the log marginal likelihood
  log_prior <- -0.5 * d * log(2 * pi) - 0.5 * sum(what^2)
  log_lik <- get_loglik(gp$lik, fhat_new, y, ...)
  log_evidence <- 0.5 * n * log(2 * pi) + 0.5 * wcov_logdet + log_prior + log_lik
  list(fmean = fhat_new, wmean = what, wcov = wcov, z = z, V = V, log_evidence = log_evidence)
}

laplace_iter.method_rbf <- function(object, gp, Z, y, fhat_old, ...) {
  laplace_iter.method_rf(object, gp, Z, y, fhat_old, ...)
}







laplace.method_full <- function(object, gp, K, y, maxiter = 100, tol = 1e-4, ...) {

  # this is the newton iteration, so iterate the second order approximation
  # to the log likelihood by updating the mean until convergence

  n <- length(y)
  fhat <- rep(0, n)
  if ("lik_gaussian" %in% class(gp$lik)) {
    maxiter <- 1
  }

  for (iter in 1:maxiter) {
    fit <- laplace_iter(object, gp, K, y, fhat, ...)
    diff <- max(abs(fhat - fit$fmean))
    fhat <- fit$fmean
    if (diff < tol) {
      break
    }
  }
  if (maxiter > 1 && iter == maxiter) {
    warning("Maximum number of iterations in Laplace reached, max delta f = ", diff)
  }
  return(fit)
}

laplace.method_fitc <- function(object, gp, Kz, Kz_chol, Kxz, D, y, 
                                maxiter = 100, tol = 1e-4, ...) {
  n <- length(y)
  fhat <- rep(0, n)
  if ("lik_gaussian" %in% class(gp$lik)) {
    maxiter <- 1
  }
  
  log_post_old <- -Inf

  for (iter in 1:maxiter) {
    
    fit <- laplace_iter(object, gp, Kz, Kz_chol, Kxz, D, y, fhat, ...)
    fhat <- fit$fmean
    diff_log_post <- fit$log_post_unnorm - log_post_old
    log_post_old <- fit$log_post_unnorm
    
    if (abs(diff_log_post) < tol) {
      break
    }
  }
  
  if (maxiter > 1 && iter == maxiter) {
    warning("Maximum number of iterations in Laplace reached, delta log post = ", diff_log_post)
  }
  return(fit)
}

laplace.method_rf <- function(object, gp, Z, y, maxiter = 100, tol = 1e-4, ...) {

  # this is the newton iteration, so iterate the second order approximation
  # to the log likelihood by updating the mean until convergence

  n <- length(y)
  fhat <- rep(0, n)
  if ("lik_gaussian" %in% class(gp$lik)) {
    maxiter <- 1
  }

  for (iter in 1:maxiter) {
    fit <- laplace_iter(object, gp, Z, y, fhat, ...)
    fhat_new <- Z %*% fit$wmean
    diff <- max(abs(fhat - fhat_new))
    fhat <- fhat_new
    if (diff < tol) {
      break
    }
  }
  if (maxiter > 1 && iter == maxiter) {
    warning("Maximum number of iterations in Laplace reached, max delta f = ", diff)
  }
  return(fit)
}

laplace.method_rbf <- function(object, gp, Z, y, maxiter = 100, tol = 1e-4, ...) {
  laplace.method_rf(object, gp, Z, y, maxiter = maxiter, tol = tol, ...)
}
