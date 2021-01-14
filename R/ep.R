




ep_iter <- function(object, ...) {
  UseMethod("ep_iter")
}

ep <- function(object, ...) {
  UseMethod("ep", object)
}



ep_iter.method_full <- function(object, gp, K, y, fmean_old, fvar_old, z_old, P_old,
                                pobs = NULL, ...) {

  # get the pseudo data first
  if (is.null(pobs)) {
    pobs <- get_pseudodata_ep(gp$lik, fmean_old, 1 / fvar_old, z_old, P_old, y, ...)
  }
  z <- pobs$z
  V <- pobs$var
  log_Z <- pobs$log_normalizing_const
  mean_cavity <- pobs$mean_cavity
  var_cavity <- pobs$var_cavity
  n <- length(z)

  if (any(is.na(V) | is.na(z))) {
    return(list(fmean_new = fmean_old, log_evidence = -Inf))
  }

  if (all(V > 0)) {

    # normal case; all pseudo-variances are positive (log-concave likelihood)

    B <- diag(n) + diag_times_dense(sqrt(1 / V), diag_times_dense(K, sqrt(1 / V), diag_right = TRUE))
    C_chol <- diag_times_dense(sqrt(V), t(chol(B)))
    alpha <- backsolve(t(C_chol), forwardsolve(C_chol, z))
    fmean_new <- as.vector(K %*% alpha)
    fcov_new <- diag_times_dense(V, backsolve(t(C_chol), forwardsolve(C_chol, K)))

    # compute the log marginal likelihood
    C_logdet <- 2 * sum(log(diag(C_chol)))
    z_invC_z <- t(z) %*% alpha
    log_evidence <- -0.5 * C_logdet - 0.5 * z_invC_z + sum(log_Z) +
      0.5 * sum(log(var_cavity + V)) + 0.5 * sum((mean_cavity - z)^2 / (var_cavity + V))
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
        object = object, gp = gp, K = K[good, good, drop = FALSE], y = y[good],
        fmean_old = NULL, fvar_old = NULL, z_old = NULL, P_old = NULL,
        pobs = list(
          z = z[good], var = V[good], log_normalizing_const = log_Z[good],
          mean_cavity = mean_cavity[good], var_cavity = var_cavity[good]
        )
      ),
      list(...)
    )
    for (argname in required_extra_args(gp$lik)) {
      args[[argname]] <- args[[argname]][good]
    }
    fit1 <- do.call(ep_iter, args)
    log_evidence1 <- fit1$log_evidence

    # compute predictive mean and covariance for f2, given posterior for y1, and
    # then use these as 'prior' mean and covariance to computer marginal likelihood
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
      return(list(fmean = Inf, log_evidence = -Inf))
    }
    alpha2 <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% (z[bad] / V[bad]))))

    aux <- backsolve(t(Lk), forwardsolve(Lk, pred_mean))
    aux <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% aux)))
    f2hat <- pred_cov %*% (alpha2 + aux)

    # logdet((inv(pred_cov + V)*V) = logdet((inv(inv(V)*pred_cov + I))
    invC2_V2_logdet <- 2 * sum(log(diag(A_chol)))

    df2_invC2_df2 <- (pred_mean - z[bad]) %*% backsolve(t(Lk), backsolve(
      t(A_chol),
      forwardsolve(A_chol, t(Lk) %*% ((pred_mean - z[bad]) / V[bad]))
    ))

    log_evidence2 <- -0.5 * df2_invC2_df2 + 0.5 * invC2_V2_logdet + sum(log_Z[bad]) +
      0.5 * sum(log((var_cavity[bad] + V[bad]) / V[bad])) +
      0.5 * sum((mean_cavity[bad] - z[bad])^2 / (var_cavity[bad] + V[bad]))
    log_evidence2 <- as.numeric(log_evidence2)

    log_evidence <- log_evidence1 + log_evidence2

    # finally, fit the model using all the observations.
    # here we use LU-decomposition instead of Cholesky for matrix K+V,
    # as it might not be positive definite
    C_lu <- Matrix::expand(Matrix::lu(K + diag(V)))
    alpha <- solve(C_lu$U, solve(C_lu$L, solve(C_lu$P, z)))
    fmean_new <- as.vector(K %*% alpha)
    fcov_new <- diag_times_dense(V, solve(C_lu$U, solve(C_lu$L, solve(C_lu$P, K))))
    return(list(
      fmean = fmean_new, fvar = diag(fcov_new),
      cavity_mean = mean_cavity, cavity_var = var_cavity,
      z = z, P = 1 / V, C_lu = C_lu, quad_ok = pobs$quad_ok,
      alpha = alpha, log_evidence = log_evidence
    ))
  }

  list(
    fmean = fmean_new, fvar = diag(fcov_new),
    cavity_mean = mean_cavity, cavity_var = var_cavity,
    z = z, P = 1 / V, C_chol = C_chol,
    alpha = alpha, log_evidence = log_evidence, quad_ok = pobs$quad_ok
  )
}

ep_iter.method_fitc <- function(object, gp, Kz, Kz_chol, Kxz, D, y,
                                fmean_old, fvar_old, z_old, P_old, pobs = NULL, ...) {

  # get the pseudo data first
  if (is.null(pobs)) {
    pobs <- get_pseudodata_ep(gp$lik, fmean_old, 1 / fvar_old, z_old, P_old, y, ...)
  }
  z <- pobs$z
  V <- pobs$var
  log_Z <- pobs$log_normalizing_const
  mean_cavity <- pobs$mean_cavity
  var_cavity <- pobs$var_cavity
  n <- length(z)
  S <- D + V

  if (any(is.na(V) | is.na(z))) {
    return(list(fmean_new = fmean_old, log_evidence = -Inf))
  }

  if (all(V > 0)) {

    # normal case; all pseudo-variances are positive (log-concave likelihood)

    C_inv <- inv_lemma_get(Kz, Kxz, S)
    alpha <- inv_lemma_solve(C_inv, z)
    inv_Kz_Kzx <- backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz)))
    aux <- inv_lemma_solve(C_inv, V, inv_Kz_Kzx, rhs_diag = TRUE)
    diag1 <- colSums(t(Kxz) * aux)
    diag2 <- inv_lemma_solve(C_inv, V, D, rhs_diag = TRUE, lhs_diag = TRUE, diag_only = TRUE)
    fvar_new <- diag1 + diag2
    fmean_new <- Kxz %*% backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz) %*% alpha)) + D * alpha

    # compute log marginal likelihood
    C_logdet <- C_inv$logdet
    z_invC_z <- t(z) %*% alpha
    log_evidence <- -0.5 * C_logdet - 0.5 * z_invC_z + sum(log_Z) +
      0.5 * sum(log(var_cavity + V)) + 0.5 * sum((mean_cavity - z)^2 / (var_cavity + V))
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
        object = object, gp = gp,
        Kz = Kz, Kz_chol = Kz_chol, Kxz = Kxz[good, , drop = FALSE], D[good], y[good],
        fmean_old = NULL, fvar_old = NULL, z_old = NULL, P_old = NULL,
        pobs = list(
          z = z[good], var = V[good], log_normalizing_const = log_Z[good],
          mean_cavity = mean_cavity[good], var_cavity = var_cavity[good]
        )
      ),
      list(...)
    )
    for (argname in required_extra_args(gp$lik)) {
      args[[argname]] <- args[[argname]][good]
    }
    fit1 <- do.call(ep_iter, args)
    log_evidence1 <- fit1$log_evidence

    # compute predictive mean and covariance for f2, given posterior for y1, and
    # then use these as 'prior' mean and covariance to computer marginal likelihood
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
      return(list(fmean = Inf, log_evidence = -Inf))
    }
    alpha2 <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% (z[bad] / V[bad]))))

    aux <- backsolve(t(Lk), forwardsolve(Lk, pred_mean))
    aux <- backsolve(t(Lk), backsolve(t(A_chol), forwardsolve(A_chol, t(Lk) %*% aux)))
    f2hat <- pred_cov %*% (alpha2 + aux)

    # logdet((inv(pred_cov + V)*V) = logdet((inv(inv(V)*pred_cov + I))
    invC2_V2_logdet <- 2 * sum(log(diag(A_chol)))

    df2_invC2_df2 <- (pred_mean - z[bad]) %*% backsolve(t(Lk), backsolve(
      t(A_chol),
      forwardsolve(A_chol, t(Lk) %*% ((pred_mean - z[bad]) / V[bad]))
    ))

    log_evidence2 <- -0.5 * df2_invC2_df2 + 0.5 * invC2_V2_logdet + sum(log_Z[bad]) +
      0.5 * sum(log((var_cavity[bad] + V[bad]) / V[bad])) +
      0.5 * sum((mean_cavity[bad] - z[bad])^2 / (var_cavity[bad] + V[bad]))
    log_evidence2 <- as.numeric(log_evidence2)

    log_evidence <- log_evidence1 + log_evidence2

    # finally, fit the model using all the observations
    C_inv <- inv_lemma_get(Kz, Kxz, S, logdet = FALSE)
    alpha <- inv_lemma_solve(C_inv, z)
    inv_Kz_Kzx <- backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz)))
    aux <- inv_lemma_solve(C_inv, V, inv_Kz_Kzx, rhs_diag = TRUE)
    diag1 <- colSums(t(Kxz) * aux)
    diag2 <- inv_lemma_solve(C_inv, V, D, rhs_diag = TRUE, lhs_diag = TRUE, diag_only = TRUE)
    fvar_new <- diag1 + diag2
    fmean_new <- Kxz %*% backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz) %*% alpha)) + D * alpha
  }

  list(
    fmean = as.vector(fmean_new), fvar = fvar_new,
    cavity_mean = mean_cavity, cavity_var = var_cavity, z = z, P = 1 / V,
    Kz = Kz, Kxz = Kxz, Kz_chol = Kz_chol, C_inv = C_inv, diag = D,
    alpha = alpha, log_evidence = log_evidence, quad_ok = pobs$quad_ok
  )
}





ep.method_full <- function(object, gp, K, y, maxiter = 300, tol = 1e-4,
                           damping = 0.9, damping_min = 0.1, ...) {

  # this is the EP iteration, so iterate the parallel EP until convergence

  n <- length(y)

  P_old <- rep(0, n)
  z_old <- rep(0, n)
  V_old <- 1 / P_old
  fmean_old <- rep(0, n)
  fvar_old <- diag(K)
  diff_old <- Inf
  damp_value <- damping

  for (iter in 1:maxiter) {
    fit <- ep_iter(object, gp, K, y, fmean_old, fvar_old, z_old, P_old,
      damping = damp_value, ...
    )

    diff_new <- max(abs(fmean_old - fit$fmean))
    if (diff_new > diff_old) {
      damp_value <- damp_value * damping
      next
    }
    damp_value <- max(damp_value, damping_min)
    fmean_old <- fit$fmean
    fvar_old <- fit$fvar
    z_old <- fit$z
    P_old <- fit$P
    diff_old <- diff_new
    if (diff_new < tol) {
      break
    }
  }

  if (maxiter > 1 && iter == maxiter) {
    fit$log_evidence <- -Inf
    warning("Maximum number of iterations in EP reached, max delta f = ", diff_new)
  }
  if (!is.null(fit$quad_ok) && !fit$quad_ok) {
    fit$log_evidence <- -Inf
  }
  return(fit)
}


ep.method_fitc <- function(object, gp, Kz, Kz_chol, Kxz, K_diag, D, y,
                           damping = 0.9, damping_min = 0.1,
                           maxiter = 100, tol = 1e-3, ...) {
  n <- length(y)

  P_old <- rep(0, n)
  z_old <- rep(0, n)
  V_old <- 1 / P_old
  fmean_old <- rep(0, n)
  fvar_old <- K_diag
  diff_old <- Inf
  damp_value <- damping

  for (iter in 1:maxiter) {
    fit <- ep_iter(object, gp, Kz, Kz_chol, Kxz, D, y,
      fmean_old, fvar_old, z_old, P_old,
      damping = damp_value, ...
    )
    diff_new <- max(abs(fmean_old - fit$fmean))
    if (diff_new > diff_old) {
      damp_value <- damp_value * damping
      next
    }
    damp_value <- max(damp_value, damping_min)
    fmean_old <- fit$fmean
    fvar_old <- fit$fvar
    z_old <- fit$z
    P_old <- fit$P
    diff_old <- diff_new
    if (diff_new < tol) {
      break
    }
  }
  if (maxiter > 1 && iter == maxiter) {
    fit$log_evidence <- -Inf
    warning("Maximum number of iterations in EP reached, max delta f = ", diff_new)
  }
  if (!is.null(fit$quad_ok) && !fit$quad_ok) {
    fit$log_evidence <- -Inf
  }
  return(fit)
}
