

laplace_iter <- function(object, ...) {
  UseMethod("laplace_iter")
}

laplace <- function(object, ...) {
  UseMethod("laplace", object)
}



laplace_iter.approx_full <- function(object, gp, K, y, fhat_old, ...) {
  
  # calculate first the new estimate for posterior mean for f
  pobs <- get_pseudodata(gp$lik, fhat_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  K_chol <- t(chol(K)) # TODO: this could be computed only once, not at every iteration
  C_chol <- t(chol(K+diag(V)))
  fhat_new <- K %*% backsolve(t(C_chol), forwardsolve(C_chol, z)) 
  
  # compute the log marginal likelihood
  K_logdet <- 2*sum(log(diag(K_chol)))
  C_logdet <- 2*sum(log(diag(C_chol)))
  V_logdet <- sum(log(V))
  aux <- forwardsolve(K_chol, fhat_new)
  log_prior <- -0.5*n*log(2*pi) - 0.5*K_logdet - 0.5*sum(aux^2)
  log_lik <- get_loglik(gp$lik, fhat_new, y, ...) 
  log_evidence <- 0.5*n*log(2*pi) + 0.5*(K_logdet - C_logdet + V_logdet) + log_prior + log_lik
  list(fmean=fhat_new, K_chol=K_chol, C_chol=C_chol, log_evidence=log_evidence)
}

laplace_iter.approx_fitc <- function(object, gp, Kz, Kz_chol, Kxz, D, y, fhat_old, ...) {
  
  # calculate first the new estimate for posterior mean for f
  pobs <- get_pseudodata(gp$lik, fhat_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  S <- D+V
  out <- solve_inv_lemma(Kz, Kxz, S, z, log_det=T)
  v <- out[[1]]
  C_logdet <- out[[2]]
  fhat_new <- Kxz %*% backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz) %*% v)) 
  # alternative way:
  aux <- (1/sqrt(S))*Kxz
  Sigma_chol <- t(chol(Kz + t(aux) %*% aux))
  alpha <- backsolve(t(Sigma_chol), forwardsolve(Sigma_chol, t(Kxz) %*% (z/S)))
  #fhat_new <- Kxz %*% alpha
  
  # compute the log marginal likelihood
  K_logdet <- 0 # TODO: this is wrong, but it cancels out in the computation
  V_logdet <- sum(log(V))
  f_invK_f <- t(fhat_new) %*% solve_inv_lemma(Kz, Kxz, S, z)
  log_prior <- -0.5*n*log(2*pi) - 0.5*K_logdet - 0.5*f_invK_f 
  log_lik <- get_loglik(gp$lik, fhat_new, y, ...) 
  log_evidence <- 0.5*n*log(2*pi) + 0.5*(K_logdet - C_logdet + V_logdet) + log_prior + log_lik
  list(fmean=fhat_new, pseudovar=V, Kz=Kz, Kxz=Kxz, Kz_chol=Kz_chol, diag=D, 
       alpha=alpha, log_evidence=log_evidence)
}


solve_inv_lemma <- function(K, U, D, b, log_det=F) {
  # solves system (D + U * K^-1 * U^t)^-1 b using inversion lemma
  # (D should be a vector giving only the diagonal)
  aux <- (1/D)*U 
  aux2 <- (1/sqrt(D))*U
  L <- t(chol(K + t(aux2) %*% aux2))
  v <- b/D - aux %*% backsolve(t(L), forwardsolve(L, t(aux) %*% b))
  if (log_det) {
    # compute also log determinant of the matrix that is to be inverted
    logdet <- sum(log(D)) - 2*sum(log(diag(chol(K)))) + 2*sum(log(diag(L)))
    return(list(v,logdet))
  }
  return(v)
}


laplace_iter.approx_rf <- function(object, gp, Z, y, fhat_old, ...) {
  
  # calculate first the new estimate for posterior mean for f
  pobs <- get_pseudodata(gp$lik, fhat_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  d <- NCOL(Z)
  wprec_prior <- diag(d)
  H <- t(Z) %*% (Z/V)
  L <- t(chol(H + wprec_prior))
  what <- backsolve(t(L), forwardsolve(L, t(Z) %*% (z/V)))
  wcov <- backsolve(t(L), forwardsolve(L, diag(d)))
  wcov_logdet <- -2*sum(log(diag(L)))
  fhat_new <- Z %*% what
  
  # compute the log marginal likelihood
  log_prior <- -0.5*d*log(2*pi) - 0.5*sum(what^2)
  log_lik <- get_loglik(gp$lik, fhat_new, y, ...) 
  log_evidence <- 0.5*n*log(2*pi) + 0.5*wcov_logdet + log_prior + log_lik
  list(wmean=what, wcov=wcov, log_evidence=log_evidence)
}

laplace_iter.approx_rbf <- function(object, gp, Z, y, fhat_old, ...) {
  laplace_iter.approx_rf(object, gp, Z, y, fhat_old, ...)
}







laplace.approx_full <- function(object, gp, K, y, maxiter=100, tol=1e-4, ...) {
  
  # this is the newton iteration, so iterate the second order approximation
  # to the log likelihood by updating the mean until convergence
  
  n <- length(y)
  fhat <- rep(0, n)
  if ('lik_gaussian' %in% class(gp$lik))
    maxiter <- 1
  
  for (iter in 1:maxiter) {
    fit <- laplace_iter(object, gp, K, y, fhat, ...)
    diff <- max(abs(fhat - fit$fmean))
    fhat <- fit$fmean
    if (diff < tol)
      break
  }
  if (maxiter > 1 && iter == maxiter)
    warning('Maximum number of iterations in Laplace reached, results can be unreliable.')
  return(fit)
}

laplace.approx_fitc <- function(object, gp, Kz, Kz_chol, Kxz, D, y, maxiter=100, tol=1e-4, ...) {
  n <- length(y)
  fhat <- rep(0, n)
  if ('lik_gaussian' %in% class(gp$lik))
    maxiter <- 1
  
  for (iter in 1:maxiter) {
    #fit <- laplace_iter(object, gp, K, y, fhat, ...)
    fit <- laplace_iter(object, gp, Kz, Kz_chol, Kxz, D, y, fhat, ...)
    diff <- max(abs(fhat - fit$fmean))
    fhat <- fit$fmean
    if (diff < tol)
      break
  }
  if (maxiter > 1 && iter == maxiter)
    warning('Maximum number of iterations in Laplace reached, results can be unreliable.')
  return(fit)
}

laplace.approx_rf <- function(object, gp, Z, y, maxiter=100, tol=1e-4, ...) {
  
  # this is the newton iteration, so iterate the second order approximation
  # to the log likelihood by updating the mean until convergence
  
  n <- length(y)
  fhat <- rep(0, n)
  if ('lik_gaussian' %in% class(gp$lik))
    maxiter <- 1
  
  for (iter in 1:maxiter) {
    fit <- laplace_iter(object, gp, Z, y, fhat, ...)
    fhat_new <- Z %*% fit$wmean
    diff <- max(abs(fhat - fhat_new))
    fhat <- fhat_new
    if (diff < tol)
      break
  }
  if (maxiter > 1 && iter == maxiter)
    warning('Maximum number of iterations in Laplace reached, results can be unreliable.')
  return(fit)
}

laplace.approx_rbf <- function(object, gp, Z, y, maxiter=100, tol=1e-4, ...) {
  laplace.approx_rf(object, gp, Z, y, maxiter=maxiter, tol=tol, ...)
}












