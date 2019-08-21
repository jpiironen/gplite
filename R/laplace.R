




approx_laplace <- function(gp, K, y, fhat_old, ...) {
  
  # calculate first the new estimate for posterior mean for f
  pobs <- get_pseudodata(gp$lik, fhat_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  K_chol <- t(chol(K)) # TODO: this might be numerically unstable, but perhaps this is not needed, since it cancels out when computing marginal likelihood?
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


approx_laplace_iterated <- function(gp, K, y, maxiter=100, tol=1e-4, ...) {
  
  # this is the newton iteration, so iterate the second order approximation
  # to the log likelihood by updating the mean until convergence
  
  n <- length(y)
  fhat <- rep(0, n)
  if ('lik_gaussian' %in% class(gp$lik))
    maxiter <- 1
  
  for (iter in 1:maxiter) {
    approx <- approx_laplace(gp, K, y, fhat, ...)
    diff <- max(abs(fhat - approx$fmean))
    fhat <- approx$fmean
    if (diff < tol)
      break
  }
  if (maxiter > 1 && iter == maxiter)
    warning('Maximum number of iterations in Laplace reached, results can be unreliable.')
  return(approx)
}







approx_linearized_laplace <- function(gp, z, y, fhat_old, ...) {
  
  # calculate first the new estimate for posterior mean for f
  pobs <- get_pseudodata(gp$lik, fhat_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  K_chol <- t(chol(K))
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


