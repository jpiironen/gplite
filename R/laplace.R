




approx_laplace <- function(gp, y, fhat_old) {
  
  pobs <- get_pseudodata(gp$lik, fhat_old, y)
  z <- pobs$z
  V <- pobs$var
  n <- length(z)
  K <- gp$K
  K_chol <- t(chol(K))
  C_chol <- t(chol(K+diag(V)))
  fhat_new <- K %*% solve(t(C_chol), solve(C_chol, z)) # new estimate for posterior mean for f
  K_logdet <- 2*sum(log(diag(K_chol)))
  C_logdet <- 2*sum(log(diag(C_chol)))
  V_logdet <- sum(log(V))
  aux <- solve(K_chol, fhat_new)
  log_prior <- -0.5*n*log(2*pi) - 0.5*K_logdet - 0.5*sum(aux^2)
  log_lik <- get_loglik(gp$lik, fhat_new, y) 
  log_evidence <- 0.5*n*log(2*pi) + 0.5*(K_logdet - C_logdet + V_logdet) + log_prior + log_lik
  list(fmean=fhat_new, K_chol=K_chol, C_chol=C_chol, log_evidence=log_evidence)
}


approx_laplace_iterated <- function(gp, y, maxiter=100, tol=1e-5, ...) {
  n <- length(y)
  fhat <- rep(0, n)
  if (is.null(gp$K))
    stop('Internal error: covariance K not yet initialized.')
  if ('lik_gaussian' %in% class(gp$lik))
    maxiter <- 1
  
  for (iter in 1:maxiter) {
    approx <- approx_laplace(gp, y, fhat)
    diff <- max(abs(fhat - approx$fmean))
    fhat <- approx$fmean
    if (diff < tol)
      break
  }
  if (maxiter > 1 && iter == maxiter)
    warning('Maximum number of iterations in Laplace reached, results can be unreliable.')
  return(approx)
}



