




ep_iter <- function(object, ...) {
  UseMethod("ep_iter")
}

ep <- function(object, ...) {
  UseMethod("ep", object)
}



ep_iter.approx_full <- function(object, gp, K, y, fmean_old, fvar_old, z_old, P_old, ...) {
  
  # calculate first the new estimate for posterior mean and variance for f
  pobs <- get_pseudodata_ep(gp$lik, fmean_old, 1/fvar_old, z_old, P_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  log_Z <- pobs$log_normalizing_const
  mean_cavity <- pobs$mean_cavity
  var_cavity <- pobs$var_cavity
  n <- length(z)
  B <- diag(n) + diag_times_dense(sqrt(1/V), diag_times_dense(K, sqrt(1/V), diag_right=T))
  C_chol <- diag_times_dense(sqrt(V), t(chol(B)))
  alpha <- backsolve(t(C_chol), forwardsolve(C_chol, z))
  fmean_new <- K %*% alpha
  fcov_new <- diag_times_dense(V, backsolve(t(C_chol), forwardsolve(C_chol, K)))
  
  # compute the log marginal likelihood
  C_logdet <- 2*sum(log(diag(C_chol)))
  z_invC_z <- t(z) %*% alpha
  log_evidence <- -0.5*C_logdet - 0.5*z_invC_z + sum(log_Z) + 
    0.5*sum(log(var_cavity+V)) + 0.5*sum((mean_cavity-z)^2/(var_cavity+V))
  log_evidence <- as.numeric(log_evidence)
  
  list(fmean=as.vector(fmean_new), fvar=diag(fcov_new), z=z, P=1/V, C_chol=C_chol, 
       alpha=alpha, log_evidence=log_evidence)
}

ep_iter.approx_fitc <- function(object, gp, Kz, Kz_chol, Kxz, D, y,
                                fmean_old, fvar_old, z_old, P_old, ...) {
  
  # calculate first the new estimate for posterior mean and variance for f
  pobs <- get_pseudodata_ep(gp$lik, fmean_old, 1/fvar_old, z_old, P_old, y, ...)
  z <- pobs$z
  V <- pobs$var
  log_Z <- pobs$log_normalizing_const
  mean_cavity <- pobs$mean_cavity
  var_cavity <- pobs$var_cavity
  n <- length(z)
  S <- D+V
  
  C_inv <- inv_lemma_get(Kz, Kxz, S)
  alpha <- inv_lemma_solve(C_inv, z)
  inv_Kz_Kzx <- backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz)))
  aux <- inv_lemma_solve(C_inv, V, inv_Kz_Kzx, rhs_diag=T)
  diag1 <- colSums(t(Kxz)*aux)
  diag2 <- inv_lemma_solve(C_inv, V, D, rhs_diag=T, lhs_diag=T, diag_only=T)
  fvar_new <- diag1 + diag2
  fmean_new <- Kxz %*% backsolve(t(Kz_chol), forwardsolve(Kz_chol, t(Kxz) %*% alpha)) + D*alpha

  # compute log marginal likelihood
  C_logdet <- C_inv$logdet
  z_invC_z <- t(z) %*% alpha
  log_evidence <- -0.5*C_logdet - 0.5*z_invC_z + sum(log_Z) + 
    0.5*sum(log(var_cavity+V)) + 0.5*sum((mean_cavity-z)^2/(var_cavity+V))
  log_evidence <- as.numeric(log_evidence)
  
  list(fmean=as.vector(fmean_new), fvar=fvar_new, z=z, P=1/V, pseudovar=V,
       Kz=Kz, Kxz=Kxz, Kz_chol=Kz_chol, C_inv=C_inv, diag=D,
       alpha=alpha, log_evidence=log_evidence)
}





ep.approx_full <- function(object, gp, K, y, maxiter=300, tol=1e-4, 
                           damping=0.9, damping_min=0.1, ...) {
  
  # this is the EP iteration, so iterate the parallel EP until convergence
  
  n <- length(y)
  
  P_old <- rep(0, n)
  z_old <- rep(0, n)
  V_old <- 1/P_old
  fmean_old <- rep(0, n)
  fvar_old <- diag(K)
  diff_old <- Inf
  damp_value <- damping
  
  for (iter in 1:maxiter) {
    fit <- ep_iter(object, gp, K, y, fmean_old, fvar_old, z_old, P_old, 
                   damping=damp_value, ...)
    
    diff_new <- max(abs(fmean_old - fit$fmean))
    if (diff_new > diff_old) {
      damp_value <- damp_value*damping
      next
    }
    damp_value <- max(damp_value, damping_min)
    fmean_old <- fit$fmean
    fvar_old <- fit$fvar
    z_old <- fit$z
    P_old <- fit$P
    diff_old <- diff_new
    if (diff_new < tol)
      break
  }
  
  if (maxiter > 1 && iter == maxiter) {
    fit$log_evidence <- -Inf
    warning('Maximum number of iterations in EP reached, max delta f = ', diff_new) 
  }
  return(fit)
}


ep.approx_fitc <- function(object, gp, Kz, Kz_chol, Kxz, K_diag, D, y, 
                           damping=0.9, damping_min=0.1,
                           maxiter=100, tol=1e-3, ...) {
  
  n <- length(y)
  
  P_old <- rep(0, n)
  z_old <- rep(0, n)
  V_old <- 1/P_old
  fmean_old <- rep(0, n)
  fvar_old <- K_diag
  diff_old <- Inf
  damp_value <- damping
  
  for (iter in 1:maxiter) {
    fit <- ep_iter(object, gp, Kz, Kz_chol, Kxz, D, y, 
                   fmean_old, fvar_old, z_old, P_old, damping=damp_value, ...)
    diff_new <- max(abs(fmean_old - fit$fmean))
    if (diff_new > diff_old) {
      damp_value <- damp_value*damping
      next
    }
    damp_value <- max(damp_value, damping_min)
    fmean_old <- fit$fmean
    fvar_old <- fit$fvar
    z_old <- fit$z
    P_old <- fit$P
    diff_old <- diff_new
    if (diff_new < tol)
      break
  }
  if (maxiter > 1 && iter == maxiter) {
    fit$log_evidence <- -Inf
    warning('Maximum number of iterations in EP reached, max delta f = ', diff_new) 
  }
  return(fit)
}












