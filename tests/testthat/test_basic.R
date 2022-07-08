context('basic')
source(file.path('helpers', 'helpers.R'))


# create some data
set.seed(1234)
n <- 20
nt <- 30
x <- rnorm(n)
xt <- rnorm(nt)
y <- x^2 + 0.2*rnorm(n)


#
cfs <- list(
  cf_const(),
  cf_lin(),
  cf_sexp(),
  cf_matern32(),
  cf_matern52(),
  cf_nn(),
  cf_periodic()
)
lik <- lik_gaussian()


# create some gps with random covariance function combinations
k <- 1
gps <- list()

# all covariance functions alone
for (i in seq_along(cfs)) {
  gp <- gp_init(cfs=cfs[[i]], lik=lik)
  gps[[k]] <- gp_fit(gp,x,y)
  k <- k+1
}

# all pairs of covariance functions
cf_comb <- combn(cfs,2)
for (i in 1:NCOL(cf_comb)) {
  gp <- gp_init(cfs=cf_comb[,i], lik=lik)
  gps[[k]] <- gp_fit(gp,x,y)
  k <- k+1
}

# add products of kernels
cf_comb <- combn(cfs,3)
for (i in 1:NCOL(cf_comb)) {
  SWO(cf <- cf_comb[[1,i]] * cf_comb[[2,i]] * cf_comb[[3,i]])
  gp <- gp_init(cfs=cf, lik=lik)
  gps[[k]] <- gp_fit(gp,x,y)
  k <- k+1
}










test_that("gp_fit: marginal likelihood is correctly calculated", {
  
  set.seed(2309)
  
  for (k in seq_along(gps)) {
    
    gp <- gps[[k]]
    param0 <- get_param(gp)
    
    nperturb <- 10
    scale_perturb <- 0.3
    for (j in 1:nperturb) {
      
      # fit the model with the new hyperparameter values
      param <- param0 + rnorm(length(param0))*scale_perturb
      gp <- set_param(gp, param)
      gp <- gp_fit(gp,x,y)
      log_evidence <- gp$fit$log_evidence
      
      # analytic marginal likelihood for Gaussian likelihood 
      # (see Rasmussen and Williams, 2006)
      K <- eval_cf(gp$cfs, x, x)
      C <- K + gp$lik$sigma^2*diag(n)
      L <- t(chol(C))
      aux <- solve(L,y)
      log_evidence_analytic <- c( -0.5*t(aux) %*% aux -sum(log(diag(L))) -n/2*log(2*pi) )
      
      expect_equal(log_evidence, log_evidence_analytic, tol=1e-3)
    }
  }
})




test_that("gp_fit: predictions are correctly calculated", {
  
  set.seed(2209)
  
  for (k in seq_along(gps)) {
    gp <- gps[[k]]
    param0 <- get_param(gp)
    
    nperturb <- 10
    scale_perturb <- 0.3
    for (j in 1:nperturb) {
      
      # fit the model with the new hyperparameter values
      param <- param0 + rnorm(length(param0))*scale_perturb
      gp <- set_param(gp, param)
      gp <- gp_fit(gp,x,y)
      pred <- gp_pred(gp,xt,var=T, jitter = 1e-6)
      
      # analytic predictive equations mean and variance
      K <- eval_cf(gp$cfs, x, x)
      Kt <- eval_cf(gp$cfs, xt, x)
      Ktt <- eval_cf(gp$cfs, xt, xt)
      C <- K + gp$lik$sigma^2*diag(n)
      L <- t(chol(C))
      alpha <- solve(t(L),solve(L,y))
      aux <- solve(L,t(Kt))
      pred_mean_analytic <- c(Kt %*% alpha)
      pred_std_analytic <- sqrt(diag(Ktt - t(aux) %*% aux))
      
      expect_equal(pred$mean, pred_mean_analytic, tol=1e-3)
      expect_equal(sqrt(pred$var), pred_std_analytic, tol=1e-3)
    }
  }
})




































