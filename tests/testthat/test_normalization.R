context('normalize')
source(file.path('helpers', 'helpers.R'))


# create some data
set.seed(1234)
n <- 20
nt <- 30
x <- 5*(rnorm(n) - 1)
x_norm <- scale(x)
xt <- 5*(rnorm(nt) - 1)
xt_norm <- (xt - attr(x_norm,"scaled:center")) / attr(x_norm,"scaled:scale")
y <- 0.1*x^2 + 0.2*rnorm(n)


#
cfs <- list(
  cf_lin(),
  cf_sexp(),
  cf_matern32(),
  cf_matern52(),
  cf_nn()
)
cfs_norm <- list(
  cf_lin(normalize = T),
  cf_sexp(normalize = T),
  cf_matern32(normalize = T),
  cf_matern52(normalize = T),
  cf_nn(normalize = T)
)
lik <- lik_gaussian()


# create some gps with random covariance function combinations
k <- 1
gps <- list()
gps_norm <- list()

# all covariance functions alone
for (i in seq_along(cfs)) {
  gp <- gp_init(cfs=cfs[[i]], lik=lik)
  gp_norm <- gp_init(cfs=cfs_norm[[i]], lik=lik)
  gps[[k]] <- gp_fit(gp,x_norm,y) 
  gps_norm[[k]] <- gp_fit(gp_norm,x,y)
  k <- k+1
}

# all pairs of covariance functions
cf_comb <- combn(cfs,2)
cf_norm_comb <- combn(cfs_norm,2)
for (i in 1:NCOL(cf_comb)) {
  gp <- gp_init(cfs=cf_comb[,i], lik=lik)
  gp_norm <- gp_init(cfs=cf_norm_comb[,i], lik=lik)
  gps[[k]] <- gp_fit(gp,x_norm,y)
  gps_norm[[k]] <- gp_fit(gp_norm,x,y)
  k <- k+1
}

# add products of kernels
cf_comb <- combn(cfs,3)
cf_norm_comb <- combn(cfs_norm,3)
for (i in 1:NCOL(cf_comb)) {
  cf <- cf_comb[[1,i]] * cf_comb[[2,i]] * cf_comb[[3,i]]
  cf_norm <- cf_norm_comb[[1,i]] * cf_norm_comb[[2,i]] * cf_norm_comb[[3,i]]
  gp <- gp_init(cfs=cf, lik=lik)
  gp_norm <- gp_init(cfs=cf_norm, lik=lik)
  gps[[k]] <- gp_fit(gp,x_norm,y)
  gps_norm[[k]] <- gp_fit(gp_norm,x,y)
  k <- k+1
}










test_that("gp_fit: marginal likelihood is correctly calculated", {
  
  set.seed(2309)
  
  for (k in seq_along(gps)) {
    
    gp <- gps[[k]]
    gp_norm <- gps_norm[[k]]
    param0 <- get_param(gp)
    
    nperturb <- 10
    scale_perturb <- 0.3
    for (j in 1:nperturb) {
      
      # fit the models with the new hyperparameter values
      param <- param0 + rnorm(length(param0))*scale_perturb
      gp <- set_param(gp, param)
      gp_norm <- set_param(gp_norm, param)
      gp <- gp_fit(gp,x_norm,y)
      gp_norm <- gp_fit(gp_norm,x,y)
      
      log_evidence1 <- gp$fit$log_evidence
      log_evidence2 <- gp_norm$fit$log_evidence
      
      expect_equal(log_evidence1, log_evidence2, tol=1e-3)
    }
  }
})




test_that("gp_fit: predictions are correctly calculated", {
  
  set.seed(2209)
  
  for (k in seq_along(gps)) {
    gp <- gps[[k]]
    gp_norm <- gps_norm[[k]]
    param0 <- get_param(gp)
    
    nperturb <- 10
    scale_perturb <- 0.3
    for (j in 1:nperturb) {
      
      # fit the models with the new hyperparameter values
      param <- param0 + rnorm(length(param0))*scale_perturb
      gp <- set_param(gp, param)
      gp_norm <- set_param(gp_norm, param)
      gp <- gp_fit(gp,x_norm,y)
      gp_norm <- gp_fit(gp_norm,x,y)
      
      pred1 <- gp_pred(gp, xt_norm, var=T, jitter = 1e-6)
      pred2 <- gp_pred(gp_norm, xt, var=T, jitter = 1e-6)
      
      
      expect_equal(pred1$mean, pred2$mean, tol=1e-5)
      expect_equal(pred1$var, pred2$var, tol=1e-5)
    }
  }
})




































