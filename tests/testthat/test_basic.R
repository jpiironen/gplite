context('basic')
source(file.path('helpers', 'helpers.R'))


# create some data
set.seed(1234)
n <- 20
nt <- 30
x <- rnorm(n)
xt <- rnorm(nt)
y <- x^2 + 0.1*rnorm(n)


# 
cfs <- list(cf_const(), 
            cf_lin(), 
            cf_sexp(),
            cf_matern32(), 
            cf_matern52(), 
            cf_nn())


# create some gps with random covariance function combinations
k <- 1
gps <- list()
for (ncf in 1:length(cfs)) {
  for (nrep in 1:3) {
    gp <- gp_init(cfs=sample(cfs, ncf), lik=lik_gaussian())
    gps[[k]] <- gp_fit(gp,x,y)
    k <- k+1
  }
}

# these are gps which are only initialized, not fitted
k <- 1
gps_notfitted <- list()
for (ncf in 1:length(cfs)) {
  for (nrep in 1:3) {
    gps_notfitted[[k]] <- gp_init(cfs=sample(cfs, ncf), lik=lik_gaussian())
    k <- k+1
  }
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
      log_evidence <- gp$log_evidence
      
      # analytic marginal likelihood for Gaussian likelihood 
      # (see Rasmussen and Williams, 2006)
      Ky <- gp$K + gp$lik$sigma^2*diag(n)
      L <- t(chol(Ky))
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
      Ktt <- eval_cf(gp$cfs, xt, xt)
      Kt <- eval_cf(gp$cfs, xt, x)
      Ky <- gp$K + gp$lik$sigma^2*diag(n)
      L <- t(chol(Ky))
      alpha <- solve(t(L),solve(L,y))
      aux <- solve(L,t(Kt))
      pred_mean_analytic <- c(Kt %*% alpha)
      pred_std_analytic <- sqrt(diag(Ktt - t(aux) %*% aux))
      
      expect_equal(pred$mean, pred_mean_analytic, tol=1e-3)
      expect_equal(sqrt(pred$var), pred_std_analytic, tol=1e-3)
    }
  }
})


# TODO: test that the mcmc and analytic prediction are the same for the gaussian model


test_that("gp_pred: error is raised (only) if model has not been refitted after 
          resetting hyperparameters", {
  
  for (k in seq_along(gps)) {
    gp0 <- gps_notfitted[[k]]
    gp <- gps[[k]]
    
    # prior predicion, should work fine
    expect_silent(gp_pred(gp0, x, draws = 1))
    expect_silent(gp_pred(gp0, x, var=T))
    expect_silent(gp_pred(gp0, x, var=F))
    
    # should work fine
    expect_silent(gp_pred(gp, x, draws = 1))
    expect_silent(gp_pred(gp, x, var=T))
    expect_silent(gp_pred(gp, x, var=F))
    
    # reset one of the hyperparameters
    param <- get_param(gp)
    param[1] <- 0.0
    gp <- set_param(gp, param)
    
    # these should raise an error
    expect_error(gp_pred(gp, x, draws=1))
    expect_error(gp_pred(gp, x, var=T))
    expect_error(gp_pred(gp, x, var=F))
    
    # refit the model
    gp1 <- gp_fit(gp,x,y)
    SWO(gp2 <- gp_sample(gp,x,y, iter=400, chains=1))
    
    # these should again work fine
    expect_silent(gp_pred(gp1, x, draws = 1))
    expect_silent(gp_pred(gp1, x, var=T))
    expect_silent(gp_pred(gp1, x, var=F))
    expect_silent(gp_pred(gp2, x, draws = 1))
    expect_silent(gp_pred(gp2, x, var=T))
    expect_silent(gp_pred(gp2, x, var=F))
  }
})




























