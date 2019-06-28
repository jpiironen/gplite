context('basic')



# create some data
set.seed(1234)
n <- 20
nt <- 30
x <- rnorm(n)
xt <- rnorm(nt)
y <- x^2 + 0.1*rnorm(n)


# initialize model
gp <- gp_init()
gp <- gp_fit(gp,x,y)



#test_that("gp_init: the model has relevant fields", {
#})


test_that("gp_fit: marginal likelihood is correctly calculated", {
  
  set.seed(2309)
  
  param0 <- get_param(gp)
  
  nperturb <- 10
  scale_perturb <- 0.3
  for (j in 1:nperturb) {
    
    # fit the model with the new hyperparameter values
    param <- param0 + rnorm(length(param0))*scale_perturb
    gp <- set_param(gp, param)
    gp <- gp_fit(gp,x,y)
    log_evidence <- gp$log_evidence
    
    # analytic marginal likelihood for Gaussian likelihood (see Rasmussen and Williams, 2006)
    Ky <- gp$K + gp$lik$sigma^2*diag(n)
    L <- t(chol(Ky))
    aux <- solve(L,y)
    log_evidence_analytic <- c( -0.5*t(aux) %*% aux -sum(log(diag(L))) -n/2*log(2*pi) )
    
    expect_equal(log_evidence, log_evidence_analytic, tol=1e-3)
  }
})




test_that("gp_fit: predictions are correctly calculated", {
  
  set.seed(2209)
  
  param0 <- get_param(gp)
  
  nperturb <- 10
  scale_perturb <- 0.3
  for (j in 1:nperturb) {
    
    # fit the model with the new hyperparameter values
    param <- param0 + rnorm(length(param0))*scale_perturb
    gp <- set_param(gp, param)
    gp <- gp_fit(gp,x,y)
    pred <- gp_pred(gp,xt,var=T)
    
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
})




# TODO: test that the mcmc and analytic prediction are the same for the gaussian model




