context('rf')
source(file.path('helpers', 'helpers.R'))


# create some data
set.seed(1234)
n <- 30
nt <- 30
x <- runif(n)*6-3
xt <- runif(nt)*6-3
f <- x^2 - 2 
trials <- sample(10, n, replace = T)



# 
cfs <- list(cf_const(), 
            cf_lin(), 
            cf_sexp(),
            cf_nn())

liks <- list(lik_gaussian(), 
             lik_binomial('logit'), lik_binomial('probit'),
             lik_betabinom('logit'), lik_betabinom('probit'))

method <- 'rf'

# create some gps 
k <- 1
gps <- list()
yval <- list()

# loop through the likelihoods
for (j in seq_along(liks)) {
  
  # all covariance functions alone
  for (i in seq_along(cfs)) {
    gps[[k]] <- gp_init(cfs=cfs[[i]], lik=liks[[j]], method=method)
    yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
    k <- k+1
  }
  
  # all pairs of covariance functions
  cf_comb <- combn(cfs,2)
  for (i in 1:NCOL(cf_comb)) {
    gps[[k]] <- gp_init(cfs=cf_comb[,i], lik=liks[[j]], method=method)
    yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
    k <- k+1
  }
}



test_that("gp_pred: error is raised (only) if model has not been refitted after 
          resetting hyperparameters", {
            
  for (k in seq_along(gps)) {
    gp0 <- gps[[k]]
    gp <- gp_fit(gps[[k]], x, yval[[k]], trials=trials)
    
    # prior predicion, should work fine
    expect_silent(gp_draw(gp0, x, draws = 1))
    expect_silent(gp_pred(gp0, x, var=T))
    expect_silent(gp_pred(gp0, x, var=F))
    
    # should work fine
    expect_silent(gp_draw(gp, x, draws = 1))
    expect_silent(gp_pred(gp, x, var=T))
    expect_silent(gp_pred(gp, x, var=F))
    
    # reset one of the hyperparameters
    param <- get_param(gp)
    param[1] <- 0.0
    gp <- set_param(gp, param)
    
    # these should raise an error
    expect_error(gp_draw(gp, x, draws=1))
    expect_error(gp_pred(gp, x, var=T))
    expect_error(gp_pred(gp, x, var=F))
    
    # refit the model
    gp1 <- gp_fit(gp,x,yval[[k]], trials=trials)
    SWO(gp2 <- gp_mcmc(gp,x,yval[[k]], trials=trials, iter=400, chains=1))
    
    # these should again work fine
    expect_silent(gp_draw(gp1, x, draws = 1))
    expect_silent(gp_pred(gp1, x, var=T))
    expect_silent(gp_pred(gp1, x, var=F))
    expect_silent(gp_draw(gp2, x, draws = 1))
  }
})




test_that("gp_pred: analytic prediction gives the same result as the sampling 
          based prediction", {
            
  for (k in seq_along(gps)) {
    
    # fit model
    gp <- gp_fit(gps[[k]], x, yval[[k]], trials=trials)
    
    # analytic prediction for f at test points
    pred <- gp_pred(gp,xt, var=T, transform=F)
    
    # sampling based prediction
    draws <- gp_draw(gp,xt,draws=1e5, transform=F)
    
    expect_equal(rowMeans(draws), pred$mean, tol=1e-2)
    expect_equal(apply(draws, 1, sd),  sqrt(pred$var), tol=1e-2)
  }
})































