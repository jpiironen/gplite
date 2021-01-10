context('sparse')
source(file.path('helpers', 'helpers.R'))


# create some data
set.seed(1234)
n <- 30
nt <- 5
x <- runif(n)*6-3
xt <- runif(nt)*6-3
f <- 2*x - 4
trials <- sample(10, n, replace = T)



# 
cfs <- list(
  cf_const(magn=0.1), 
  cf_lin(magn=0.1), 
  cf_sexp(magn=0.1),
  cf_nn(magn=0.1),
  cf_periodic()
) 

liks <- list(
  lik_gaussian(),
  lik_bernoulli('logit'),
  lik_binomial('logit'),
  lik_betabinom('logit'),
  lik_poisson()
)

if (exists('GPLITE_TEST_EXTENSIVE') && GPLITE_TEST_EXTENSIVE) {
  extra_liks <- list(
    lik_bernoulli('probit'),
    lik_binomial('probit'),
    lik_betabinom('probit')
  )
  liks <- c(liks, extra_liks)
}


methods <- list(
  method_rf(num_basis = 20),
  method_fitc(num_inducing = 10)
)

# create some gps 
k <- 1
gps <- list()
yval <- list()


# loop through the methods
for (m in seq_along(methods)) {
  
  # loop through the likelihoods
  for (j in seq_along(liks)) {
    
    # all covariance functions alone
    for (i in seq_along(cfs)) {
      gps[[k]] <- gp_init(cfs=cfs[[i]], lik=liks[[j]], method=methods[[m]])
      yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
      k <- k+1
    }
    
    if (exists('GPLITE_TEST_EXTENSIVE') && GPLITE_TEST_EXTENSIVE) {
      
      # additional combinations if extensive tests are desired
      
      # all pairs of covariance functions
      cf_comb <- combn(cfs,2)
      for (i in 1:NCOL(cf_comb)) {
        gps[[k]] <- gp_init(cfs=cf_comb[,i], lik=liks[[j]], method=methods[[m]])
        yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
        k <- k+1
      }
      
      ## add products of kernels
      cf_comb <- combn(cfs,3)
      for (i in 1:NCOL(cf_comb)) {
        cf <- cf_comb[[1,i]] * cf_comb[[2,i]] * cf_comb[[3,i]]
        if ('cf_const' %in% sapply(cf$cfs, class) ||
            'cf_lin' %in% sapply(cf$cfs, class) )
          next
        gps[[k]] <- gp_init(cfs=cf, lik=liks[[j]], method=methods[[m]])
        yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
        k <- k+1
      }
    }
  }
}





test_that("gp_pred: error is raised (only) if model has not been refitted after 
          resetting hyperparameters", {
            
  for (k in seq_along(gps)) {

    #print(paste0('k = ', k))
    gp0 <- gps[[k]]
    gp <- gp_fit(gp0, x, yval[[k]], trials=trials)
    
    if ('method_rf' %in% class(gp0$method)) {
      # prior prediction, should work fine
      expect_silent(gp_draw(gp0, x, draws = 1))
      expect_silent(gp_pred(gp0, x, var=T))
      expect_silent(gp_pred(gp0, x, var=F))
    }
    
    if ('method_fitc' %in% class(gp0$method)) {
      # these should raise error because inducing points not set yet
      expect_error(gp_draw(gp0, x, draws = 1))
      expect_error(gp_pred(gp0, x, var=T))
      # this should work
      expect_silent(gp_pred(gp0, x, var=F))
    }
    
    
    # should work fine
    expect_silent(gp_draw(gp, x, draws = 1))
    expect_silent(gp_pred(gp, x, var=T))
    expect_silent(gp_pred(gp, x, var=F))
    
    # reset one of the hyperparameters
    param <- get_param(gp)
    param[1] <- -2.0
    gp <- set_param(gp, param)
    
    # these should raise an error
    expect_error(gp_draw(gp, x, draws=1))
    expect_error(gp_pred(gp, x, var=T))
    expect_error(gp_pred(gp, x, var=F))
    
    # refit the model
    gp1 <- gp_fit(gp,x,yval[[k]], trials=trials)
    
    # these should again work fine
    expect_silent(gp_draw(gp1, x, draws = 1))
    expect_silent(gp_pred(gp1, x, var=T))
    expect_silent(gp_pred(gp1, x, var=F))
    
    
  }
})




test_that("gp_pred: analytic prediction gives the same result as the sampling 
          based prediction", {
            
  for (k in seq_along(gps)) {
    
    # fit model
    gp <- gp_fit(gps[[k]], x, yval[[k]], trials=trials)
    
    # analytic prediction for f at test points
    pred <- gp_pred(gp,xt, var=T)
    
    # sampling based prediction
    draws <- gp_draw(gp,xt,draws=1e5, transform=F, seed=4321)
    
    expect_equal(rowMeans(draws), pred$mean, tol=1e-2)
    expect_equal(apply(draws, 1, sd),  sqrt(pred$var), tol=1e-2)
  }
})































