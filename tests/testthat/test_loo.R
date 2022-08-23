context('loo')
source(file.path('helpers', 'helpers.R'))


# create some data
set.seed(1234)
n <- 30
x <- runif(n)*6-3
f <- 2*x - 4
trials <- sample(10, n, replace = T)



# 
cfs <- list(
  cf_const(), 
  cf_lin(), 
  cf_sexp(),
  #cf_matern32(), 
  #cf_matern52(), 
  cf_nn(),
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
  method_full(),
  method_rf(num_basis = 20),
  method_fitc(num_inducing = 10)
)

approx <- list(
  approx_laplace(),
  approx_ep(damping = 0.8)
)

# create some gps 
k <- 1
gps <- list()
yval <- list()

# loop through the methods
for (m in seq_along(methods)) {
  
  # loop through approxmations
  for (a in seq_along(approx)) {
    
    if ((class(approx[[a]])[1] == 'approx_ep') & (class(methods[[m]])[1] == 'method_rf'))
      # no ep-implementation for rf yet
      next
    
    # loop through the likelihoods
    for (j in seq_along(liks)) {
      
      if ("approx_ep" %in% class(approx[[a]]) & "lik_gaussian" %in% class(liks[[j]]))
        # ep not allowed for gaussian likelihood
        next
      
      # all covariance functions alone
      for (i in seq_along(cfs)) {
        gps[[k]] <- gp_init(cfs=cfs[[i]], lik=liks[[j]], method=methods[[m]], approx=approx[[a]])
        yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
        k <- k+1
      }
      
      if (exists('GPLITE_TEST_EXTENSIVE') && GPLITE_TEST_EXTENSIVE) {
        
        # additional combinations if extensive tests are desired
        
        # all pairs of covariance functions
        cf_comb <- combn(cfs,2)
        for (i in 1:NCOL(cf_comb)) {
          gps[[k]] <- gp_init(cfs=cf_comb[,i], lik=liks[[j]], method=methods[[m]], approx=approx[[a]])
          yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
          k <- k+1
        }
        
        ## add products of kernels
        cf_comb <- combn(cfs,3)
        for (i in 1:NCOL(cf_comb)) {
          SWO(cf <- cf_comb[[1,i]] * cf_comb[[2,i]] * cf_comb[[3,i]])
          if ('cf_const' %in% sapply(cf$cfs, class) ||
              'cf_lin' %in% sapply(cf$cfs, class) )
            next
          gps[[k]] <- gp_init(cfs=cf, lik=liks[[j]], method=methods[[m]], approx=approx[[a]])
          yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
          k <- k+1
        }
        
      }
    }
  }
}




test_that("gp_loo: error is raised if model is not fitted yet, and fitted model does not give erros", {
            
  for (k in seq_along(gps)) {
    gp0 <- gps[[k]]
    SWO(gp <- gp_fit(gps[[k]], x, yval[[k]], trials=trials))
    
    # not fitted, should give error
    expect_error(gp_loo(gp0, x, trials=trials))
    
    # fitted, should work fine
    expect_silent(loo <- gp_loo(gp, x, yval[[k]], trials=trials, quad_order = 7))
    expect_false(is.na(loo$loo))
    expect_false(any(is.na(loo$loos)))
    
    # reset one of the hyperparameters
    param <- get_param(gp)
    param[1] <- 0.0
    gp <- set_param(gp, param)
    
    # this should raise an error
    expect_error(gp_loo(gp, x, yval[[k]], trials=trials))
    
    # refit the model
    SWO(gp1 <- gp_fit(gp, x, yval[[k]], trials=trials))
    
    # these should again work fine
    expect_silent(loo <- gp_loo(gp1, x, yval[[k]], trials=trials, quad_order = 7))
    expect_false(is.na(loo$loo))
    expect_false(any(is.na(loo$loos)))
  }
})





























