context('param')
source(file.path('helpers', 'helpers.R'))


# create some data
set.seed(1234)
n <- 30
nt <- 30
x <- runif(n)*6-3
xt <- runif(nt)*6-3
f <- 2*x - 2 
trials <- sample(10, n, replace = T)



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

liks <- list(
  lik_gaussian(), 
  lik_bernoulli('logit'),
  lik_bernoulli('probit'),
  lik_binomial('logit'), 
  lik_binomial('probit'),
  lik_betabinom('logit'), 
  lik_betabinom('probit'),
  lik_poisson()
)


# create some gps 
k <- 1
gps <- list()
yval <- list()

# loop through the likelihoods
for (j in seq_along(liks)) {
  
  # all covariance functions alone
  for (i in seq_along(cfs)) {
    gps[[k]] <- gp_init(cfs=cfs[[i]], lik=liks[[j]])
    yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
    k <- k+1
  }
  
  # all pairs of covariance functions
  cf_comb <- combn(cfs,2)
  for (i in 1:NCOL(cf_comb)) {
    gps[[k]] <- gp_init(cfs=cf_comb[,i], lik=liks[[j]])
    yval[[k]] <- generate_target(gps[[k]], f, trials=trials)
    k <- k+1
  }
  
  # add products of kernels
  cf_comb <- combn(cfs,3)
  for (i in 1:NCOL(cf_comb)) {
    SWO(cf <- cf_comb[[1,i]] * cf_comb[[2,i]] * cf_comb[[3,i]])
    gps[[k]] <- gp_init(cfs=cf, lik=liks[[j]])
    k <- k+1
  }
}



test_that("set_param: setting the parameters manually has the expected effect", {
            
  for (k in seq_along(gps)) {
    
    gp <- gps[[k]]
    
    # set the parameters to new values and check that the result is correct
    param0 <- get_param(gp)
    param <- runif(length(param0))
    names(param) <- names(param0)
    expect_equal(get_param(set_param(gp, param)), param)
  }
})


