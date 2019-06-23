
#rstan_options(auto_write = TRUE)

# TODO: FOR SOME REASON THIS DOES NOT WORK IF RSTAN IS NOT LOADED IN THE SESSION

PATH <- 'stan/'
stanmodels <- list(
  gaussian = rstan::stan_model(paste0(PATH, 'gp_full_gaussian.stan')),
  binomial_logit = rstan::stan_model(paste0(PATH, 'gp_full_binomial_logit.stan')),
  binomial_probit = rstan::stan_model(paste0(PATH, 'gp_full_binomial_probit.stan'))
)
rm(PATH)