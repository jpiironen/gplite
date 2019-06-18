
#rstan_options(auto_write = TRUE)

STANFILE_PATH <- '../stan/'
stanmodels <- list(
  gaussian = stan_model(paste0(STANFILE_PATH, 'gp_full_gaussian.stan')),
  binomial_logit = stan_model(paste0(STANFILE_PATH, 'gp_full_binomial_logit.stan')),
  binomial_probit = stan_model(paste0(STANFILE_PATH, 'gp_full_binomial_probit.stan'))
)
rm(STANFILE_PATH)