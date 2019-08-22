

data {
  int<lower=0> n;
  int<lower=0> trials[n];
  real<lower=0> phi; // dispersion parameter
  int<lower=0> y[n];
  int<lower=0> link;
}

parameters {
  vector[n] f;
}

transformed parameters {
  vector[n] mu;
  vector[n] alpha;
  vector[n] beta;
  if (link == 0)
    mu = inv_logit(f);
  else if (link == 1)
    mu = Phi_approx(f);
  alpha = mu/phi;
  beta = (1-mu)/phi;
}

model {
  target += beta_binomial_lpmf(y | trials, alpha, beta);
}

