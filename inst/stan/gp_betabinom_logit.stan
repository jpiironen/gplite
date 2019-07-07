

data {
  int<lower=0> n;
  int<lower=0> trials[n];
  real<lower=0> phi; // dispersion parameter
  matrix[n,n] L; // cholesky of the covariance matrix (lower triangular) 
  int<lower=0> y[n];
}

parameters {
  vector[n] f_white;
}

transformed parameters {
  vector[n] f;
  vector[n] mu;
  vector[n] alpha;
  vector[n] beta;
  f = L*f_white;
  mu = inv_logit(f);
  alpha = mu/phi;
  beta = (1-mu)/phi;
}

model {
  target += normal_lpdf(f_white | 0, 1);
  target += beta_binomial_lpmf(y | trials, alpha, beta);
}

