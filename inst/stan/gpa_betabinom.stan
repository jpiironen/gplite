

data {
  int<lower=0> n;
  int<lower=0> m;
  int<lower=0> trials[n];
  real<lower=0> phi; // dispersion parameter
  matrix[n,m] Z; 
  int<lower=0> y[n];
  int<lower=0> link;
}

parameters {
  vector[m] w;
}

transformed parameters {
  vector[n] mu;
  vector[n] alpha;
  vector[n] beta;
  if (link == 0)
    mu = inv_logit(Z*w);
  else if (link == 1)
    mu = Phi_approx(Z*w);
  alpha = mu/phi;
  beta = (1-mu)/phi;
}

model {
  target += normal_lpdf(w | 0, 1);
  target += beta_binomial_lpmf(y | trials, alpha, beta);
}

