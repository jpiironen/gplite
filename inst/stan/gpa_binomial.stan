

data {
  int<lower=0> n;
  int<lower=0> m;
  int<lower=0> trials[n];
  matrix[n,m] Z; 
  int<lower=0> y[n];
  int<lower=0> link;
}

parameters {
  vector[m] w;
}

transformed parameters {
  vector[n] mu;
  if (link == 0)
    mu = inv_logit(Z*w);
  else if (link == 1)
    mu = Phi_approx(Z*w);
}

model {
  target += normal_lpdf(w | 0, 1);
  target += binomial_lpmf(y | trials, mu);
}

