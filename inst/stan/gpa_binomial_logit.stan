

data {
  int<lower=0> n;
  int<lower=0> m;
  int<lower=0> trials[n];
  matrix[n,m] Z; 
  int<lower=0> y[n];
}

parameters {
  vector[m] w;
}

model {
  target += normal_lpdf(w | 0, 1);
  target += binomial_logit_lpmf(y | trials, Z*w);
}

