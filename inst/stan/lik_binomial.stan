
  
data {
  int<lower=0> n;
  int<lower=0> trials[n];
  int<lower=0> y[n];
  int<lower=0> link;
}

parameters {
  vector[n] f;
}

transformed parameters {
  vector[n] mu;
  if (link == 0)
    mu = inv_logit(f);
  else if (link == 1)
    mu = Phi_approx(f);
}

model {
  target += binomial_lpmf(y | trials, mu);
}
  
