
  
data {
  int<lower=0> n;
  int<lower=0> trials[n];
  matrix[n,n] L; // cholesky of the covariance matrix (lower triangular) 
  int<lower=0> y[n];
  int<lower=0> link;
}

parameters {
  vector[n] f_white;
}

transformed parameters {
  vector[n] f;
  vector[n] mu;
  f = L*f_white;
  if (link == 0)
    mu = inv_logit(f);
  else if (link == 1)
    mu = Phi_approx(f);
}

model {
  target += normal_lpdf(f_white | 0, 1);
  target += binomial_lpmf(y | trials, mu);
}
  
