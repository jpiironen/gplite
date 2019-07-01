
  
data {
  int<lower=0> n;
  int<lower=0> trials[n];
  matrix[n,n] L; // cholesky of the covariance matrix (lower triangular) 
  int<lower=0> y[n];
}

parameters {
  vector[n] f_white;
}

transformed parameters {
  vector[n] f;
  f = L*f_white;
}

model {
  target += normal_lpdf(f_white | 0, 1);
  target += binomial_lpmf(y | trials, Phi_approx(f));
}
  
