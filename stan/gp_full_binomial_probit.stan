
  
data {
  int<lower=0> n;
  int<lower=0> trials[n];
  matrix[n,n] K; 
  int<lower=0> y[n];
}

transformed data {
  vector[n] zeros;
  zeros = rep_vector(0, n);
}

parameters {
  vector[n] f;
}

model {
  target += multi_normal_lpdf(f | zeros, K);
  target += binomial_lpmf(y | trials, Phi_approx(f));
}
  
