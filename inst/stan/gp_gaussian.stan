

data {
  int<lower=0> n;
  matrix[n,n] L; // cholesky of the covariance matrix (lower triangular) 
  vector[n] y;
  real<lower=0> sigma;
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
  target += normal_lpdf(y | f, sigma);
}

