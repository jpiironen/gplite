

data {
  int<lower=0> n;
  matrix[n,n] K; 
  vector[n] y;
  real<lower=0> sigma; // noise std
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
  target += normal_lpdf(y | f, sigma);
}

