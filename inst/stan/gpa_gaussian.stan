

data {
  int<lower=0> n;
  int<lower=0> m;
  matrix[n,m] Z; 
  vector[n] y;
  real<lower=0> sigma;
}

parameters {
  vector[m] w;
}

model {
  target += normal_lpdf(w | 0, 1);
  target += normal_lpdf(y | Z*w, sigma);
}


