

data {
  int<lower=0> n;
  int<lower=0> m; // number of inducing points
  int<lower=0> trials[n];
  matrix[n,m] Kxz_U_inv; // product of Kxz and inverse of upper cholesky of the Kuu
  vector[n] D; // FITC diagonal correction
  int<lower=0> y[n];
  int<lower=0> link;
}

parameters {
  vector[m] u_white;
  vector[n] e; 
}

transformed parameters {
  vector[n] f;
  vector[n] mu;
  f = Kxz_U_inv * u_white + sqrt(D) .* e;
  if (link == 0)
    mu = inv_logit(f);
  else if (link == 1)
    mu = Phi_approx(f);
}

model {
  target += normal_lpdf(u_white | 0, 1);
  target += normal_lpdf(e | 0, 1);
  target += binomial_lpmf(y | trials, mu);
}

