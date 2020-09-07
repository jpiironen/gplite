#include <RcppArmadillo.h>
#include <cmath>



/** Evaluates the squared exponential covariance function between given input matrices.
*/
  // [[Rcpp::export]]
arma::mat cf_sexp_c(
  arma::mat x1, // input matrix 1
  arma::mat x2, // input matrix 2
  double lscale, // length-scale
  double magn, // magnitude
  bool diag_only=false // whether to evaluate the diagonal only
) 
{
  if (diag_only) {
    size_t n = std::min(x1.n_rows, x2.n_rows);
    arma::mat K(n,1);
    size_t i;
    for (i=0; i<n; i++) {
      K.row(i) = cf_sexp_c(x1.row(i), x2.row(i), lscale, magn, false);
    }
    return(K);
  }
  
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double lscale2 = lscale*lscale;
  double magn2 = magn*magn;
  double r2;
  size_t i,j;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      r2 = sum(square(x1.row(i)-x2.row(j)));
      K(i,j) = magn2*exp(-0.5*r2/lscale2);
    }
  }
  return(K);
}


/** Evaluates the Matern nu=3/2 covariance function between given input matrices.
 */
// [[Rcpp::export]]
arma::mat cf_matern32_c(
  arma::mat x1, // input matrix 1
  arma::mat x2, // input matrix 2
  double lscale, // length-scale
  double magn, // magnitude
  bool diag_only=false // whether to evaluate the diagonal only
)
{
  if (diag_only) {
    size_t n = std::min(x1.n_rows, x2.n_rows);
    arma::mat K(n,1);
    size_t i;
    for (i=0; i<n; i++) {
      K.row(i) = cf_matern32_c(x1.row(i), x2.row(i), lscale, magn, false);
    }
    return(K);
  }
  
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double r;
  double magn2 = magn*magn;
  size_t i,j;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      r = sqrt(sum(square(x1.row(i)-x2.row(j))));
      K(i,j) = magn2*(1 + sqrt(3)*r/lscale) * exp(-sqrt(3)*r/lscale);
    }
  }
  return(K);
}


/** Evaluates the Matern nu=5/2 covariance function between given input matrices.
 */
// [[Rcpp::export]]
arma::mat cf_matern52_c(
  arma::mat x1, // input matrix 1
  arma::mat x2, // input matrix 2
  double lscale, // length-scale
  double magn, // magnitude
  bool diag_only=false // whether to evaluate the diagonal only
)
{
  if (diag_only) {
    size_t n = std::min(x1.n_rows, x2.n_rows);
    arma::mat K(n,1);
    size_t i;
    for (i=0; i<n; i++) {
      K.row(i) = cf_matern52_c(x1.row(i), x2.row(i), lscale, magn, false);
    }
    return(K);
  }
  
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double r;
  double lscale2 = lscale*lscale;
  double magn2 = magn*magn;
  size_t i,j;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      r = sqrt(sum(square(x1.row(i)-x2.row(j))));
      K(i,j) = magn2*(1 + sqrt(5)*r/lscale + 5*r*r/(3*lscale2)) * exp(-sqrt(5)*r/lscale);
    }
  }
  return(K);
}


/** Evaluates the neural network covariance function between given input matrices.
 */
// [[Rcpp::export]]
arma::mat cf_nn_c(
    arma::mat x1, // input matrix 1
    arma::mat x2, // input matrix 2
    double sigma0, // bias std
    double sigma, // weight std
    double magn, // magnitude
    bool diag_only = false // whether to evaluate diagonal only
)  
{
  
  if (diag_only) {
    size_t n = std::min(x1.n_rows, x2.n_rows);
    arma::mat K(n,1);
    size_t i;
    for (i=0; i<n; i++) {
      K.row(i) = cf_nn_c(x1.row(i), x2.row(i), sigma0, sigma, magn, false);
    }
    return(K);
  }
  
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
  double magn2 = magn*magn;
  double sigma0_sq = sigma0*sigma0*sigma*sigma;
  double a,b,c;
  size_t i,j;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      a = sigma0_sq + sum((x1.row(i)*sigma) % (x2.row(j)*sigma));
      b = sigma0_sq + sum(square(x1.row(i)*sigma));
      c = sigma0_sq + sum(square(x2.row(j)*sigma));
      K(i,j) = magn2*(2.0/M_PI)*asin( 2*a/sqrt((1+2*b)*(1+2*c)) );
    }
  }
  return(K);
}





