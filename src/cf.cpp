#include <RcppArmadillo.h>



/** Evaluates the squared exponential covariance function between given input matrices.
*/
  // [[Rcpp::export]]
arma::mat cf_sexp_c(arma::mat x1, // input matrix 1
                    arma::mat x2, // input matrix 2
                    double lscale, // length-scale
                    double magn)  // magnitude
{
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
arma::mat cf_matern32_c(arma::mat x1, // input matrix 1
                        arma::mat x2, // input matrix 2
                        double lscale, // length-scale
                        double magn)  // magnitude
{
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
arma::mat cf_matern52_c(arma::mat x1, // input matrix 1
                        arma::mat x2, // input matrix 2
                        double lscale, // length-scale
                        double magn)  // magnitude
{
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





