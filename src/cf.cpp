#include <RcppArmadillo.h>



/** Evaluates the squared exponential covariance function between given input matrices.
*/
  // [[Rcpp::export]]
arma::mat cf_sexp_c(arma::mat x1, // input matrix 1
                    arma::mat x2, // input matrix 2
                    double lscale2, // length-scale squared
                    double magn2)  // magnitude squared
{
  size_t n1 = x1.n_rows;
  size_t n2 = x2.n_rows;
  arma::mat K(n1,n2);
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