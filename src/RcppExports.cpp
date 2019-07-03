// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// cf_sexp_c
arma::mat cf_sexp_c(arma::mat x1, arma::mat x2, double lscale, double magn);
RcppExport SEXP _gplite_cf_sexp_c(SEXP x1SEXP, SEXP x2SEXP, SEXP lscaleSEXP, SEXP magnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type lscale(lscaleSEXP);
    Rcpp::traits::input_parameter< double >::type magn(magnSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_sexp_c(x1, x2, lscale, magn));
    return rcpp_result_gen;
END_RCPP
}
// cf_matern32_c
arma::mat cf_matern32_c(arma::mat x1, arma::mat x2, double lscale, double magn);
RcppExport SEXP _gplite_cf_matern32_c(SEXP x1SEXP, SEXP x2SEXP, SEXP lscaleSEXP, SEXP magnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type lscale(lscaleSEXP);
    Rcpp::traits::input_parameter< double >::type magn(magnSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_matern32_c(x1, x2, lscale, magn));
    return rcpp_result_gen;
END_RCPP
}
// cf_matern52_c
arma::mat cf_matern52_c(arma::mat x1, arma::mat x2, double lscale, double magn);
RcppExport SEXP _gplite_cf_matern52_c(SEXP x1SEXP, SEXP x2SEXP, SEXP lscaleSEXP, SEXP magnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type lscale(lscaleSEXP);
    Rcpp::traits::input_parameter< double >::type magn(magnSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_matern52_c(x1, x2, lscale, magn));
    return rcpp_result_gen;
END_RCPP
}
// cf_nn_c
arma::mat cf_nn_c(arma::mat x1, arma::mat x2, double sigma0, double sigma, double magn);
RcppExport SEXP _gplite_cf_nn_c(SEXP x1SEXP, SEXP x2SEXP, SEXP sigma0SEXP, SEXP sigmaSEXP, SEXP magnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type magn(magnSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_nn_c(x1, x2, sigma0, sigma, magn));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4gp_binomial_logit_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gp_binomial_probit_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gp_gaussian_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gpa_binomial_logit_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gpa_binomial_probit_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gpa_gaussian_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_gplite_cf_sexp_c", (DL_FUNC) &_gplite_cf_sexp_c, 4},
    {"_gplite_cf_matern32_c", (DL_FUNC) &_gplite_cf_matern32_c, 4},
    {"_gplite_cf_matern52_c", (DL_FUNC) &_gplite_cf_matern52_c, 4},
    {"_gplite_cf_nn_c", (DL_FUNC) &_gplite_cf_nn_c, 5},
    {"_rcpp_module_boot_stan_fit4gp_binomial_logit_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gp_binomial_logit_mod, 0},
    {"_rcpp_module_boot_stan_fit4gp_binomial_probit_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gp_binomial_probit_mod, 0},
    {"_rcpp_module_boot_stan_fit4gp_gaussian_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gp_gaussian_mod, 0},
    {"_rcpp_module_boot_stan_fit4gpa_binomial_logit_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gpa_binomial_logit_mod, 0},
    {"_rcpp_module_boot_stan_fit4gpa_binomial_probit_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gpa_binomial_probit_mod, 0},
    {"_rcpp_module_boot_stan_fit4gpa_gaussian_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gpa_gaussian_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_gplite(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
