// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_combn_sum
double fast_combn_sum(arma::vec vals, int n, arma::vec n_k);
RcppExport SEXP _ClusROC_fast_combn_sum(SEXP valsSEXP, SEXP nSEXP, SEXP n_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n_k(n_kSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_combn_sum(vals, n, n_k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ClusROC_fast_combn_sum", (DL_FUNC) &_ClusROC_fast_combn_sum, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ClusROC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
