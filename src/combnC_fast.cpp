/*
 * This code generate all combinations of n elements, taken m at a time
 * This is fast version of combn(), and was modified from source code of
 * package Rfast <https://CRAN.R-project.org/package=Rfast>
 */

#define ARMA_NO_DEBUG
#define ARMA_USE_BLAS

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rinternals.h>
#include <R.h>

using namespace Rcpp;
using namespace arma;

static void combn_mat(arma::vec& vals, const int n, const unsigned int start_idx,
                      std::vector<double>& combn_data, double*& combn_col) {
  if (!n) {
    for (unsigned int i = 0; i < combn_data.size(); ++i) {
      *combn_col++ = combn_data[i];
    }
    return;
  }
  for (unsigned int i = start_idx; i <= (vals.size() - n); ++i) {
    combn_data.at(combn_data.size() - n) = vals(i);
    combn_mat(vals, n - 1, i + 1, combn_data, combn_col);
  }
}

// [[Rcpp::export]]
double fast_combn_sum(arma::vec vals, int n, arma::vec n_k) {
  const unsigned int nrows = n;
  const unsigned int ncols = std::round(R::choose(vals.size(), n));
  std::vector<double> combn_data(nrows);
  const unsigned int start_idx = 0;
  NumericMatrix combn_ds;
  static double* combn_col;
  combn_ds = PROTECT(Rf_allocMatrix(REALSXP, nrows, ncols));
  combn_col = REAL(combn_ds);
  combn_mat(vals, n, start_idx, combn_data, combn_col);
  UNPROTECT(1);
  //
  double res = 0;
  for (unsigned int i = 0; i < ncols; ++i) {
    vec temp = combn_ds(_,i);
    res += n_k(temp[0] - 1) * n_k(temp[1] - 1) * n_k(temp[2] - 1);
  }
  return res;
}
