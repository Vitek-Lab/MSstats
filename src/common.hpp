#pragma once
#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
arma::mat mult(const arma::mat& A, const arma::mat& B);
NumericVector unlist(List l);
NumericVector grep(String pattern, CharacterVector character_vector);