#include <RcppArmadillo.h>
using namespace Rcpp;


//[[Rcpp::depends(RcppArmadillo)]]
arma::mat mult(const arma::mat& A, const arma::mat& B) {
    return A * B;
}


NumericVector unlist(List l) {
    Rcpp::Environment base("package:base");
    Rcpp::Function unlist = base["unlist"];
    return(unlist(l));
}


NumericVector grep(String pattern, CharacterVector character_vector) {
    Rcpp::Environment base("package:base");
    Rcpp::Function grep = base["grep"];
    Rcpp::NumericVector result = grep(Rcpp::_["pattern"] = pattern,
                                      Rcpp::_["x"]  = character_vector);
    if (result.length() > 0 && !is_true(all(is_na(result)))) {
        result = result - 1;
    } else {
        result = -1;
    }
    return result;
}
