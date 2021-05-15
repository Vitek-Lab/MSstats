#include <RcppArmadillo.h>
#include "common.h"
using namespace Rcpp;

// [[Rcpp::export]]
List get_estimable_fixed_random(const List& parameters, const arma::vec& contrast) {
    arma::mat cf = parameters["cf"];
    arma::mat coefs_m = cf.cols(0, 0);
    arma::mat contrast_m = reshape(contrast, 1, contrast.n_elem);

    arma::mat estimated = mult(contrast_m, coefs_m);
    NumericMatrix vcv = parameters["vcv"];
    arma::mat vcv_m = as<arma::mat>(vcv);
    arma::mat vc_square = mult(mult(contrast_m, vcv_m), trans(contrast_m));
    double vc = sqrt(vc_square(0, 0));
    int df =  parameters["df"];
    double tval = estimated(0, 0) / vc;
    double prob = 2 * (1 - R::pt(fabs(tval), df, 1, 0));

    List result = List::create(Named("logFC") = estimated(0, 0),
                               Named("SE") = vc,
                               Named("Tvalue") = tval,
                               Named("DF") = df,
                               Named("pvalue") = prob
    );
    return(result);
}
