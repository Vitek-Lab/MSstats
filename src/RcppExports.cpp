// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_estimable_fixed_random
List get_estimable_fixed_random(const List& parameters, const arma::vec& contrast);
RcppExport SEXP _MSstats_get_estimable_fixed_random(SEXP parametersSEXP, SEXP contrastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type contrast(contrastSEXP);
    rcpp_result_gen = Rcpp::wrap(get_estimable_fixed_random(parameters, contrast));
    return rcpp_result_gen;
END_RCPP
}
// make_contrast_run_quant
NumericVector make_contrast_run_quant(DataFrame input, NumericVector coefs, NumericVector contrast_matrix, NumericMatrix counts, bool is_labeled, bool is_reference);
RcppExport SEXP _MSstats_make_contrast_run_quant(SEXP inputSEXP, SEXP coefsSEXP, SEXP contrast_matrixSEXP, SEXP countsSEXP, SEXP is_labeledSEXP, SEXP is_referenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type input(inputSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type contrast_matrix(contrast_matrixSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< bool >::type is_labeled(is_labeledSEXP);
    Rcpp::traits::input_parameter< bool >::type is_reference(is_referenceSEXP);
    rcpp_result_gen = Rcpp::wrap(make_contrast_run_quant(input, coefs, contrast_matrix, counts, is_labeled, is_reference));
    return rcpp_result_gen;
END_RCPP
}
// get_linear_summary
NumericVector get_linear_summary(const DataFrame& input, const NumericVector& coefs, const NumericMatrix& counts, const bool is_labeled);
RcppExport SEXP _MSstats_get_linear_summary(SEXP inputSEXP, SEXP coefsSEXP, SEXP countsSEXP, SEXP is_labeledSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type input(inputSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_labeled(is_labeledSEXP);
    rcpp_result_gen = Rcpp::wrap(get_linear_summary(input, coefs, counts, is_labeled));
    return rcpp_result_gen;
END_RCPP
}
// median_polish_summary
NumericVector median_polish_summary(NumericMatrix x, double eps, int maxiter);
RcppExport SEXP _MSstats_median_polish_summary(SEXP xSEXP, SEXP epsSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(median_polish_summary(x, eps, maxiter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSstats_get_estimable_fixed_random", (DL_FUNC) &_MSstats_get_estimable_fixed_random, 2},
    {"_MSstats_make_contrast_run_quant", (DL_FUNC) &_MSstats_make_contrast_run_quant, 6},
    {"_MSstats_get_linear_summary", (DL_FUNC) &_MSstats_get_linear_summary, 4},
    {"_MSstats_median_polish_summary", (DL_FUNC) &_MSstats_median_polish_summary, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSstats(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
