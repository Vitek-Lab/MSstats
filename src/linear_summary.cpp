#include <RcppArmadillo.h>
#include "common.h"
using namespace Rcpp;


NumericVector get_feature_props(DataFrame input, NumericVector contrast_matrix,
                                NumericMatrix counts) {
    LogicalVector logical_run = contrast_matrix == 1;
    int select_run = which_max(NumericVector(logical_run)); 
    NumericVector run_counts = counts(select_run, _);
    NumericVector result = run_counts[seq(1, run_counts.length() - 1)] / sum(run_counts);
    return(result);
} 


NumericVector get_intercept(CharacterVector coef_names) {
    NumericVector intercept_name = grep("Intercept", coef_names);
    CharacterVector temp_intercept = coef_names[intercept_name];

    NumericVector intercept(temp_intercept.length(), 1.0);
    if (is_true(all(is_na(temp_intercept))) || (temp_intercept.length() == 0)) {
        intercept = NumericVector(0.0);
    } else { 
        intercept.names() = temp_intercept;
    } 
    return(intercept);
} 

NumericVector get_features(CharacterVector coef_names, NumericVector find_features,
                           NumericVector find_colon, NumericVector contrast_matrix,
                           NumericMatrix counts, DataFrame input) {
    NumericVector just_features = setdiff(find_features, find_colon);
    CharacterVector temp_feature = coef_names[just_features];
    
    NumericVector feature(0);
    if (temp_feature.length() != 0) {
        feature = get_feature_props(input, contrast_matrix, counts);
        feature.attr("names") = temp_feature;
    } 
    return(feature);
} 

NumericVector get_run(CharacterVector coef_names,
                      NumericVector find_runs, NumericVector find_colon,
                      NumericVector contrast_matrix, bool label, int n_runs) {
    NumericVector just_runs = setdiff(find_runs, find_colon);
    NumericVector run(0);
    if (just_runs.length() != 0) {
        CharacterVector temp_run = coef_names[just_runs];
        if (label) {
            run = rep(1 / n_runs, temp_run.length());
        } else {
            run = contrast_matrix[seq(1, contrast_matrix.length() - 1)];
        }
        run.attr("names") = temp_run;
    }
    return(run);
}


NumericVector get_ref(const CharacterVector& coef_names, const NumericVector& find_ref, 
                      const NumericVector& contrast_matrix, const DataFrame& input,
                      const bool is_reference) {
    NumericVector ref(0);
    if ((find_ref.length() != 0) & !(find_ref[0] == -1)) {
        if (is_reference) {
            CharacterVector temp_ref = coef_names[find_ref];
            ref = rep(0.0, find_ref.length());
            ref.attr("names") = temp_ref;
        } else {
            CharacterVector temp_ref = coef_names[find_ref];
            CharacterVector refs = input["ref"];
            CharacterVector unique_refs = unique(refs);
            int n_refs = refs.length();
            ref = contrast_matrix[seq(0, n_refs - 2)];
            ref.attr("names") = temp_ref;
        }
    }
    return(ref);
}


NumericVector get_feature_run(NumericVector find_runs, NumericVector find_features,
                              CharacterVector coef_names, NumericMatrix& counts) {
    NumericVector find_run_feature = intersect(find_runs, find_features);
    int n_rows = counts.nrow();
    NumericVector rf(0);
    if ((find_run_feature.length() != 0) & !(find_run_feature[0] == -1)) {
        CharacterVector temp_rf = coef_names[find_run_feature];
        rf = rep(1 / n_rows, temp_rf.length());
        rf.attr("names") = temp_rf;
    }
    return(rf);
}


double get_quant(const NumericVector& coefs, const NumericVector& contrast) {
    NumericMatrix coefs_tmp(coefs.length(), 1, coefs.begin());
    NumericMatrix contrast_tmp(1, contrast.length(), contrast.begin());
    arma::mat contrast_m = as<arma::mat>(contrast_tmp);
    arma::mat coefs_m = as<arma::mat>(coefs_tmp);
    arma::mat result = mult(contrast_m, coefs_m);
    return (result(0, 0));
}


NumericVector combine_contrast(bool is_reference, NumericVector intercept, 
                               NumericVector features, NumericVector runs, 
                               NumericVector refs, NumericVector run_features) {
    if (is_reference) {
        NumericVector contrast = unlist(List::create(intercept, features, runs, 
                                                     refs, run_features));
        return(contrast);
    } else {
        NumericVector contrast = unlist(List::create(intercept, features, runs, 
                                                     refs));
        return(contrast);
    }
}

// [[Rcpp::export]]
NumericVector make_contrast_run_quant(DataFrame input,
                                      NumericVector coefs,
                                      NumericVector contrast_matrix,
                                      NumericMatrix counts,
                                      bool is_labeled,
                                      bool is_reference = false) {
    CharacterVector all_runs = input["RUN"];
    CharacterVector unique_runs = unique(all_runs);
    int n_runs = unique_runs.length();
    
    CharacterVector coef_names = coefs.attr("names");
    NumericVector find_features = grep("FEATURE", coef_names);
    NumericVector find_colon = grep(":", coef_names);
    NumericVector find_runs = grep("RUN", coef_names);
    NumericVector find_ref = grep("ref", coef_names);

    NumericVector intercept = get_intercept(coef_names);
    NumericVector features = get_features(coef_names, find_features, find_colon,
                                          contrast_matrix, counts, input);
    NumericVector runs = get_run(coef_names, find_runs, find_colon, 
                                 contrast_matrix, is_labeled, n_runs);
    NumericVector refs = get_ref(coef_names, find_ref, contrast_matrix, input,
                                 is_reference);
    NumericVector run_features = get_feature_run(find_runs, find_features,
                                                 coef_names, counts);

    NumericVector contrast = combine_contrast(is_reference,
                                              intercept, features, runs,
                                              refs, run_features);
    if (contrast.length() == coefs.length()) {
        contrast = contrast[!is_na(coefs)];
    } else {
        CharacterVector names = contrast.attr("names");
        contrast = contrast[!is_na(names)];
    }
    return(contrast);
}


// [[Rcpp::export]]
NumericVector get_linear_summary(const DataFrame& input,
                                 const NumericVector& coefs,
                                 const NumericMatrix& counts,
                                 const bool is_labeled) {
    CharacterVector runs = input["RUN"];
    CharacterVector unique_runs = unique(runs);
    CharacterVector ref = {"ref"};
    unique_runs = setdiff(unique_runs, ref);
    int num_runs = unique_runs.length();
    NumericVector log_intensities(num_runs);
    
    for (int run_id = 0; run_id < num_runs; ++run_id) {
        NumericVector contrast_matrix(num_runs);
        contrast_matrix[run_id] = 1;
        NumericVector contrast = make_contrast_run_quant(input,
                                                         coefs,
                                                         contrast_matrix,
                                                         counts,
                                                         is_labeled);
        double quantified = get_quant(coefs, contrast);
        log_intensities(run_id) = quantified;
    }
    
    if (is_labeled) {
        NumericVector contrast_matrix(num_runs);
        contrast_matrix[num_runs - 1] = 1;
        NumericVector contrast = make_contrast_run_quant(input, coefs, 
                                                         contrast_matrix, 
                                                         counts, true, true);
        log_intensities.push_back(get_quant(coefs, contrast));
    }
    
    return(log_intensities);
}
