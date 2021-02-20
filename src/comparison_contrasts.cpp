#include <RcppArmadillo.h>
#include "common.hpp"
using namespace Rcpp;

int get_index(LogicalVector x) {
    for (int i = 0; i < x.length(); ++i) {
        bool val = x[i];
        if (val) {
            return(i);
        } else {
            continue;
        }
    }
    return(-1);
}

LogicalVector find_matches(CharacterVector x, String y) {
    LogicalVector result(x.length());
    for (int i = 0; i < result.length(); ++i) {
        result[i] = (x[i] == y);
    }
    return(result);
}

// TODO: template
IntegerVector stl_sort(IntegerVector x) {
    IntegerVector y = clone(x);
    std::sort(y.begin(), y.end());
    return y;
}

NumericVector stl_sort(NumericVector x) {
    NumericVector y = clone(x);
    std::sort(y.begin(), y.end());
    return y;
}


NumericVector get_intercept_contrast(CharacterVector coef_names) {
    NumericVector intercept_name = grep("Intercept", coef_names);
    CharacterVector temp_intercept = coef_names[intercept_name];

    NumericVector intercept(temp_intercept.length());
    if (is_true(all(is_na(temp_intercept))) || (temp_intercept.length() == 0)) {
        intercept = NumericVector(0);
    } else {
        intercept.names() = temp_intercept;
    }
    return(intercept);
}


// [[Rcpp::export]]
NumericVector get_patient_seq(const DataFrame& patients, const NumericVector& patient_counts,
                              const CharacterVector& coef_names) {
    NumericVector find_subjects = grep("SUBJECT", coef_names);
    NumericVector find_nested = grep(":|NESTED", coef_names);
    find_subjects = setdiff(find_subjects, find_nested);
    bool any_positive = false;
    if (find_subjects.length() > 0) {
        find_subjects = stl_sort(find_subjects);
        LogicalVector is_positive = (find_subjects > 0);
        any_positive = sum(is_positive) > 0;
    }
    if (any_positive) {
        CharacterVector names = coef_names[find_subjects];
        bool is_nonzero = false;
        CharacterVector subjects = patients["SUBJECT"];
        CharacterVector groups = patients["GROUP"];
        NumericVector values = patients["Value"];
        IntegerVector match_subjects = match(subjects, coef_names);
        IntegerVector unique_matches = unique(match_subjects);
        unique_matches = unique_matches[!is_na(unique_matches)];
        unique_matches = stl_sort(unique_matches);
        NumericVector patient_seq(names.length());
        for (int i = 0; i < unique_matches.length(); ++i) {
            int index = unique_matches[i] - 1;
            double value = values[index];
            if (value != 0) {
                String group = groups[index];
                int count = patient_counts[group];
                double patient_value = value / count;
                patient_seq[index] += patient_value;
            }
        }
        patient_seq.attr("names") = names;
        return(patient_seq);
    } else {
        return(NumericVector(0));
    }
}


// [[Rcpp::export]]
NumericVector get_group(const IntegerVector groups,
                        const NumericMatrix& contrast_matrix,
                        const CharacterVector& coef_names) {
    NumericVector find_group = grep("GROUP", coef_names);
    NumericVector find_colon = grep(":", coef_names);
    find_group = setdiff(find_group, find_colon);
    find_group = stl_sort(find_group);
    if (find_group.length() > 0 && !(find_group[0] == -1)) {
        CharacterVector names = coef_names[find_group];
        NumericVector temp_contrast = contrast_matrix(0, _);
        temp_contrast = temp_contrast[groups - 1];
        NumericVector result = temp_contrast;
        result.erase(0);
        result.attr("names") = names;
        return(result);
    } else {
        return(NumericVector(0));
    }
}


// [[Rcpp::export]]
NumericVector get_interaction_seq(const DataFrame& patients,
                                  const NumericVector& patient_counts,
                                  const CharacterVector& coef_names) {
    NumericVector find_subjects = grep("SUBJECT", coef_names);
    NumericVector find_groups = grep("GROUP", coef_names);
    NumericVector find_interaction = intersect(find_subjects, find_groups);

    if (find_interaction.length() > 0 && !(find_groups[0] == -1)) {
        find_interaction = stl_sort(find_interaction);
        CharacterVector names = coef_names[find_interaction];
        CharacterVector groups = patients["GROUP"];
        CharacterVector subjects = patients["SUBJECT"];
        NumericVector interaction_seq(names.length());
        CharacterVector interaction_labels(patients.nrow());
        for (int i = 0; i < patients.nrow(); ++i) {
            String group = groups[i];
            String subject = subjects[i];
            group += String(":");
            group += subject;
            interaction_labels[i] = group;
        }
        NumericVector values = patients["Value"];
        bool is_matching = false;

        LogicalVector matches(names.length());
        for (int i = 0; i < values.length(); ++i) {
            String label = interaction_labels[i];
            matches = find_matches(names, label);
            is_matching = is_true(any(matches));
            bool is_nonzero = (values[i] != 0);
            if (is_matching && is_nonzero) {
                int index = get_index(matches);
                double value = values[i];
                String group = groups[index];
                int count = patient_counts[group];
                double patient_value = value / count;
                interaction_seq[index] += patient_value;
            }
        }
        interaction_seq.attr("names") = names;
        return(interaction_seq);
    } else {
        return(NumericVector(0));
    }
}

// [[Rcpp::export]]
NumericVector get_contrast_free(const DataFrame& input,
                                const DataFrame& patients,
                                const NumericVector& patient_counts,
                                const NumericMatrix& contrast_matrix,
                                const NumericVector& coefs) {
    CharacterVector coef_names = coefs.attr("names");
    IntegerVector groups = input["GROUP"];
    groups = unique(groups);
    groups = stl_sort(groups);

    NumericVector intercept = get_intercept_contrast(coef_names);
    NumericVector group = get_group(groups, contrast_matrix, coef_names);
    NumericVector subject = get_patient_seq(patients, patient_counts, coef_names);
    NumericVector interaction = get_interaction_seq(patients, patient_counts, coef_names);

    NumericVector contrast = unlist(List::create(intercept, group, subject, interaction));
    contrast = contrast[!is_na(coefs)];
    return(contrast);
}


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
