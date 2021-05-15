#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector median_polish_summary(NumericMatrix x, double eps = 0.01,
                                    int maxiter = 10L) {
    int num_cols = x.ncol();
    int num_rows = x.nrow();
    NumericVector current_col_medians(num_cols);
    NumericVector current_row_medians(num_rows);
    NumericMatrix y = clone(x);
    double old_sum = 0.0;
    double overall_effect = 0.0;
    NumericVector row_effects = NumericVector(num_rows);
    NumericVector column_effects = NumericVector(num_cols);
    
    double overall_change = 0.0;
    double row_change = 0.0;
    double new_sum = 0.0;

    for (int iter = 0; iter < maxiter; ++iter) {
        for (int row_id = 0; row_id < num_rows; ++row_id) {
            current_row_medians[row_id] = median(NumericVector(y(row_id, _)), true);
        }
        for (int row_id = 0; row_id < num_rows; ++row_id) {
            y(row_id, _) = y(row_id, _) - current_row_medians[row_id];
        }
        row_effects = row_effects + current_row_medians;
        
        overall_change = median(NumericVector(column_effects), true);
        column_effects = column_effects - overall_change;
        overall_effect = overall_effect + overall_change;
        
        for (int col_id = 0; col_id < num_cols; ++col_id) {
            current_col_medians[col_id] = median(NumericVector(y(_, col_id)),
                                                 true);
        }
        for (int col_id = 0; col_id < num_cols; ++col_id) {
            y(_, col_id) = y(_, col_id) - current_col_medians[col_id];
        }
        column_effects = column_effects + current_col_medians;
        
        row_change = median(NumericVector(row_effects), true);
        row_effects = row_effects - row_change;
        
        overall_effect = overall_effect + row_change;
        new_sum = sum(abs(na_omit(y)));
        
        double diff = new_sum - old_sum;
        diff = diff < 0 ? -1 * diff : diff;
        if (new_sum == 0.0 || diff < eps * new_sum) {
            break;
        }
        old_sum = new_sum;
    }
    NumericVector result = overall_effect + row_effects;
    return(result);
}
