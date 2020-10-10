#' Handle censored missing values 
#' 
#' @param input `data.table` in MSstats data format
#' @param summary_method summarization method (`summaryMethod` parameter to `dataProcess`)
#' @param impute if TRUE, missing values are supposed to be imputed 
#' (`MBimpute` parameter to `dataProcess`)
#' @param missing_symbol `censoredInt` parameter to `dataProcess`
#' @param censored_cutoff `cutoffCensored` parameter to `dataProcess`
#' 
#' @export
#'  
MSstatsHandleMissing = function(input, summary_method, impute, 
                                missing_symbol, censored_cutoff) {
    INTENSITY = LABEL = ABUNDANCE = NULL
    
    if (summary_method == "TMP" & impute) {
        is_labeled = nlevels(input$LABEL) == 2 # label <- nlevels(input$LABEL) == 2
        input$censored = FALSE
        ## if intensity = 1, but abundance > cutoff after normalization, it also should be censored.
        if (!is.null(censored_cutoff)) {
            quantiles = input[!is.na(INTENSITY) & INTENSITY > 1 & LABEL == "L", 
                              quantile(ABUNDANCE, 
                                       prob = c(0.01, 0.25, 0.5, 0.75, 
                                                censored_cutoff), 
                                       na.rm = TRUE)]
            iqr = quantiles[4] - quantiles[2]
            multiplier = (quantiles[5] - quantiles[4]) / iqr
            cutoff_lower = (quantiles[2] - multiplier * iqr) 
            input$censored = !is.na(input$INTENSITY) & 
                input$LABEL == "L" &
                input$ABUNDANCE < cutoff_lower
            if (cutoff_lower <= 0 & !is.null(missing_symbol) & missing_symbol == "0") {
                # zero_one_filter = (!is.na(input$INTENSITY) & input$INTENSITY == 1) |
                #     (!is.na(input$ABUNDANCE & input$ABUNDANCE) <= 0)
                zero_one_filter = !is.na(input$ABUNDANCE & input$ABUNDANCE) <= 0
                input$censored = ifelse(zero_one_filter, TRUE, input$censored)
            }
            if (!is.null(missing_symbol) & missing_symbol == "NA") {
                input$censored = ifelse(is.na(input$INTENSITY), TRUE, 
                                        input$censored)
            }
            
            if (!is_labeled) {
                msg = paste('** Log2 intensities under cutoff =', 
                            format(cutoff_lower, digits = 5), 
                            ' were considered as censored missing values.')
                msg_2 = paste("** Log2 intensities =", missing_symbol, "were considered as censored missing values.")
            } else {
                msg = paste('** Log2 endogenous intensities under cutoff =', 
                            format(cutoff_lower, digits = 5), 
                            ' were considered as censored missing values.')
                getOption("MSstatsMsg")("INFO", msg)
                getOption("MSstatsMsg")("INFO", msg)
            }
            getOption("MSstatsMsg")("INFO", msg)
            getOption("MSstatsMsg")("INFO", msg_2)
        } else {
            if (missing_symbol == '0') {
                input$censored = input$LABEL == "L" & 
                    !is.na(input$INTENSITY) &
                    (input$INTENSITY == 1 | input$ABUNDANCE <= 0)
            } else if (missing_symbol == 'NA') {
                input$censored = input$LABEL == "L" & is.na(input$ABUNDANCE)
            }
        }
    }
    input
}


# .getMin = function(abundance, nonmissing) {
#     just_nonmissing = abundance[nonmissing]
#     if (length(just_nonmissing) > 0) {
#         0.99*min(, na.rm = TRUE)
#     } else {
#         NA
#     }
# }

.getMin = function(abundance, nonmissing) {
    0.99*min(abundance[nonmissing], na.rm = TRUE)
}

#' Set censored values based on minimum in run/feature/run or feature
#' @param input `data.table` in MSstats format
#' @param cutoff_base cutoffCensored parameter to `dataProcess`
#' @param censored_symbol censoredInt parameter to `dataProcess`
#' @param remove50missing if TRUE, features with at least 50% missing values
#' will be removed
#' @keywords internal
.setCensoredByThreshold = function(input, cutoff_base, censored_symbol,
                                   remove50missing) {
    ABUNDANCE = newABUNDANCE = perc_nm = RUN = FEATURE = NULL
    
    if (censored_symbol == "NA") {
        input[, nonmissing_all := !is.na(newABUNDANCE)]
    } else if (censored_symbol == "0") {
        input[, nonmissing_all := !is.na(newABUNDANCE) & input$newABUNDANCE != 0]
    }
    
    # remove feature with 1 measurement
    input[, nonmissing_all := ifelse(total_features > 1 & n_obs <= 1, FALSE, nonmissing_all)] 
    
    if (cutoff_base == "minFeature") {
        grouping_vars = c("PROTEIN", "FEATURE", "LABEL")
    } else if (cutoff_base == "minRun") {
        grouping_vars = c("PROTEIN", "RUN", "LABEL")
    } 
    
    if (cutoff_base %in% c("minFeature", "minRun")) {
        input[n_obs > 1 & n_obs_run > 0, 
              ABUNDANCE_cut := .getMin(newABUNDANCE, nonmissing_all),
              by = grouping_vars]
        input[, newABUNDANCE := ifelse(!nonmissing_all & censored, 
                                       ABUNDANCE_cut, newABUNDANCE)]
    } else {
        feature_cutoffs = input[nonmissing_filter, 
                                list(ABUNDANCE_cut_fea = 0.99*min(newABUNDANCE)),
                                by = c("PROTEIN", "FEATURE", "LABEL")]
        
        if (remove50missing & censored_symbol == "0") {
            n_features = data.table::uniqueN(input$FEATURE)
            missing_runs = input[, list(perc_nm = .N / n_features), 
                                 by = c("PROTEIN", "RUN")]
            missing_runs = missing_runs[perc_nm < 0.5, ]
            if (nrow(missing_runs) > 0) {
                input = input[!(RUN %in% unique(missing_runs$RUN)), ]
            }
        }
        
        run_cutoffs = input[nonmissing_filter, 
                            list(ABUNDANCE_cut_run = 0.99*min(newABUNDANCE)),
                            by = c("PROTEIN", "RUN", "LABEL")]
        cutoffs = merge(feature_cutoffs, run_cutoffs,
                        by = c("PROTEIN", "LABEL"),
                        allow.cartesian = TRUE, sort = FALSE)
        cutoffs$final_cutoff = sapply(1:nrow(cutoffs), 
                                      function(i) min(cutoffs$ABUNDANCE_cut_fea[i],
                                                      cutoffs$ABUNDANCE_cut_run[i]))
        # input[, n_feat_run := data.table::uniqueN(.SD), by = "PROTEIN",
        #       .SDcols = c("FEATURE", "RUN")]
        # input[, ABUNDANCE := ifelse(n_feat_run > 1,
        #                             ifelse(nonmissing, ABUNDANCE, ))]
        # if (data.table::uniqueN(input[, list(FEATURE, RUN)]) > 1) {
        #     input$ABUNDANCE = ifelse(nonmissing_filter,
        #                              input$ABUNDANCE, cutoffs$final_cutoff)
        # } else {
        #     if (censored_symbol == "0") {
        #         input$ABUNDANCE = ifelse(nonmissing_filter, input$ABUNDANCE, 
        #                                  cutoffs$ABUNDANCE_cut_fea)
        #     }
        # }
    }
}


#' Identify non-missing values
#' @param input `data.table` in MSstats format
#' @param impute if TRUE, missing values are supposed to be imputed
#' @param censored_symbol `censoredInt` parameter to dataProcess
#' @keywords internal
.getNonMissingFilter = function(input, impute, censored_symbol) {
    if (impute) {
        if (!is.null(censored_symbol)) {
            if (censored_symbol == "0") {
                nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE) & input$newABUNDANCE != 0
            } else if (censored_symbol == "NA") {
                nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE)
            }  
        } 
    } else {
        nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE) & input$newABUNDANCE != 0
    }
    nonmissing_filter
}
