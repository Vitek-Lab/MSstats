#' Check if data represents repeated measurements design
#' 
#' @param summarization_output output of the dataProcess function
#' 
#' @return logical, TRUE if data represent repeated measurements design
#' 
#' @details This extracts information required by the group comparison workflow
#' 
#' @export
#' 
#' @examples
#' QuantData1 <- dataProcess(SRMRawData, use_log_file = FALSE)
#' checkRepeatedDesign(QuantData1)
#' 
checkRepeatedDesign = function(summarization_output) {
    SUBJECT = GROUP = NULL
    
    input = as.data.table(summarization_output$ProteinLevelData)
    subject_by_group = table(input[, list(SUBJECT, GROUP)])
    subject_appearances = apply(subject_by_group, 1, function(x) sum(x > 
                                                                         0))
    repeated = any(subject_appearances > 1)
    if (repeated) {
        msg = "Time course design of experiment - okay"
    }
    else {
        msg = "Case control design of experiment - okay"
    }
    getOption("MSstatsLog")("INFO", msg)
    repeated
}

#' Get information about number of measurements for each group
#' 
#' @param summarization_output output of the dataProcess function
#' 
#' @return data.table
#' 
#' @details This function extracts information required to compute percentages
#' of missing and imputed values in group comparison.
#' 
#' @export
#' 
#' @examples
#' QuantData <- dataProcess(DDARawData, use_log_file = FALSE)
#' samples_info <- getSamplesInfo(QuantData)
#' samples_info
#' 
getSamplesInfo = function(summarization_output) {
    RUN = NULL
    
    summarized = as.data.table(summarization_output$ProteinLevelData)
    summarized[, list(NumRuns = data.table::uniqueN(RUN)),
               by = "GROUP"]
}


#' Prepare data for a single protein for group comparison
#' @param single_protein data.table
#' @keywords internal
.prepareSingleProteinForGC = function(single_protein) {
    ABUNDANCE = GROUP = SUBJECT = RUN = NULL
    
    data.table::setnames(single_protein,
                         c("LogIntensities"),
                         c("ABUNDANCE"),
                         skip_absent = TRUE)
    single_protein = single_protein[!is.na(ABUNDANCE)]
    single_protein[, GROUP := factor(GROUP)]
    single_protein[, SUBJECT := factor(SUBJECT)]
    single_protein[, RUN := factor(RUN)]
    single_protein
}


#' Fit model and perform group comparison for a single protein
#' @param input data.table of summarized data
#' @param contrast_matrix contrast matrix
#' @param has_tech_replicates if TRUE, there are technical replicates
#' @param is_single_subject if TRUE, experiment consists of a single subject
#' @param repeated if TRUE, experiment consists of repeated measurements
#' @param groups unique labels for experimental conditions
#' @param samples_info number of runs per group
#' @param save_fitted_models if TRUE, fitted model will be saved. If FALSE,
#' it will be replaced by NULL
#' @param has_imputed if TRUE, missing values have been imputed by dataProcess
#' @importFrom stats resid fitted
#' @keywords internal
.fitModelSingleProtein = function(input, contrast_matrix, has_tech_replicates,
                                  is_single_subject, repeated, groups,
                                  samples_info,
                                  save_fitted_models, has_imputed) {
    GROUP = SUBJECT = NULL
    
    input[, GROUP := factor(GROUP)]
    input[, SUBJECT := factor(SUBJECT)]
    
    protein = unique(input$Protein)
    n_groups = nlevels(input$GROUP)
    if (n_groups == 1) {
        result = .getEmptyComparison(input, contrast_matrix,
                                     groups, protein)
        residuals = NA
        fitted_values = NA
        fit = NULL
    } else {
        fitted_model = .fitModelForGroupComparison(input, repeated,
                                                   is_single_subject,
                                                   has_tech_replicates)
        result = .getAllComparisons(input, fitted_model, contrast_matrix,
                                    groups, protein)
        result = .countMissingPercentage(contrast_matrix,
                                         input, result, samples_info,
                                         has_imputed)
        
        if (inherits(fitted_model[["full_fit"]], "lm")) {
            residuals = fitted_model[["full_fit"]][["residuals"]]
            fitted_values = fitted_model[["full_fit"]][["fitted.values"]]
        } else {
            residuals = resid(fitted_model[["full_fit"]])
            fitted_values = fitted(fitted_model[["full_fit"]])
        }
        fit = fitted_model[["full_fit"]]
    }
    
    if (!save_fitted_models) {
        fit = NULL
    }
    
    input[, residuals := residuals]
    input[, fitted := fitted_values]
    list(result, fit)
}


#' Choose a model type (fixed/mixed effects) and fit it for a single protein
#' @inheritParams .fitModelSingleProtein
#' @keywords internal
.fitModelForGroupComparison = function(input, repeated, is_single_subject,
                                       has_tech_replicates) {
    if (!repeated) {
        if (!has_tech_replicates | is_single_subject) {
            full_fit = lm(ABUNDANCE ~ GROUP, data = input)
            df_full = full_fit[["df.residual"]]
        } else {
            full_fit = suppressMessages(try(
                lme4::lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data = input),
                TRUE
            ))
            df_full = suppressMessages(try(
                lm(ABUNDANCE ~ GROUP + SUBJECT, data = input)$df.residual,
                TRUE
            ))
        }
    } else {
        ## time-course
        if (is_single_subject) {
            full_fit = lm(ABUNDANCE ~ GROUP,
                          data = input)
            df_full = full_fit$df.residual
        } else {
            ## no single subject
            if (!has_tech_replicates) {
                full_fit = suppressMessages(try(
                    lme4::lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data = input),
                    TRUE
                ))
                df_full = suppressMessages(try(
                    lm(ABUNDANCE ~ GROUP + SUBJECT, data = input)$df.residual,
                    TRUE))
            } else {
                full_fit = suppressMessages(try(
                    lme4::lmer(ABUNDANCE ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT),
                               data = input),
                    TRUE
                ))
                df_full = suppressMessages(try(
                    lm(ABUNDANCE ~ GROUP + SUBJECT + GROUP:SUBJECT,
                       data = input)$df.residual,
                    TRUE
                ))
            }
        }
    }
    list(full_fit = full_fit,
         df_full = df_full)
}


#' Get params (coefficients, covariance matrix, degrees of freedom) from a model
#' @param fitted_model object of class lm or lmerMod
#' @importFrom stats vcov
#' @importFrom methods is
#' @keywords internal
.getModelParameters = function(fitted_model) {
    if (is(fitted_model[["full_fit"]], "lm")) {
        model_summary = summary(fitted_model[["full_fit"]])
        cf = model_summary[["coefficients"]]
        vcv = model_summary[["cov.unscaled"]] * (model_summary[["sigma"]] ^ 2)
    } else {
        cf = as.matrix(lme4::fixef(fitted_model[["full_fit"]]))
        vcv = as.matrix(vcov(fitted_model[["full_fit"]]))
    }
    list(cf = cf, vcv = vcv, df = fitted_model[["df_full"]])
}


#' Comparison output when there are measurements only in a single condition
#' @param input summarized data
#' @param contrast_matrix contrast matrix
#' @param groups unique labels of experimental conditions
#' @param protein name of a protein
#' @keywords internal
.getEmptyComparison = function(input, contrast_matrix, groups, protein) {
    all_comparisons = lapply(seq_len(nrow(contrast_matrix)), function(row_id) {
        ith_comparison = contrast_matrix[row_id, , drop = FALSE]
        
        if (any(groups[ith_comparison != 0] %in% unique(input$GROUP))) {
            msg = paste("*** error: results of protein", protein,
                        "for comparison", row.names(ith_comparison),
                        "are NA because there are measurements",
                        "only in a single group")
            getOption("MSstatsLog")("INFO", msg)
            
            if (ith_comparison[ith_comparison != 0 &
                               (groups %in% unique(input$GROUP))] > 0) {
                list(
                    logFC = Inf,
                    issue = "oneConditionMissing"
                )
            } else {
                list(
                    logFC = -Inf,
                    issue = "oneConditionMissing"
                )
            }
        } else {
            msg = paste("*** error: results of protein", protein,
                        "for comparison", row.names(ith_comparison),
                        "are NA because there are no measurements",
                        "in both conditions.")
            getOption("MSstatsLog")("INFO", msg)
            
            list(logFC = NA,
                 issue = "completeMissing")
        }
    })
    empty_result = data.table::rbindlist(all_comparisons, fill = TRUE)
    empty_result = cbind(empty_result,
                         data.table::data.table(
                             Protein = protein, 
                             Label = row.names(contrast_matrix),
                             SE = NA, Tvalue = NA,
                             DF = NA, pvalue = NA
                         ))
    empty_result
}


#' Get all comparisons for a single protein and a contrast matrix
#' @param input summarized data
#' @param fitted_model model fitted by the .fitModelForGroupComparison function
#' @param contrast_matrix contrast matrix
#' @param groups unique labels of experimental conditions
#' @param protein name of a protein
#' @keywords internal
.getAllComparisons = function(input, fitted_model, contrast_matrix,
                              groups, protein) {
    empty_conditions = setdiff(groups, unique(input$GROUP))
    parameters = .getModelParameters(fitted_model)
    fit = fitted_model[["full_fit"]]
    coefs = parameters$cf[, 1]
    
    all_comparisons = vector("list", nrow(contrast_matrix))
    for (row_id in seq_len(nrow(contrast_matrix))) {
        ith_contrast = contrast_matrix[row_id, , drop = FALSE]
        if (length(empty_conditions) != 0) {
            result = .handleEmptyConditions(input, fit, ith_contrast,
                                            groups, parameters, protein,
                                            empty_conditions, coefs)
        } else {
            result = .handleSingleContrast(input, fit, ith_contrast, groups,
                                           parameters, protein, coefs)
        }
        all_comparisons[[row_id]] = result
    }
    data.table::rbindlist(all_comparisons, fill = TRUE)
}


#' Handle contrast when some of the conditions are missing
#' @param input summarized data
#' @param contrast single row of a contrast matrix
#' @param groups unique labels of experimental conditions
#' @param parameters parameters extracted from the model 
#' @param protein name of a protein
#' @param empty_conditions labels of empty conditions
#' @param coefs coefficient of the fitted model
#' @keywords internal
.handleEmptyConditions = function(input, fit, contrast,
                                  groups, parameters, protein,
                                  empty_conditions, coefs) {
    count_diff_pos = intersect(colnames(contrast)[contrast != 0 & contrast > 0],
                               empty_conditions)
    count_diff_neg = intersect(colnames(contrast)[contrast != 0 & contrast < 0],
                               empty_conditions)
    flag_issue_pos = length(count_diff_pos) != 0
    flag_issue_neg = length(count_diff_neg) != 0
    
    if (any(c(flag_issue_pos, flag_issue_neg))) {
        if (flag_issue_pos & flag_issue_neg) {
            issue = "completeMissing"
            logFC = NA
        } else if (flag_issue_pos & !flag_issue_neg) {
            issue_side = count_diff_pos
            logFC = -Inf
            issue = "oneConditionMissing"
        } else if (!flag_issue_pos & flag_issue_neg) {
            issue_side = count_diff_neg
            logFC = Inf
            issue = "oneConditionMissing"
        }
        result = list(Protein = protein,
                      logFC = logFC,
                      Label = row.names(contrast),
                      SE = NA, Tvalue = NA, DF = NA, pvalue = NA,
                      issue = issue)
    } else {
        result = .handleSingleContrast(input, fit, contrast, groups,
                                       parameters, protein, coefs)
    }
    result
}


#' Group comparison for a single contrast
#' @inheritParams .handleEmptyConditions
#' @keywords internal
.handleSingleContrast = function(input, fit, contrast, groups,
                                 parameters, protein, coefs) {
    groups = sort(groups)
    contrast_values = .getContrast(input, contrast, coefs, groups)
    parameters$cf = parameters$cf[names(contrast_values), , drop = FALSE]
    parameters$vcv = parameters$vcv[names(contrast_values), names(contrast_values)]
    result = get_estimable_fixed_random(parameters, contrast_values)
    if (is.null(result)) {
        result = list(Protein = protein,
                      Label = row.names(contrast),
                      logFC = NA, SE = NA, Tvalue = NA,
                      DF = NA, pvalue = NA, issue = NA)
    } else {
        result$Protein = protein
        result$Label = row.names(contrast)
        result$issue = NA
    }
    result
}


#' Create a contrast for a model with only group as a fixed effect
#' @param input summarized data for a single protein
#' @param contrast_matrix row of a contrast_matrix
#' @param coefs coefficients of a linear model (named vector)
#' @param groups unique group labels
#' @keywords internal
.getContrast = function(input, contrast, coefs, groups) {
    coef_names = names(coefs)
    intercept = grep("Intercept", coef_names, value = TRUE)
    if (length(intercept) > 0) {
        intercept_term = rep(0, length(intercept))
        names(intercept_term) = intercept
    } else {
        intercept_term = NULL
    }
    group = grep("GROUP", coef_names, value = TRUE)
    interaction = grep(":", coef_names, value = TRUE)
    group = setdiff(group, interaction)
    if (length(group) > 0) {
        group_term = contrast[, as.character(groups[groups %in% unique(input$GROUP)])]
        names(group_term) = paste0("GROUP", names(group_term))
        group_term = group_term[-1]
    } else {
        group_term = NULL
    }
    contrast = c(intercept_term, group_term)
    contrast[!is.na(coefs)]
}


#' Count percentage of missing values in given conditions
#' @param contrast_matrix contrast matrix
#' @param summarized data.table summarized by the dataProcess function
#' @param result result of groupComparison
#' @param samples_info number of runs per group
#' @param has_imputed if TRUE, missing values have been imputed by dataProcess
#' @keywords internal
.countMissingPercentage = function(contrast_matrix, summarized, 
                                   result, samples_info, has_imputed) {
    TotalGroupMeasurements = NumMeasuredFeature = NumImputedFeature = NULL
    NumFea = NumRuns = totalN = NULL
        
        counts = summarized[,
                            list(totalN = unique(TotalGroupMeasurements),
                                 NumMeasuredFeature = sum(NumMeasuredFeature, 
                                                          na.rm = TRUE),
                                 NumImputedFeature = sum(NumImputedFeature, 
                                                         na.rm = TRUE)),
                            by = "GROUP"]
    
    empty_conditions = setdiff(samples_info$GROUP, unique(counts$GROUP))
    if (length(empty_conditions) !=0) {
        counts = merge(samples_info, counts, by = "GROUP", all.x = TRUE)
        counts[, NumFea := totalN / NumRuns ] # calculate number of features, to get the expected number of measurements in missing conditions
        nofea = max(ceiling(counts$NumFea), na.rm = TRUE) # it should be integer, just in case of double, use ceiling to get interger
        counts[is.na(totalN), NumMeasuredFeature := 0]
        counts[is.na(totalN), NumImputedFeature := 0]
        counts[is.na(totalN), totalN := NumRuns * nofea]
    }
    
    missing_vector = numeric(nrow(contrast_matrix))
    imputed_vector = numeric(nrow(contrast_matrix))
    for (i in seq_len(nrow(contrast_matrix))) {
        conditions = contrast_matrix[i, ] != 0
        missing_percentage = 1 - sum(counts$NumMeasuredFeature[conditions],
                                     na.rm = TRUE) / sum(counts$totalN[conditions],
                                                         na.rm = TRUE)
        if (has_imputed) {
            imputed_percentage = sum(counts$NumImputedFeature[conditions],
                                     na.rm = TRUE) / sum(counts$totalN[conditions],
                                                         na.rm = TRUE)
            imputed_vector[i] = imputed_percentage
        }
        missing_vector[i] = missing_percentage
    }
    result$MissingPercentage = missing_vector
    if (has_imputed) {
        result$ImputationPercentage = imputed_vector
    }
    result
}

