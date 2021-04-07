#' Prepare data for a single protein for group comparison
#' @param single_protein data.table
#' @keywords internal
.prepareSingleProteinForGC = function(single_protein) {
    ABUNDANCE = GROUP = SUBJECT = SUBJECT_NESTED = RUN = NULL
    
    data.table::setnames(single_protein,
                         c("LogIntensities"),
                         c("ABUNDANCE"),
                         skip_absent = TRUE)
    single_protein = single_protein[!is.na(ABUNDANCE)]
    single_protein[, GROUP := factor(GROUP)]
    single_protein[, SUBJECT := factor(SUBJECT)]
    single_protein[, SUBJECT_NESTED := factor(SUBJECT_NESTED)]
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
#' @param save_fitted_models if TRUE, fitted model will be saved. If FALSE,
#' it will be replaced by NULL
#' @param processed data.table of processed data
#' @param has_imputed if TRUE, missing values have been imputed by dataProcess
#' @keywords internal
.fitModelSingleProtein = function(input, contrast_matrix, has_tech_replicates,
                                  is_single_subject, repeated, groups,
                                  save_fitted_models, processed, has_imputed) {
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
        result = .countMissingPercentage(contrast_matrix, processed,
                                         input, result, has_imputed)
        
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
                lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data = input),
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
                    lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data = input),
                    TRUE
                ))
                df_full = suppressMessages(try(
                    lm(ABUNDANCE ~ GROUP + SUBJECT, data = input)$df.residual,
                    TRUE))
            } else {
                full_fit = suppressMessages(try(
                    lmer(ABUNDANCE ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT),
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
#' @keywords internal
.getModelParameters = function(fitted_model) {
    if (class(fitted_model[["full_fit"]]) == "lm") {
        model_summary = summary(fitted_model[["full_fit"]])
        cf = model_summary[["coefficients"]]
        vcv = model_summary[["cov.unscaled"]] * (model_summary[["sigma"]] ^ 2)
    } else {
        cf = as.matrix(fixef(fitted_model[["full_fit"]]))
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
    all_comparisons = lapply(1:nrow(contrast_matrix), function(row_id) {
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
    for (row_id in 1:nrow(contrast_matrix)) {
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
#' @param ith_contrast single row of a contrast matrix
#' @param groups unique labels of experimental conditions
#' @param parameters parameters extracted from the model 
#' @param protein name of a protein
#' @param empty_conditions labels of empty conditions
#' @param coefs coefficient of the fitted model
#' @keywords internal
.handleEmptyConditions = function(input, fit, ith_contrast,
                                  groups, parameters, protein,
                                  empty_conditions, coefs) {
    count_diff_pos = intersect(groups[ith_contrast != 0 & ith_contrast > 0],
                               empty_conditions)
    count_diff_neg = intersect(groups[ith_contrast != 0 & ith_contrast < 0],
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
                      Label = row.names(ith_contrast),
                      SE = NA, Tvalue = NA, DF = NA, pvalue = NA,
                      issue = issue)
    } else {
        result = .handleSingleContrast(input, fit, ith_contrast, groups,
                                       parameters, protein, coefs)
    }
    result
}


#' Group comparison for a single contrast
#' @inheritParams .handleEmptyConditions
#' @keywords internal
.handleSingleContrast = function(input, fit, contrast_matrix, groups,
                                 parameters, protein, coefs) {
    contrast = .getContrast(input, contrast_matrix, coefs)
    result = get_estimable_fixed_random(parameters, contrast)
    if (is.null(result)) {
        result = list(Protein = protein,
                      Label = row.names(contrast_matrix),
                      logFC = NA, SE = NA, Tvalue = NA,
                      DF = NA, pvalue = NA, issue = NA)
    } else {
        result$Protein = protein
        result$Label = row.names(contrast_matrix)
        result$issue = NA
    }
    result
}


#' Create a contrast for a model with only group as a fixed effect
#' @param input summarized data for a single protein
#' @param contrast_matrix row of a contrast_matrix
#' @param coefs coefficients of a linear model (named vector)
#' @keywords internal
.getContrast = function(input, contrast_matrix, coefs) {
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
        # THE LINE BELOW WILL REQUIRE CHANGE WHEN SWITCHING TO V4
        group_term = contrast_matrix[, as.numeric(levels(input$GROUP))]
        group_term = group_term[-1]
        names(group_term) = group
    } else {
        group_term = NULL
    }
    contrast = c(intercept_term, group_term)
    contrast[!is.na(coefs)]
}


#' Count percentage of missing values in given conditions
#' @param contrast_matrix contrast matrix
#' @param processed data.table processed by the dataProcess function
#' @param summarized data.table summarized by the dataProcess function
#' @param result result of groupComparison
#' @param has_imputed if TRUE, missing values have been imputed by dataProcess
#' @keywords internal
.countMissingPercentage = function(contrast_matrix, processed, 
                                   summarized, result, has_imputed) {
    counts = processed[LABEL == "L", 
                       .(totalN = .N,
                         NumMeasuredFeature = sum(!is.na(ABUNDANCE) & ABUNDANCE > 0,
                                                  na.rm = TRUE),
                         NumImputedFeature = sum(censored & (!is.na(ABUNDANCE) & ABUNDANCE > 0),
                                                 na.rm=  TRUE)), 
                       by = "GROUP_ORIGINAL"]
    missing_vector = numeric(nrow(contrast_matrix))
    imputed_vector = numeric(nrow(contrast_matrix))
    for (i in 1:nrow(contrast_matrix)) {
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
