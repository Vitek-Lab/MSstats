#' Check if a protein can be summarized with TMP
#' @param input data.table
#' @param remove50missing if TRUE, proteins with more than 50% missing values
#' in all runs will not be summarized
#' @return data.table
#' @keywords internal 
.isSummarizable = function(input, remove50missing) {
    n_obs_run = RUN = NULL
    
    if (all(is.na(input$newABUNDANCE) | input$newABUNDANCE == 0)) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN),
                    "because all measurements are missing or censored.")
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg)
        return(NULL)
    }
    
    if (all(is.na(input$n_obs) | input$n_obs == 0)) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                    "because all measurements are missing or censored.")
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg)
        return(NULL)
    } 
    
    if (all(input$n_obs == 1 | is.na(input$n_obs))) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                    "because features have only one measurement across MS runs.")
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg)
        return(NULL)
    }
    
    if (all(is.na(input$newABUNDANCE) | input$newABUNDANCE == 0) | nrow(input) == 0) {
        msg = paste("After removing features which has only 1 measurement,",
                    "Can't summarize for protein", unique(input$PROTEIN), 
                    "because all measurements are missing or censored.")
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg)
        return(NULL)
    }
    
    missing_runs = setdiff(unique(input$RUN), 
                           unique(input[n_obs_run == 0 | is.na(n_obs_run), RUN]))
    if (length(missing_runs) > 0 & length(intersect(missing_runs, as.character(unique(input$RUN))))) { 
        input = input[n_obs_run > 0 & !is.na(n_obs_run), ]
    }
    
    if (remove50missing) {
        if (all(input$prop_features <= 0.5 | is.na(input$prop_features))) {
            msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                        "because all runs have more than 50% missing values and",
                        "are removed with the option, remove50missing=TRUE.")
            getOption("MSstatsMsg")("INFO", msg)
            getOption("MSstatsLog")("INFO", msg)
            return(NULL)
        }
    }
    input
}


#' Fit Tukey median polish
#' @param input data.table with data for a single protein
#' @param is_labeled_reference logical, if TRUE, H channel is used as a
#'   normalization reference (SRM experiment): L abundances are adjusted by
#'   subtracting the H value and adding back the H median, and only L results
#'   are returned. If FALSE (e.g. protein turnover), each label is summarized
#'   independently and results for all labels are returned.
#' @inheritParams MSstatsSummarizeWithSingleCore
#' @return data.table
#' @keywords internal
.runTukey = function(input, is_labeled_reference, censored_symbol, remove50missing) {
    Protein = RUN = newABUNDANCE = NULL

    if (nlevels(input$FEATURE) > 1) {
        tmp_result = .fitTukey(input, is_labeled_reference)
    } else {
        if (is_labeled_reference) {
            tmp_result = .adjustLRuns(input, TRUE)
            tmp_result[, LABEL := "L"]
        } else {
            tmp_result = input[, list(LABEL, RUN, LogIntensities = newABUNDANCE)]
        }
    }
    tmp_result[, Protein := unique(input$PROTEIN)]
    tmp_result
}


#' Fit tukey median polish for a data matrix
#' @inheritParams .runTukey
#' @return data.table
#' @keywords internal
.fitTukey = function(input, is_labeled_reference) {
    LABEL = RUN = newABUNDANCE = NULL

    features = as.character(unique(input$FEATURE))
    wide = data.table::dcast(LABEL + RUN ~ FEATURE, data = input,
                             value.var = "newABUNDANCE", keep = TRUE)
    tmp_fitted = median_polish_summary(as.matrix(wide[, features, with = FALSE]))
    wide[, newABUNDANCE := tmp_fitted]
    tmp_result = wide[, list(LABEL, RUN, newABUNDANCE)]

    if (is_labeled_reference) {
        tmp_result = .adjustLRuns(tmp_result)
        tmp_result[, list(LABEL = "L", RUN, LogIntensities = newABUNDANCE)]
    } else {
        tmp_result[, list(LABEL, RUN, LogIntensities = newABUNDANCE)]
    }
}


#' Adjust summarized abundance based on the heavy channel
#' @param input data.table
#' @param rename if TRUE, rename the output column to LogIntensities
#' @return data.table
#' @importFrom stats median
#' @keywords internal
.adjustLRuns = function(input, rename = FALSE) {
    LABEL = newABUNDANCE = RUN = newABUNDANCE.h = NULL
    
    h_runs = input[LABEL == "H", list(RUN, newABUNDANCE)]
    h_median = median(input[LABEL == "H", newABUNDANCE], na.rm = TRUE)
    input = input[LABEL == "L"]
    input = merge(input[, list(RUN, newABUNDANCE)], h_runs, by = "RUN", suffixes = c("", ".h"))
    input[, newABUNDANCE := newABUNDANCE - newABUNDANCE.h + h_median]
    if (rename) {
        input[, list(RUN, LogIntensities = newABUNDANCE)]
    } else {
        input[, list(RUN, newABUNDANCE)]
    }
}


#' Get a logical vector for non-missing values to calculate summary statistics
#' @inheritParams .runTukey
#' @return data.table
#' @keywords internal
.getNonMissingFilterStats = function(input, censored_symbol) {
    if (!is.null(censored_symbol)) {
        nonmissing_filter = !is.na(input$newABUNDANCE) & !input$censored
    } else {
        nonmissing_filter = !is.na(input$INTENSITY)
    }
    nonmissing_filter = nonmissing_filter & input$n_obs_run > 0 & input$n_obs > 1
    nonmissing_filter
}


#' Fit a linear model
#' @param input data.table
#' @param is_single_feature logical, if TRUE, data has single feature
#' @param is_labeled logical, if TRUE, data comes from a labeled experiment
#' @param equal_variances logical, if TRUE, equal variances are assumed
#' @return lm or merMod
#' @importFrom stats lm
#' @keywords internal
.fitLinearModel = function(input, is_single_feature, is_labeled,
                           equal_variances) {
    if (all(!is.na(input$weights))){
        weight_input = input$weights
    } else {
        weight_input = NULL
    }
    
    if (!is_labeled) {
        if (is_single_feature) {
            linear_model = lm(newABUNDANCE ~ RUN , data = input, 
                              weights=weight_input)
        } else {
            linear_model = lm(newABUNDANCE ~ FEATURE + RUN, 
                              data = input, weights=weight_input)
        }
    } else {
        if (is_single_feature) {
            linear_model = lm(ABUNDANCE ~ RUN + ref_covariate, data = input)
        } else {
            linear_model = lm(ABUNDANCE ~ FEATURE + RUN + ref_covariate, data = input)
        }
    }
    if (!equal_variances) {
        linear_model = .updateUnequalVariances(input = input,
                                               fit = linear_model,
                                               num_iter = 1)
    }
    linear_model$model = NULL
    linear_model
}


#' Adjust model for unequal variances
#' @param input data.table
#' @param fit lm
#' @param num_iter number of iterations
#' @importFrom lme4 lmer
#' @importFrom stats loess resid lm formula
#' @return merMod
#' @keywords internal
.updateUnequalVariances = function(input, fit, num_iter) {
    weight = NULL
    input = as.data.frame(input)
    for (i in seq_len(num_iter)) {
        if (i == 1) {
            input[["abs.resids"]] = abs(fit$residuals)
            input[["fitted"]] = fit$fitted.values
        }
        fit.loess = loess(abs.resids ~ fitted, data = input)
        input[["loess.fitted"]] = fitted(fit.loess)
        ## loess fitted values are predicted sd
        input[["weight"]] = 1 / (input[["loess.fitted"]] ^ 2)
        input[["abs.resids"]] = NULL
        ## re-fit using weight
        wls.fit = lm(formula(fit), data = input, weights = weight)
        input[["abs.resids"]] = abs(wls.fit$residuals)
        input[["loess.fitted"]] = NULL
        input[["weight"]] = NULL
    }
    wls.fit
}


#' Check if data has less than two features
#' @param input data.table
#' @return logical
#' @keywords internal
.checkSingleFeature = function(input) {
    data.table::uniqueN(input$FEATURE) < 2
}
