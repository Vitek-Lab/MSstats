#' Feature selection before feature-level data summarization
#' 
#' @param input data.table
#' @param method "all" / "highQuality", "topN"
#' @param top_n number of features to use for "topN" method
#' @param minimum number of quality features for "highQuality" method
#' 
#' @return data.table
#' 
#' @export
#' 
MSstatsSelectFeatures = function(input, method, top_n = NULL, min_feature_count = NULL) {
    if (method == "all") {
        msg = "** Use all features that the dataset originally has."
    } else if (method == "highQuality") {
        msg = "** Flag uninformative feature and outliers by feature selection algorithm."
        features_quality = .selectHighQualityFeatures(input, min_feature_count)
        input = merge(input, features_quality, all.x = TRUE,
                      by.x = c("LABEL", "PROTEIN", "FEATURE", "originalRUN"),
                      by.y = c("label", "protein", "feature", "run"))
        input$feature_quality = ifelse(is.na(input$feature_quality),
                                                  "Uninformative", input$feature_quality) # is this OK?
        input$is_outlier = ifelse(is.na(input$is_outlier),
                                             TRUE, input$is_outlier)
    } else if (method %in% c("top3", "topN")) {
        msg = paste0("** Use top", top_n, " features that have highest average of log2(intensity) across runs.")
        input = .selectTopFeatures(input, top_n)
    }
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input
}

.selectTopFeatures = function(input, top_n) {
    ABUNDANCE = MeanAbundance = remove = FEATURE = NULL
    
    mean_by_feature = input[ABUNDANCE > 0, 
                            list(MeanAbundance = mean(ABUNDANCE, na.rm = TRUE)),
                            by = c("PROTEIN", "FEATURE")]
    mean_by_feature[, feature_rank := rank(-MeanAbundance), by = "PROTEIN"]
    mean_by_feature = mean_by_feature[feature_rank <= top_n, ]
    input[, remove := !(FEATURE %in% mean_by_feature$FEATURE)]
    input
}


#' @keywords internal
.selectHighQualityFeatures = function(input, min_feature_count) {
    PROTEIN = PEPTIDE = FEATURE = originalRUN = ABUNDANCE = is_censored = NULL
    
    cols = c("PROTEIN", "PEPTIDE", "FEATURE", "originalRUN", "LABEL", 
             "ABUNDANCE", "censored")
    cols = intersect(cols, colnames(input))
    input = input[, cols, with = FALSE]
    if (!("censored" %in% cols)) {
        input$censored = FALSE
    } 
    data.table::setnames(input, "censored", "is_censored")
    input = input[, list(protein = as.character(PROTEIN),
                         peptide = as.character(PEPTIDE),
                         feature = as.character(FEATURE),
                         run = as.character(originalRUN),
                         label = as.character(LABEL),
                         log2inty = ifelse(!(is.na(ABUNDANCE) | is_censored),
                                           ABUNDANCE, NA),
                         is_censored)]
    input[, is_obs := !(is.na(log2inty) | is_censored)]
    input[, is_censored := NULL]
    
    features_quality = data.table::rbindlist(lapply(split(input, input$label),
                                                    .flagUninformativeSingleLabel,
                                                    min_feature_count = min_feature_count))
    features_quality
}


.flagUninformativeSingleLabel = function(input, min_feature_count) {
    
    min_feature_count = 3
    
    if (nrow(input) == 0) {
        return(NULL)
    }
    if (unique(input$label) == "H") {
        input = input[log2inty > 0, ] # Alternatively, filter them out in the next steps
    }
    
    input[, n_observed := sum(is_obs), by = c("protein", "feature")]
    input[, any_observed := any(is_obs), by = c("protein", "run")]
    input[, singleton := n_observed <= 1]
    input[, few_observed := n_observed <= min_feature_count, ]
    input[, unrep := !any_observed, ]
    
    .addOutlierCutoff(input)
    .addCoverageInfo(input)
    .addNInformativeInfo(input, min_feature_count)
    .addModelInformation(input)
    input = .addModelVariances(input)
    input = .addNoisyFlag(input)
    
    input$feature_quality = ifelse(
        !input$unrep & !input$is_lowcvr & !input$is_noisy & !input$few_observed, 
        "Informative", "Uninformative"
    )
    input$is_outlier = ifelse(input$label == "H" & input$log2inty <= 0,
                              TRUE, input$is_outlier)
    input = unique(input[, list(label, protein, feature, run, feature_quality, is_outlier)])
    input
}

.addOutlierCutoff = function(input, quantile_order = 0.01) {
    input[, 
          min_obs := .calculateOutlierCutoff(.SD, quantile_order), 
          by = "protein",
          .SDcols = c("run", "feature", "is_obs",
                      "unrep", "singleton", "few_observed")]
}

.calculateOutlierCutoff = function(input, quantile_order = 0.01) {
    n_runs = data.table::uniqueN(input[!unrep & !few_observed & !singleton, run])
    n_features = data.table::uniqueN(input[!unrep & !few_observed & !singleton, feature])
    n = input[!unrep & !few_observed & !singleton, sum(is_obs, na.rm = TRUE)]
    qbinom(quantile_order, n_runs,
           n / (n_features * n_runs))
}


.addCoverageInfo = function(input) {
    input[, is_lowcvr := .flagLowCoverage(.SD), by = c("protein", "feature"),
          .SDcols = c("is_obs", "min_obs")]
}

.flagLowCoverage = function(input) {
    sum(input$is_obs) < unique(input$min_obs)
}

.addNInformativeInfo = function(input, min_feature_count) {
    input[, n_informative := .countInformative(.SD), by = "protein",
          .SDcols = c("feature", "is_lowcvr")]
    input[, has_three_informative := n_informative >= min_feature_count]
    # input[, n_informative := NULL]
}

#' @importFrom data.table uniqueN
.countInformative = function(input) {
    data.table::uniqueN(input[!(is_lowcvr), feature])
}


.addModelInformation = function(input) {
    input[, c("model_residuals", "df_resid", "var_resid") := .calculateProteinVariance(.SD),
          by = "protein", .SDcols = c("protein", "log2inty", "run", "feature", "is_lowcvr")]
    
}

.calculateProteinVariance = function(input) {
    robust_model = tryCatch(.fitHuber(input[(!is_lowcvr), ]),
                            error = function(e) NULL,
                            warning = function(w) NULL) 
    if (is.null(robust_model)) {
        list(NA_real_, NA_real_, NA_real_)
    } else {
        model_residuals = rep(NA_real_, nrow(input))
        model_residuals[!input$is_lowcvr & !is.na(input$log2inty)] = residuals(robust_model)
        list(as.numeric(model_residuals),
             rep(summary(robust_model)$df[2], nrow(input)),
             rep(summary(robust_model)$sigma ^ 2, nrow(input)))
    }
}


#' Wrapper to fit robust linear model for one protein
#' @importFrom MASS rlm
#' @keywords internal
.fitHuber = function(input) {
    MASS::rlm(log2inty ~ run + feature, data = input, scale.est = "Huber")
}

#' @importFrom limma squeezeVar
.addModelVariances = function(input) {
    model_variances = unique(input[, list(protein, df_resid, var_resid)])
    model_variances = model_variances[!is.na(df_resid), ]
    
    if (nrow(model_variances) > 0) {
        eb_fit = limma::squeezeVar(model_variances$var_resid, model_variances$df_resid, 
                                   robust = TRUE)
        model_variances$var_resid_eb = eb_fit$var.post
        model_variances$s_resid_eb = sqrt(eb_fit$var.post)
        model_variances = model_variances[, list(protein, s_resid_eb)]
        input = merge(input, model_variances, by = "protein", all.x = TRUE)
    } else {
        input$s_resid_eb = NA_real_
    }
    input
}

.addNoisyFlag = function(input) {
    feature_vars = .getFeatureVariances(input)
    if (nrow(feature_vars) > 0) {
        input = merge(input, feature_vars, by = "feature")
        input[, is_outlier := .addOutlierInformation(.SD), by = "feature"]
        input[, is_noisy := svar_feature > quantile(svar_ref, 0.05, na.rm = TRUE)]
        input
    } else {
        input[, is_outlier := NA]
        input[, is_noisy := NA]
        input
    }
}


.getFeatureVariances = function(input, tolerance = 3) {
    remove_outliers = unique(input$label) == "L"
    if (remove_outliers) {
        input[, num_filter := abs(model_residuals / s_resid_eb) > tolerance]
    } else {
        input[, num_filter := rep(TRUE, .N)]
    }
    input = input[!(num_filter) & !is.na(model_residuals)]
    input[, resid_null := log2inty - mean(log2inty, na.rm = TRUE)]
    input[, n_runs := data.table::uniqueN(run), by = "feature"]
    
    sums_of_squares = input[, list(
        svar_feature = sum(model_residuals ^ 2, na.rm = TRUE) / (n_runs - 1) / unique(s_resid_eb ^ 2),
        svar_ref = sum(resid_null ^ 2, na.rm = TRUE) / (n_runs - 1) / unique(s_resid_eb) ^ 2), by = "feature"]
    unique(sums_of_squares)
}


.addOutlierInformation = function(input, tol = 3, keep_run = FALSE) {
    input$result = abs(input$model_residuals / input$s_resid_eb) > tol
    if (keep_run) {
        input[, all_missing := all(result | is.na(result)), by = "run"]
        input[, result := ifelse(all_missing, FALSE, result)]
    }
    if (is.na(unique(input$s_resid_eb))) {
        input$result = rep(FALSE, nrow(input))
    }
    input$result
}

