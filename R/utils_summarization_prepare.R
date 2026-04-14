#' Prepare feature-level data for protein-level summarization
#' 
#' @param input feature-level data processed by dataProcess subfunctions
#' @param method summarization method - `summaryMethod` parameter of the dataProcess function
#' @param impute if TRUE, censored missing values will be imputed - `MBimpute`
#' parameter of the dataProcess function
#' @param censored_symbol censored missing value indicator - `censoredInt` 
#' parameter of the dataProcess function
#' @param remove_uninformative_feature_outlier if TRUE, features labeled as 
#' outlier of uninformative by the MSstatsSelectFeatures function will not be 
#' used in summarization
#' 
#' @return data.table
#' 
#' @export
#' 
#' @examples 
#' raw = DDARawData 
#' method = "TMP"
#' cens = "NA"
#' impute = TRUE
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input = MSstatsPrepareForDataProcess(raw, 2, NULL)
#' head(input)
#' 
MSstatsPrepareForSummarization = function(input, method, impute, censored_symbol,
                                          remove_uninformative_feature_outlier) {
    ABUNDANCE = feature_quality = is_outlier = PROTEIN = NULL
    
    if (!"ANOMALYSCORES" %in% colnames(input)) {
        input[, ANOMALYSCORES := NA]
    }
    
    add_ref_covariate = "is_labeled_ref" %in% colnames(input) && any(input$is_labeled_ref, na.rm = TRUE)
    if (add_ref_covariate) {
      input[, ref_covariate := factor(ifelse(LABEL == "L", RUN, 0))]
    }
    
    if (is.element("remove", colnames(input))) {
        input = input[!(remove)]
    }
    
    if (remove_uninformative_feature_outlier & 
        is.element("feature_quality", colnames(input))) {
        input[, ABUNDANCE := ifelse(feature_quality == "Uninformative", 
                                    NA, ABUNDANCE)]
        input[, ABUNDANCE := ifelse(is_outlier, NA, ABUNDANCE)]
        msg = "** Filtered out uninformative features and outliers."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    input = .prepareSummary(input, impute, censored_symbol, add_ref_covariate)
    input[, PROTEIN := factor(PROTEIN)]
    input
}


#' Get feature-level data to be used in the MSstatsSummarizationOutput function
#' 
#' @param input data.table processed by dataProcess subfunctions
#' 
#' @return data.table processed by dataProcess subfunctions
#' 
#' @export
#' 
#' @examples 
#' raw = DDARawData 
#' method = "TMP"
#' cens = "NA"
#' impute = TRUE
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input = MSstatsPrepareForDataProcess(raw, 2, NULL)
#' input = MSstatsNormalize(input, "EQUALIZEMEDIANS")
#' input = MSstatsMergeFractions(input)
#' input = MSstatsHandleMissing(input, "TMP", TRUE, "NA", 0.999)
#' input_all = MSstatsSelectFeatures(input, "all") # all features
#' input_5 = MSstatsSelectFeatures(data.table::copy(input), 
#' "topN", top_n = 5) # top 5 features
#' 
#' proc1 = getProcessed(input_all)
#' proc2 = getProcessed(input_5)
#' 
#' proc1
#' proc2
#' 
getProcessed = function(input) {
    remove = NULL
    
    if (is.element("remove", colnames(input))) {
        if (all(!(input$remove))) {
            NULL
        } else {
            input[(remove)]
        }
    } else {
        NULL
    }
}


#' Prepare feature-level data for summarization
#' @param input data.table
#' @param impute logical
#' @param censored_symbol "0"/"NA"
#' @param is_labeled_reference logical, if TRUE the H channel is a normalization
#'   reference (SRM) and grouping keys do not include LABEL; if FALSE (e.g.
#'   protein turnover) LABEL is added to grouping keys so each label is
#'   processed independently.
#' @return data.table
#' @keywords internal
.prepareSummary = function(input, impute, censored_symbol,
                           is_labeled_reference = FALSE) {
    censored = feature_quality = newABUNDANCE = cen = nonmissing = n_obs = NULL
    n_obs_run = total_features = FEATURE = prop_features = NULL
    remove50missing = ABUNDANCE = NULL
    
    label_by = if (is_labeled_reference) character(0) else "LABEL"
    
    if (impute & !is.null(censored_symbol)) {
        if (is.element("feature_quality", colnames(input))) {
            input[, censored := ifelse(feature_quality == "Informative",
                                       censored, FALSE)]
        }
        if (censored_symbol == "0") {
            input[, newABUNDANCE := ifelse(censored, 0, ABUNDANCE)]
        } else if (censored_symbol == "NA") {
            input[, newABUNDANCE := ifelse(censored, NA, ABUNDANCE)]
        }
        input[, cen := ifelse(censored, 0, 1)]
    } else {
        input[, newABUNDANCE := ABUNDANCE]
    }
    
    input[, nonmissing := .getNonMissingFilter(input, impute, censored_symbol)]
    input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE", label_by)]
    input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)]
    input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN", label_by)]
    
    input[, total_features := uniqueN(FEATURE), by = c("PROTEIN", label_by)]
    input[, prop_features := sum(nonmissing) / total_features,
          by = c("PROTEIN", "RUN", label_by)]
    
    if (is.element("cen", colnames(input))) {
        if (any(input[["cen"]] == 0)) {
            .setCensoredByThreshold(input, censored_symbol, remove50missing)
        }
    }
    
    input
}
