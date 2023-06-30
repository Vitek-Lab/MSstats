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
    
    label = data.table::uniqueN(input$LABEL) == 2
    if (label) {
        input[, ref := factor(ifelse(LABEL == "L", RUN, 0))]
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
    
    input = .prepareSummary(input, method, impute, censored_symbol)
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
#' @param method "TMP" / "linear"
#' @param impute logical
#' @param censored_symbol "0"/"NA"
#' @return data.table
#' @keywords internal
.prepareSummary = function(input, method, impute, censored_symbol) {
    if (method == "TMP") {
        input = .prepareTMP(input, impute, censored_symbol)
    } else {
        input = .prepareLinear(input, FALSE, censored_symbol)
    }
    input
}


#' Prepare feature-level data for linear summarization
#' @inheritParams .prepareSummary
#' @return data.table
#' @keywords internal
.prepareLinear = function(input, impute, censored_symbol) {
    newABUNDANCE = ABUNDANCE = nonmissing = n_obs = n_obs_run = NULL
    total_features = FEATURE = prop_features = NULL
    
    input[, newABUNDANCE := ABUNDANCE]
    input[, nonmissing := .getNonMissingFilter(.SD, impute, censored_symbol)]
    input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE")]
    # remove feature with 1 measurement
    input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)] 
    input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]
    
    input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
    input[, prop_features := sum(nonmissing) / total_features,
          by = c("PROTEIN", "RUN")] 
    input
}


#' Prepare feature-level data for TMP summarization
#' @inheritParams .prepareSummary
#' @return data.table
#' @keywords internal
.prepareTMP = function(input, impute, censored_symbol) {
    censored = feature_quality = newABUNDANCE = cen = nonmissing = n_obs = NULL
    n_obs_run = total_features = FEATURE = prop_features = NULL
    remove50missing = ABUNDANCE = NULL
    
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
    input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE")]
    input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)] 
    input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]
    
    input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
    input[, prop_features := sum(nonmissing) / total_features,
          by = c("PROTEIN", "RUN")] 
    
    if (is.element("cen", colnames(input))) {
        if (any(input[["cen"]] == 0)) {
            .setCensoredByThreshold(input, censored_symbol, remove50missing)
        }
    }

    input
}
