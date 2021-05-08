
MSstatsPrepareForSummarization = function(input, method, impute, censored_symbol,
                                          remove_uninformative_feature_outlier) {
    ABUNDANCE = feature_quality = is_outlier = NULL
    
    input = .removeSingleLabelFeatures(input)
    label = data.table::uniqueN(input$LABEL) == 2
    if (label) {
        input$ref = factor(ifelse(input$LABEL == "L", 
                                  input$RUN[input$LABEL == "L"], 0))
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


getProcessed = function(input) {
    if (is.element("remove", colnames(input))) {
        input(remove)
    } else {
        NULL
    }
}


.prepareSummary = function(input, method, impute, censored_symbol) {
    if (method == "TMP") {
        input = .prepareTMP(input, impute, censored_symbol)
    } else {
        input = .prepareLinear(input, FALSE, censored_symbol)
    }
    input
}


.prepareLinear = function(input, impute, censored_symbol) {
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
.prepareTMP = function(input, impute, censored_symbol) {
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
    # remove feature with 1 measurement
    input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)] 
    input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]
    
    input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
    input[, prop_features := sum(nonmissing) / total_features,
          by = c("PROTEIN", "RUN")] 
    
    if (any(input$cen == 0)) {
        .setCensoredByThreshold(input, censored_symbol, remove50missing)
    }
    input
}
