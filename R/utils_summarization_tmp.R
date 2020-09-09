.summarizeTukey = function(input, impute, cutoff_base, censored_symbol, 
                           remove50missing, original_scale = FALSE, 
                           n_threads = NULL) {
    proteins = unique(input$PROTEIN)
    n_proteins = length(proteins)
    summarized_results = vector("list", n_proteins)
    for (protein_id in seq_along(proteins)) {
        single_protein = input[PROTEIN == proteins[protein_id], ]
        summarized_results[[protein_id]] = .summarizeTukeySingleProtein(
            single_protein, impute, cutoff_base, censored_symbol, 
            remove50missing, original_scale, n_threads)
    }
    summarized_results = data.table::rbindlist(summarized_results)
    summarized_results
}

.summarizeTukeySingleProtein = function(
    input, impute, cutoff_base, censored_symbol, remove50missing, original_scale,
    n_threads = NULL) {
    # msg = paste("Getting the summarization by Tukey's median polish",
    #             "per subplot for protein", unique(input$PROTEIN))
    # # TODO: add info about i of n proteins
    # getOption("MSstatsMsg")("INFO", msg)
    input$FEATURE = factor(input$FEATURE)	
    
    if (impute & !is.null(censored_symbol)) {
        if (censored_symbol == "0") {
            input$ABUNDANCE = ifelse(input$censored, 0, input$ABUNDANCE)
        } else if (censored_symbol == "NA") {
            input$ABUNDANCE = ifelse(input$censored, NA, input$ABUNDANCE)
        }
        input$cen = ifelse(input$censored, 0, 1)
    }
    
    if (all(is.na(input$ABUNDANCE) | input$ABUNDANCE == 0)) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN),
                    "because all measurements are NAs.")
        getOption("MSstatsMsg")("INFO", msg)
        return(NULL)
    }
    
    nonmissing_filter = .getNonMissingFilter(input, impute, censored_symbol)
    count_by_feature = input[nonmissing_filter, list(n_obs = .N), by = "FEATURE"]
    count_by_feature = count_by_feature[n_obs > 0, ]
    
    input = input[FEATURE %in% count_by_feature$FEATURE, ]
    if (nrow(input) == 0) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                    "because all measurements are NAs.")
        getOption("MSstatsMsg")("INFO", msg)
        return(NULL)
    } else {
        input[, FEATURE := factor(FEATURE)]
    }
    
    single_feature = count_by_feature[n_obs == 1, as.character(unique(FEATURE))]
    input = input[!(FEATURE %in% single_feature), ]
    if (nrow(input) == 0) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                    "because features have only one measurement across MS runs.")
        getOption("MSstatsMsg")("INFO", msg)
        return(NULL)
    } else {
        input[, FEATURE := factor(FEATURE)]
    }
    
    if (all(is.na(input$ABUNDANCE) | input$ABUNDANCE == 0) | nrow(input) == 0) {
        msg = paste("After removing features which has only 1 measurement,",
                    "Can't summarize for protein", unique(input$PROTEIN), 
                    "because all measurements are NAs.")
        getOption("MSstatsMsg")("INFO", msg)
        return(NULL)
    }
    
    nonmissing_filter = .getNonMissingFilter(input, impute, censored_symbol)
    counts = input[nonmissing_filter, list(n_obs = .N), by = "RUN"]
    counts = counts[n_obs > 0]
    
    missing_runs = setdiff(unique(input$RUN), counts$RUN)
    if (length(missing_runs) > 0 & length(intersect(missing_runs, as.character(unique(input$RUN))))) { 
        # Condition is hard to read here
        input = input[RUN %in% counts$RUN, ]
        input[, RUN := factor(RUN)]
    }
    
    if (remove50missing) {
        nonmissing_filter = .getNonMissingFilter(input, TRUE, censored_symbol)
        n_features = data.table::uniqueN(input[nonmissing_filter, FEATURE])
        n_runs = data.table::uniqueN(input$RUN)
        prop_features = input[nonmissing_filter, 
                              list(prop_features = data.table::uniqueN(FEATURE) / n_features),
                              by = "RUN"] # RUN or RUN, LABEL?
        prop_features = prop_features[prop_features <= 0.5, unique(RUN)]
        
        if (length(prop_features) == n_runs) {
            msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                        "because all runs have more than 50% NAs and",
                        "are removed with the option, remove50missing=TRUE.")
            getOption("MSstatsMsg")("INFO", msg)
            return(NULL)
        }
        input = input[!(RUN %in% unique(prop_features)), ]
        input[, RUN := factor(RUN)]
    }
    
    if (any(input$cen == 0)) {
        input = .setCensoredByThreshold(input, cutoff_base, censored_symbol, remove50missing)
        if (impute) {
            survival_fit = .fitSurvival(input[LABEL == "L", ])
            input$ABUNDANCE = ifelse(input$censored & input$LABEL == "L",
                                     predict(survival_fit, newdata = input, 
                                             type = "response"), input$ABUNDANCE)
        }
    }
    
    input = input[!is.na(ABUNDANCE), ]
    is_labeled = nlevels(input$LABEL) > 1
    result = .runTukey(input, is_labeled, censored_symbol, remove50missing,
                       original_scale, log_base = 2)
    result
}

.runTukey = function(input, is_labeled, censored_symbol, remove50missing,
                     original_scale = FALSE, log_base = 2) {
    features = as.character(unique(input$FEATURE))
    
    if (nlevels(input$FEATURE) > 1) { ## for more than 1 features
        tmp_result = .fitTukey(input, original_scale, features, log_base)
    } else { 
        if (is_labeled) {
            tmp_result = .adjustLRuns(input, TRUE)
        } else {
            tmp_result = input[input$LABEL == "L", list(RUN, LogIntensities = ABUNDANCE)]
        }
    }
    tmp_result[, Protein := unique(input$PROTEIN)]
    
    nonmissing_filter = .getNonMissingFilterStats(input, censored_symbol)
    stats = input[nonmissing_filter, 
                  list(NumMeasuredFeature = data.table::uniqueN(FEATURE)), by = "RUN"]
    stats[, MissingPercentage := 1 - NumMeasuredFeature / data.table::uniqueN(features)]
    stats[, more50missing := MissingPercentage >= 0.5]
    if (!is.null(censored_symbol)) {
        if (censored_symbol == "NA") {
            nonmissing_filter = input$LABEL == "L" & !is.na(input$INTENSITY) & !is.na(input$ABUNDANCE)
        } else {
            nonmissing_filter = input$LABEL == "L" & !is.na(input$INTENSITY) & input$censored
        }
        imputed = input[!nonmissing_filter,
                        list(NumImputedFeature = data.table::uniqueN(FEATURE)),
                        by = "RUN"]
        stats = merge(stats, imputed, by = "RUN", all.x = TRUE, sort = FALSE)
        stats[, NumImputedFeature := ifelse(is.na(NumImputedFeature), 0, 
                                            NumImputedFeature)]
    }
    tmp_result = merge(tmp_result, stats, by = "RUN", all.x = TRUE, sort = FALSE)
    tmp_result
}

.fitTukey = function(input, original_scale, features, log_base) {
    features = as.character(unique(input$FEATURE))
    wide = data.table::dcast(LABEL + RUN ~ FEATURE, data = input,
                             value.var = "ABUNDANCE", keep = TRUE)
    if (original_scale) {
        wide[, features] = wide[, lapply(.SD, function(x) log_base^x), .SDcols = features]
    }
    tmp_fit = medpolish(wide[, features, with = FALSE], na.rm = TRUE, trace.iter = FALSE)
    tmp_result = data.table::data.table(LABEL = wide$LABEL,
                                        RUN = wide$RUN,
                                        ABUNDANCE = tmp_fit$overall + tmp_fit$row)
    if (original_scale) {
        tmpt_result[, ABUNDANCE := log(ABUNDANCE, base_log)]
    }
    
    if (data.table::uniqueN(input$LABEL) == 2) {
        tmp_result = .adjustLRuns(tmp_result)
    }
    tmp_result[, list(RUN, LogIntensities = ABUNDANCE)]
}

.adjustLRuns = function(input, rename = FALSE) {
    h_runs = input[LABEL == "H", list(RUN, ABUNDANCE)]
    h_median = median(input[LABEL == "H", ABUNDANCE], na.rm = TRUE)
    input = input[LABEL == "L"]
    input = merge(input[, list(RUN, ABUNDANCE)], h_runs, by = "RUN", suffixes = c("", ".h"))
    input[, ABUNDANCE := ABUNDANCE - ABUNDANCE.h + h_median]
    if (rename) {
        input[, list(RUN, LogIntensities = ABUNDANCE)]
    } else {
        input[, list(RUN, ABUNDANCE)]
    }
}

.getNonMissingFilterStats = function(input, censored_symbol) {
    if (!is.null(censored_symbol)) {
        if (censored_symbol == "NA") {
            nonmissing_filter = input$LABEL == "L" & !is.na(input$INTENSITY)
        } else {
            nonmissing_filter = input$LABEL == "L" & !is.na(input$INTENSITY)
            nonmissing_filter = nonmissing_filter & input$INTENSITY > 1
            nonmissing_filter = nonmissing_filter & !is.na(input$ABUNDANCE) & input$ABUNDANCE > 0
        }
    } else {
        nonmissing_filter = input$LABEL == "L" & !is.na(input$INTENSITY)
    }
    nonmissing_filter
}
