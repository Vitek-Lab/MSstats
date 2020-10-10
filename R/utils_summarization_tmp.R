#' Summarization with Tukey median polish
#' @importFrom utils setTxtProgressBar
.summarizeTukey = function(input, impute, cutoff_base, censored_symbol, 
                           remove50missing, original_scale = FALSE, 
                           n_threads = NULL) {
    if (impute & !is.null(censored_symbol)) {
        if (censored_symbol == "0") {
            input[, newABUNDANCE := ifelse(censored, 0, ABUNDANCE)]
        } else if (censored_symbol == "NA") {
            input[, newABUNDANCE := ifelse(censored, NA, ABUNDANCE)]
        }
        input[, cen := ifelse(censored, 0, 1)]
    }
    
    input[, nonmissing := .getNonMissingFilter(input, impute, censored_symbol)]
    input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE")]
    # remove feature with 1 measurement
    input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)] 
    input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]
    
    input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
    input[, 
          prop_features := sum(nonmissing) / total_features,
          by = c("PROTEIN", "RUN")] 
    
    if (any(input$cen == 0)) {
        .setCensoredByThreshold(input, cutoff_base, censored_symbol, remove50missing)
    }
    
    if (impute) {
        input[n_obs_run > 0 & n_obs > 1, predicted := .addSurvivalPredictions(.SD),
              by = "PROTEIN"]
        input[, predicted := ifelse(censored & (LABEL == "L"), predicted, NA)]
        input[, newABUNDANCE := ifelse(censored & LABEL == "L",
                                       predicted, newABUNDANCE)]
    }
    
    input[, NonMissingStats := .getNonMissingFilterStats(.SD, censored_symbol)]
    input[, NumMeasuredFeature := sum(NonMissingStats), 
          by = c("PROTEIN", "RUN")]
    input[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
    input[, more50missing := MissingPercentage >= 0.5]
    if (!is.null(censored_symbol)) {
        if (censored_symbol == "NA") {
            #input[, nonmissing_orig := LABEL == "L" & !is.na(INTENSITY) & !is.na(ABUNDANCE)]
            input[, nonmissing_orig := LABEL == "L" & !censored]
        } else {
            if (is.element("censored", colnames(input))) {
                input[, nonmissing_orig := LABEL == "L" & !censored]
            } else {
                input[, nonmissing_orig := LABEL == "L" & !is.na(INTENSITY)]
            }
        }
        if (impute) {
            input[, NumImputedFeature := sum(!nonmissing_orig),
                  by = c("PROTEIN", "RUN")]
        } else {
            input[, NumImputedFeature := 0]
        }
    }
    
    proteins = unique(input$PROTEIN)
    n_proteins = length(proteins)
    summarized_results = vector("list", n_proteins)
    
    pb = utils::txtProgressBar(min = 1, max = n_proteins, style = 3)
    for (protein_id in seq_along(proteins)) {
        single_protein = input[PROTEIN == proteins[protein_id]]
        single_protein = single_protein[(n_obs > 1 & !is.na(n_obs)) &
                                            (n_obs_run > 0 & !is.na(n_obs_run))]
        single_protein[, RUN := factor(RUN)]
        single_protein[, FEATURE := factor(FEATURE)]
        
        summarized_results[[protein_id]] = .summarizeTukeySingleProtein(
            single_protein[, list(PROTEIN, LABEL, RUN, FEATURE, newABUNDANCE,
                                  n_obs, n_obs_run, prop_features)],
            impute, cutoff_base, censored_symbol,
            remove50missing, original_scale, n_threads)
        setTxtProgressBar(pb, protein_id)
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
    
    if (all(is.na(input$newABUNDANCE) | input$newABUNDANCE == 0)) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN),
                    "because all measurements are missing or censored.")
        getOption("MSstatsMsg")("INFO", msg)
        return(NULL)
    }
    
    if (all(is.na(input$n_obs) | input$n_obs == 0)) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                    "because all measurements are missing or censored.")
        getOption("MSstatsMsg")("INFO", msg)
        return(NULL)
    } 
    
    if (all(input$n_obs == 1 | is.na(input$n_obs))) {
        msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                    "because features have only one measurement across MS runs.")
        getOption("MSstatsMsg")("INFO", msg)
        return(NULL)
    }
    
    if (all(is.na(input$newABUNDANCE) | input$newABUNDANCE == 0) | nrow(input) == 0) {
        msg = paste("After removing features which has only 1 measurement,",
                    "Can't summarize for protein", unique(input$PROTEIN), 
                    "because all measurements are missing or censored.")
        getOption("MSstatsMsg")("INFO", msg)
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
            return(NULL)
        }
    }
    
    input[, RUN := factor(RUN)]
    input[, FEATURE := factor(FEATURE)]
    
    input = input[!is.na(newABUNDANCE), ]
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
            tmp_result = input[input$LABEL == "L", list(RUN, LogIntensities = newABUNDANCE)]
        }
    }
    tmp_result[, Protein := unique(input$PROTEIN)]
    tmp_result
}

.fitTukey = function(input, original_scale, features, log_base) {
    features = as.character(unique(input$FEATURE))
    wide = data.table::dcast(LABEL + RUN ~ FEATURE, data = input,
                             value.var = "newABUNDANCE", keep = TRUE)
    if (original_scale) {
        wide[, features] = wide[, lapply(.SD, function(x) log_base^x), 
                                .SDcols = features]
    }
    tmp_fit = medpolish(wide[, features, with = FALSE], na.rm = TRUE, trace.iter = FALSE)
    wide[, newABUNDANCE := tmp_fit$overall + tmp_fit$row]
    tmp_result = wide[, list(LABEL, RUN, newABUNDANCE)]
    if (original_scale) {
        tmpt_result[, newABUNDANCE := log(newABUNDANCE, base_log)]
    }
    
    if (data.table::uniqueN(input$LABEL) == 2) {
        tmp_result = .adjustLRuns(tmp_result)
    }
    tmp_result[, list(RUN, LogIntensities = newABUNDANCE)]
}

.adjustLRuns = function(input, rename = FALSE) {
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

.getNonMissingFilterStats = function(input, censored_symbol) {
    if (!is.null(censored_symbol)) {
        if (censored_symbol == "NA") {
            nonmissing_filter = input$LABEL == "L" & !input$censored
        } else {
            nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE) & !input$censored 
        }
    } else {
        nonmissing_filter = input$LABEL == "L" & !is.na(input$INTENSITY)
    }
    nonmissing_filter = nonmissing_filter & input$n_obs_run > 0 & input$n_obs > 1
    nonmissing_filter
}

