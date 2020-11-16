#' Summarization with Tukey median polish for a full dataset
#' @param input `data.table` processed by the `MSstatsSelectFeatures` function.
#' @param impute only for summaryMethod = "TMP" and censored_symbol = 'NA' or '0'. 
#' TRUE (default) imputes "NA" or "0" (depending on censored_symbol option) by 
#' Accelated failure model. FALSE uses the values assigned by cutoffCensored.
#' @param cutoff_base cutoff value for censoring (with censored_symbol = "NA"
#' or censored_symbol = "0"). "minFeature"/"minRun"/"minFeatureNRun".
#' @param censored_symbol Missing values are censored or at random. 
#' "NA" (default) assumes that all NAs in Intensity column are censored. 
#' "0" uses zero intensities as censored intensity. In this case, NA intensities 
#' are missing at random. The output from Skyline should use "0". 
#' Null assumes that all NA intensites are randomly missing.
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the runs 
#' which have more than 50% missing values. FALSE is default.
#' @param original_scale if TRUE, data will be summarized on the original scale.
#' @param n_threads currently ignored. In the future, number of threads used
#' for parallel processing.
#' @importFrom utils setTxtProgressBar
#' @return data.table
#' @keywords internal
.summarizeTukey = function(input, impute, cutoff_base, censored_symbol, 
                           remove50missing, original_scale = FALSE, 
                           n_threads = NULL) {
    cen = censored = ABUNDANCE = FEATURE = LABEL = more50missing = NULL
    INTENSITY = PROTEIN = n_obs = n_obs_run = RUN = NULL
    
    if (impute & !is.null(censored_symbol)) {
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
        .setCensoredByThreshold(input, cutoff_base, censored_symbol, remove50missing)
    }
    
    proteins = unique(input$PROTEIN)
    n_proteins = length(proteins)
    cols = intersect(colnames(input), c("newABUNDANCE", "cen", "RUN",
                                        "FEATURE", "ref"))

    summarized_results = vector("list", n_proteins)
    survival_predictions = vector("list", n_proteins)
    pb = utils::txtProgressBar(min = 0, max = n_proteins, style = 3)
    for (protein_id in seq_along(proteins)) {
        single_protein = input[PROTEIN == proteins[protein_id]]
        single_protein = single_protein[(n_obs > 1 & !is.na(n_obs)) &
                                            (n_obs_run > 0 & !is.na(n_obs_run))]
        single_protein[, RUN := factor(RUN)]
        single_protein[, FEATURE := factor(FEATURE)]
        if (impute) {
            survival_fit = .fitSurvival(single_protein[LABEL == "L", cols,
                                                       with = FALSE])
            single_protein[, predicted := predict(survival_fit,
                                                  newdata = .SD)]
            single_protein[, predicted := ifelse(censored & (LABEL == "L"), predicted, NA)]
            single_protein[, newABUNDANCE := ifelse(censored & LABEL == "L",
                                           predicted, newABUNDANCE)]
            survival_predictions[[protein_id]] = single_protein[, c(cols, "predicted"),
                                                                with = FALSE]
        }
        
        summarized_results[[protein_id]] = .summarizeTukeySingleProtein(
            single_protein[, list(PROTEIN, LABEL, RUN, FEATURE, newABUNDANCE,
                                  n_obs, n_obs_run, prop_features)],
            impute, cutoff_base, censored_symbol,
            remove50missing, original_scale, n_threads)
        setTxtProgressBar(pb, protein_id)
    }
    predicted_survival = data.table::rbindlist(survival_predictions)
    if (impute) {
        input = merge(input[, colnames(input) != "newABUNDANCE", with = FALSE], 
                      predicted_survival,
                      by = setdiff(cols, "newABUNDANCE"),
                      all.x = TRUE)
    }
    input[, NonMissingStats := .getNonMissingFilterStats(.SD, censored_symbol)]
    input[, NumMeasuredFeature := sum(NonMissingStats), 
          by = c("PROTEIN", "RUN")]
    input[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
    input[, more50missing := MissingPercentage >= 0.5]
    if (!is.null(censored_symbol)) {
        if (is.element("censored", colnames(input))) {
            input[, nonmissing_orig := LABEL == "L" & !censored]
        } else {
            input[, nonmissing_orig := LABEL == "L" & !is.na(INTENSITY)]
        }
        input[, nonmissing_orig := ifelse(is.na(newABUNDANCE), TRUE, nonmissing_orig)]
        if (impute) {
            input[, NumImputedFeature := sum(!nonmissing_orig),
                  by = c("PROTEIN", "RUN")]
        } else {
            input[, NumImputedFeature := 0]
        }
    }
    
    summarized_results = data.table::rbindlist(summarized_results)
    list(input, summarized_results)
}

#' Tukey median polish for a single protein
#' @inheritParams .summarizeTukey
#' @return data.table
#' @keywords internal
.summarizeTukeySingleProtein = function(
    input, impute, cutoff_base, censored_symbol, remove50missing, original_scale,
    n_threads = NULL) {
    newABUNDANCE = RUN = FEATURE = NULL
    
    input = .isSummarizable(input, remove50missing)
    if (is.null(input)) {
        return(NULL)
    } else {
        input[, RUN := factor(RUN)]
        input[, FEATURE := factor(FEATURE)]
        input = input[!is.na(newABUNDANCE), ]
        is_labeled = nlevels(input$LABEL) > 1
        result = .runTukey(input, is_labeled, censored_symbol, remove50missing,
                           original_scale, log_base = 2)
        result
    }
}


.isSummarizable = function(input, remove50missing) {
    n_obs_run = RUN = NULL
    
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
    input
}

#' Fit Tukey median polish
#' @param input data.table with data for a single protein
#' @param is_labeled logical, if TRUE, data is coming from an SRM experiment
#' @inheritParams .summarizeTukey
#' @param log_base base of the logarithm function used for ABUNDANCE column
#' @return data.table
#' @keywords internal
.runTukey = function(input, is_labeled, censored_symbol, remove50missing,
                     original_scale = FALSE, log_base = 2) {
    features = as.character(unique(input$FEATURE))
    Protein = RUN = newABUNDANCE = NULL
    
    if (nlevels(input$FEATURE) > 1) {
        tmp_result = .fitTukey(input, original_scale, log_base)
    } else { 
        if (is_labeled) {
            tmp_result = .adjustLRuns(input, TRUE)
        } else {
            tmp_result = input[input$LABEL == "L", 
                               list(RUN, LogIntensities = newABUNDANCE)]
        }
    }
    tmp_result[, Protein := unique(input$PROTEIN)]
    tmp_result
}


#' Fit tukey median polish for a data matrix
#' @inheritParams .runTukey
#' @return data.table
#' @keywords internal
.fitTukey = function(input, original_scale, log_base) {
    LABEL = RUN = newABUNDANCE = base_log = NULL
    
    features = as.character(unique(input$FEATURE))
    wide = data.table::dcast(LABEL + RUN ~ FEATURE, data = input,
                             value.var = "newABUNDANCE", keep = TRUE)
    if (original_scale) {
        wide[, features] = wide[, lapply(.SD, function(x) log_base^x), 
                                .SDcols = features]
    }
    tmp_fitted = median_polish_summary(as.matrix(wide[, features, with = FALSE]))
    wide[, newABUNDANCE := tmp_fitted]
    tmp_result = wide[, list(LABEL, RUN, newABUNDANCE)]
    if (original_scale) {
        tmp_result[, newABUNDANCE := log(newABUNDANCE, base_log)]
    }
    
    if (data.table::uniqueN(input$LABEL) == 2) {
        tmp_result = .adjustLRuns(tmp_result)
    }
    tmp_result[, list(RUN, LogIntensities = newABUNDANCE)]
}


#' Adjust summarized abundance based on the heavy channel
#' @param input data.table
#' @param rename if TRUE, rename the output column to LogIntensities
#' @return data.table
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
        if (censored_symbol == "NA") {
            nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE) & !input$censored
        } else {
            nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE) & !input$censored 
        }
    } else {
        nonmissing_filter = input$LABEL == "L" & !is.na(input$INTENSITY)
    }
    nonmissing_filter = nonmissing_filter & input$n_obs_run > 0 & input$n_obs > 1
    nonmissing_filter
}
