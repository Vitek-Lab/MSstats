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
#' @return data.table
#' @keywords internal
.runTukey = function(input, is_labeled, censored_symbol, remove50missing) {
    Protein = RUN = newABUNDANCE = NULL
    
    if (nlevels(input$FEATURE) > 1) {
        tmp_result = .fitTukey(input)
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
.fitTukey = function(input) {
    LABEL = RUN = newABUNDANCE = NULL
    
    features = as.character(unique(input$FEATURE))
    wide = data.table::dcast(LABEL + RUN ~ FEATURE, data = input,
                             value.var = "newABUNDANCE", keep = TRUE)
    tmp_fitted = median_polish_summary(as.matrix(wide[, features, with = FALSE]))
    wide[, newABUNDANCE := tmp_fitted]
    tmp_result = wide[, list(LABEL, RUN, newABUNDANCE)]

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
