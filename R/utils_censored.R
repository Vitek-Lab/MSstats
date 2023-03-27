#' Handle censored missing values 
#' 
#' @param input `data.table` in MSstats data format
#' @param summary_method summarization method (`summaryMethod` parameter to `dataProcess`)
#' @param impute if TRUE, missing values are supposed to be imputed 
#' (`MBimpute` parameter to `dataProcess`)
#' @param missing_symbol `censoredInt` parameter to `dataProcess`
#' @param censored_cutoff `maxQuantileforCensored` parameter to `dataProcess`
#' 
#' @importFrom stats quantile
#' 
#' @export
#' 
#' @return data.table
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
#' head(input)
#'  
MSstatsHandleMissing = function(input, summary_method, impute, 
                                missing_symbol, censored_cutoff) {
    INTENSITY = LABEL = ABUNDANCE = censored = NULL
    
    if ((summary_method == "TMP" & impute) & !is.null(missing_symbol)) {
        input$censored = FALSE
        ## if intensity = 1, but abundance > cutoff after normalization, it also should be censored.
        if (!is.null(censored_cutoff)) {
            quantiles = input[!is.na(INTENSITY) & INTENSITY > 1 & LABEL == "L", 
                              quantile(ABUNDANCE, 
                                       prob = c(0.01, 0.25, 0.5, 0.75, 
                                                censored_cutoff), 
                                       na.rm = TRUE)]
            iqr = quantiles[4] - quantiles[2]
            multiplier = (quantiles[5] - quantiles[4]) / iqr
            cutoff_lower = (quantiles[2] - multiplier * iqr) 
            input$censored = !is.na(input$INTENSITY) & 
                input$LABEL == "L" &
                input$ABUNDANCE < cutoff_lower
            if (cutoff_lower <= 0 & !is.null(missing_symbol) & missing_symbol == "0") {
                zero_one_filter = !is.na(input$ABUNDANCE) & input$ABUNDANCE <= 0
                input$censored = ifelse(zero_one_filter, TRUE, input$censored)
            }
            if (!is.null(missing_symbol) & missing_symbol == "NA") {
                input$censored = ifelse(is.na(input$INTENSITY), TRUE, 
                                        input$censored)
            }
            
            msg = paste('** Log2 intensities under cutoff =', 
                        format(cutoff_lower, digits = 5), 
                        ' were considered as censored missing values.')
            msg_2 = paste("** Log2 intensities =", missing_symbol, "were considered as censored missing values.")
            
            getOption("MSstatsMsg")("INFO", msg)
            getOption("MSstatsMsg")("INFO", msg_2)
            
            getOption("MSstatsLog")("INFO", msg)
            getOption("MSstatsLog")("INFO", msg_2)
            
        } else {
            if (missing_symbol == '0') {
                input$censored = input$LABEL == "L" & 
                    !is.na(input$INTENSITY) &
                    (input$INTENSITY == 1 | input$ABUNDANCE <= 0)
            } else if (missing_symbol == 'NA') {
                input$censored = input$LABEL == "L" & is.na(input$ABUNDANCE)
            }
        }
        input[, censored := ifelse(LABEL == "H", FALSE, censored)]
    } else {
        input$censored = FALSE
    }
    input
}


#' Set censored values based on minimum in run/feature/run or feature
#' @param input `data.table` in MSstats format
#' @param censored_symbol censoredInt parameter to `dataProcess`
#' @param remove50missing if TRUE, features with at least 50% missing values
#' will be removed
#' @keywords internal
.setCensoredByThreshold = function(input, censored_symbol, remove50missing) {
  total_features = n_obs = newABUNDANCE = n_obs_run = censored = NULL
  nonmissing_all = ABUNDANCE_cut = NULL
  
  if (censored_symbol == "NA") {
    input[, nonmissing_all := !is.na(newABUNDANCE)]
  } else if (censored_symbol == "0") {
    input[, nonmissing_all := !is.na(newABUNDANCE) & input$newABUNDANCE != 0]
  }
  
  input[, nonmissing_all := ifelse(total_features > 1 & n_obs <= 1, 
                                   FALSE, nonmissing_all)] 
  grouping_vars = c("PROTEIN", "FEATURE", "LABEL")
  input[n_obs > 1 & n_obs_run > 0, 
        ABUNDANCE_cut := .getMin(newABUNDANCE, nonmissing_all),
        by = grouping_vars]
  input[, any_censored := any(censored & n_obs > 1 & n_obs_run > 0),
        by = "PROTEIN"]
  if (censored_symbol == "NA") {
    input[, newABUNDANCE := ifelse(!nonmissing_all & censored & is.finite(ABUNDANCE_cut) & any_censored, 
                                   ABUNDANCE_cut, newABUNDANCE)]
  } else if (censored_symbol == "0") {
    input[, newABUNDANCE := ifelse(!nonmissing_all & newABUNDANCE == 0 & is.finite(ABUNDANCE_cut) & any_censored, 
                                   ABUNDANCE_cut, newABUNDANCE)]
  }
}


#' Utility function: get 0.99 * minimum of non-missing values
#' @param abundance abundances values
#' @param nonmissing logical vector
#' @keywords internal
.getMin = function(abundance, nonmissing) {
    0.99 * min(abundance[nonmissing], na.rm = TRUE)
}


#' Identify non-missing values
#' @param input `data.table` in MSstats format
#' @param impute if TRUE, missing values are supposed to be imputed
#' @param censored_symbol `censoredInt` parameter to dataProcess
#' @keywords internal
.getNonMissingFilter = function(input, impute, censored_symbol) {
    if (impute) {
        if (!is.null(censored_symbol)) {
            if (censored_symbol == "0") {
                nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE) & input$newABUNDANCE != 0
            } else if (censored_symbol == "NA") {
                nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE)
            }  
        } 
    } else {
        nonmissing_filter = input$LABEL == "L" & !is.na(input$newABUNDANCE) & input$newABUNDANCE != 0
    }
    nonmissing_filter
}
