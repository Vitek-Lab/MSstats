#' Post-processing output from MSstats summarization
#' 
#' @param input `data.table` in MSstats format
#' @param summarized output of the `MSstatsSummarize` function
#' @param processed output of MSstatsSelectFeatures
#' @param method name of the summarization method
#' (`summaryMethod` parameter to `dataProcess`)
#' @param impute if TRUE, censored missing values were imputed
#' (`MBimpute` parameter to `dataProcess`)
#' @param censored_symbol censored missing value indicator 
#' (`censoredInt` parameter to `dataProcess`)
#' 
#' @return list that consists of the following elements:
#' \itemize{
#' \item{FeatureLevelData}{ - feature-level data after processing} 
#' \item{ProteinLevelData}{ - protein-level (summarized) data}
#' \item{SummaryMethod}{ (string) - name of summarization method that was used}
#' }
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
#' input = MSstatsSelectFeatures(input, "all")
#' processed = getProcessed(input)
#' input = MSstatsPrepareForSummarization(input, method, impute, cens, FALSE)
#' input_split = split(input, input$PROTEIN)
#' summarized = MSstatsSummarize(input_split, method, impute, cens, FALSE, TRUE)
#' output = output = MSstatsSummarizationOutput(input, summarized, processed,
#' method, impute, cens)
#' 
MSstatsSummarizationOutput = function(input, summarized, processed, 
                                      method, impute, censored_symbol) {
    LABEL = TotalGroupMeasurements = GROUP = Protein = RUN = NULL
    
    input = .finalizeInput(input, summarized, method, impute, censored_symbol)
    summarized = lapply(summarized, function(x) x[[1]])
    summarized = data.table::rbindlist(summarized)
    if (inherits(summarized, "try-error")) {
        msg = paste("*** error : can't summarize per subplot with ", 
                    method, ".")
        getOption("MSstatsLog")("ERROR", msg)
        getOption("MSstatsMsg")("ERROR", msg)
        rqall = NULL
        rqmodelqc = NULL
        workpred = NULL
    } else {
        input[LABEL == "L", TotalGroupMeasurements := uniqueN(.SD),
              by = c("PROTEIN", "GROUP"), 
              .SDcols = c("FEATURE", "originalRUN")]
        cols = intersect(c("PROTEIN", "originalRUN", "RUN", "GROUP",
                           "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", 
                           "TotalGroupMeasurements",
                           "NumMeasuredFeature", "MissingPercentage", 
                           "more50missing", "NumImputedFeature"),
                         colnames(input))
        merge_col = ifelse(is.element("RUN", colnames(summarized)), 
                           "RUN", "SUBJECT_ORIGINAL")
        lab = unique(input[LABEL == "L", cols, with = FALSE])
        if (nlevels(input$LABEL) > 1) {
            lab = lab[GROUP != 0]
        }
        lab = lab[, colnames(lab) != "GROUP", with = FALSE]
        rqall = merge(summarized, lab, by.x = c(merge_col, "Protein"),
                      by.y = c(merge_col, "PROTEIN"))
        data.table::setnames(rqall, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL"),
                             c("GROUP", "SUBJECT"), skip_absent = TRUE)
        
        rqall$GROUP = factor(as.character(rqall$GROUP))
        rqall$Protein = factor(rqall$Protein)
        rqmodelqc = summarized$ModelQC
    }
    
    if (is.element("RUN", colnames(rqall)) & !is.null(rqall)) {
        rqall = rqall[order(Protein, as.numeric(as.character(RUN))), ]
        rownames(rqall) = NULL
    }
    output_cols = intersect(c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE",
                              "LABEL", "GROUP", "RUN", "SUBJECT", "FRACTION",
                              "originalRUN", "censored", "INTENSITY", "ABUNDANCE",
                              "newABUNDANCE", "predicted", "feature_quality", 
                              "is_outlier", "remove"), colnames(input))
    input = input[, output_cols, with = FALSE]
    
    if (is.element("remove", colnames(processed))) {
        processed = processed[(remove), 
                              intersect(output_cols, 
                                        colnames(processed)), with = FALSE]
        input = rbind(input, processed, fill = TRUE)
    }
    list(FeatureLevelData = as.data.frame(input), 
         ProteinLevelData = as.data.frame(rqall), 
         SummaryMethod = method)
    
}


#' Add summary statistics to dataProcess output
#' @param input feature-level data
#' @param summarized protein-level data (list)
#' @param method summary method
#' @param impute if TRUE, censored missing values were imputed
#' @param censored_symbol censored missing value indicator
#' @keywords internal
.finalizeInput = function(input, summarized, method, impute, censored_symbol) {
    if (method == "TMP") {
        input = .finalizeTMP(input, censored_symbol, impute, summarized)
    } else {
        input = .finalizeLinear(input, censored_symbol)
    }
    input
}


#' Summary statistics for output of TMP-based summarization
#' @inheritParams .finalizeInput
#' @keywords internal
.finalizeTMP = function(input, censored_symbol, impute, summarized) {
    NonMissingStats = NumMeasuredFeature = MissingPercentage = LABEL = NULL
    total_features = more50missing = nonmissing_orig = censored = NULL
    INTENSITY = newABUNDANCE = NumImputedFeature = NULL
    
    survival_predictions = lapply(summarized, function(x) x[[2]])
    predicted_survival = data.table::rbindlist(survival_predictions)
    if (impute) {
        cols = intersect(colnames(input), c("newABUNDANCE", 
                                            "cen", "RUN",
                                            "FEATURE", "ref"))
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
            input[, NumImputedFeature := sum(LABEL == "L" & !nonmissing_orig),
                  by = c("PROTEIN", "RUN")]
        } else {
            input[, NumImputedFeature := 0]
        }
    }
    input
}


#' Summary statistics for linear model-based summarization
#' @inheritParams .finalizeInput
#' @keywords internal
.finalizeLinear = function(input, censored_symbol) {
    NonMissingStats = NumMeasuredFeature = MissingPercentage = NULL
    total_features = more50missing = nonmissing_orig = LABEL = NULL
    censored = INTENSITY = newABUNDANCE = NumImputedFeature = NULL
    
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
        input[, NumImputedFeature := 0]
    }
    input
}
