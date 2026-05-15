#' Post-processing output from MSstats summarization
#' 
#' @param input `data.table` in MSstats format
#' @param summarized output of the `MSstatsSummarizeWithSingleCore` function
#' @param processed output of MSstatsSelectFeatures
#' @param method name of the summarization method
#' (`summaryMethod` parameter to `dataProcess`)
#' @param impute if TRUE, censored missing values were imputed
#' (`MBimpute` parameter to `dataProcess`)
#' @param censored_symbol censored missing value indicator 
#' (`censoredInt` parameter to `dataProcess`)
#' 
#' @return A list with the following elements:
#'   \item{FeatureLevelData}{Feature-level data after processing.}
#'   \item{ProteinLevelData}{Protein-level (summarized) data.}
#'   \item{SummaryMethod}{String: name of the summarization method used.}
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
#' summarized = MSstatsSummarizeWithSingleCore(input, method, impute, cens, FALSE, TRUE, 100)
#' output = output = MSstatsSummarizationOutput(input, summarized, processed,
#' method, impute, cens)
#' 
MSstatsSummarizationOutput = function(input, summarized, processed,
                                      method, impute, censored_symbol) {
    LABEL = TotalGroupMeasurements = GROUP = Protein = RUN = NULL

    predicted_survival = data.table::rbindlist(lapply(summarized, function(x) x[[2]]),
                                                fill = TRUE)
    for (i in seq_along(summarized)) summarized[[i]][[2]] = NULL
    input = .finalizeInput(input, predicted_survival, method, impute, censored_symbol)
    rm(predicted_survival)
    protein_summaries = lapply(summarized, function(x) x[[1]])
    rm(summarized)
    summarized = data.table::rbindlist(protein_summaries, fill = TRUE)
    rm(protein_summaries)

    if (inherits(summarized, "try-error")) {
        msg = paste("*** error : can't summarize per subplot with ",
                    method, ".")
        getOption("MSstatsLog")("ERROR", msg)
        getOption("MSstatsMsg")("ERROR", msg)
        rqall = NULL
        rqmodelqc = NULL
        workpred = NULL
    } else {
        input[, TotalGroupMeasurements := uniqueN(.SD),
              by = c("PROTEIN", "GROUP", "LABEL"),
              .SDcols = c("FEATURE", "originalRUN")]
        cols = intersect(c("PROTEIN", "originalRUN", "RUN", "GROUP",
                           "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "LABEL",
                           "TotalGroupMeasurements",
                           "NumMeasuredFeature", "MissingPercentage",
                           "more50missing", "NumImputedFeature"),
                         colnames(input))
        merge_col = ifelse(is.element("RUN", colnames(summarized)),
                           "RUN", "SUBJECT_ORIGINAL")
        lab = unique(input[, cols, with = FALSE])
        lab = lab[, colnames(lab) != "GROUP", with = FALSE]
        rqall = merge(summarized, lab, by.x = c(merge_col, "Protein", "LABEL"),
                      by.y = c(merge_col, "PROTEIN", "LABEL"))
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
                              "is_outlier", "remove", "is_labeled_ref"), colnames(input))
    drop_cols = setdiff(colnames(input), output_cols)
    for (col in drop_cols) data.table::set(input, j = col, value = NULL)

    if (is.element("remove", colnames(processed))) {
        processed = processed[(remove),
                              intersect(output_cols,
                                        colnames(processed)), with = FALSE]
        input = rbind(input, processed, fill = TRUE)
    }
    data.table::setDF(input)
    data.table::setDF(rqall)
    list(FeatureLevelData = input,
         ProteinLevelData = rqall,
         SummaryMethod = method)

}


#' Add summary statistics to dataProcess output
#' @param input feature-level data
#' @param summarized protein-level data (list)
#' @param method summary method
#' @param impute if TRUE, censored missing values were imputed
#' @param censored_symbol censored missing value indicator
#' @keywords internal
.finalizeInput = function(input, predicted_survival, method, impute, censored_symbol) {
    # if (method == "TMP") {
    input = .finalizeTMP(input, censored_symbol, impute, predicted_survival)
    # } else {
    #     input = .finalizeLinear(input, censored_symbol)
    # }
    input
}


#' Summary statistics for output of TMP-based summarization
#' @inheritParams .finalizeInput
#' @keywords internal
.finalizeTMP = function(input, censored_symbol, impute, predicted_survival) {
    NonMissingStats = NumMeasuredFeature = MissingPercentage = LABEL = NULL
    total_features = more50missing = nonmissing_orig = censored = NULL
    INTENSITY = newABUNDANCE = NumImputedFeature = NULL

    if (impute) {
        join_cols = intersect(intersect(colnames(input),
                                        colnames(predicted_survival)),
                              c("cen", "RUN", "FEATURE", "ref_covariate",
                                "LABEL"))
        data.table::set(input, j = "newABUNDANCE", value = NULL)
        idx = predicted_survival[input, on = join_cols, which = TRUE,
                                 mult = "first"]
        data.table::set(input, j = "newABUNDANCE",
                        value = predicted_survival$newABUNDANCE[idx])
        data.table::set(input, j = "predicted",
                        value = predicted_survival$predicted[idx])
    }
    input[, NonMissingStats := .getNonMissingFilterStats(.SD, censored_symbol)]
    input[, NumMeasuredFeature := sum(NonMissingStats),
          by = c("PROTEIN", "RUN", "LABEL")]
    input[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
    input[, more50missing := MissingPercentage >= 0.5]
    if (!is.null(censored_symbol)) {
        if (is.element("censored", colnames(input))) {
            input[, nonmissing_orig := !censored]
        } else {
            input[, nonmissing_orig := !is.na(INTENSITY)]
        }
        input[is.na(newABUNDANCE), nonmissing_orig := TRUE]
        if (impute) {
            input[, NumImputedFeature := sum(!nonmissing_orig),
                  by = c("PROTEIN", "RUN", "LABEL")]
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
          by = c("PROTEIN", "RUN", "LABEL")]
    input[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
    input[, more50missing := MissingPercentage >= 0.5]
    if (!is.null(censored_symbol)) {
        if (is.element("censored", colnames(input))) {
            input[, nonmissing_orig := !censored]
        } else {
            input[, nonmissing_orig := !is.na(INTENSITY)]
        }
        input[is.na(newABUNDANCE), nonmissing_orig := TRUE]
        input[, NumImputedFeature := 0]
    }
    input
}
