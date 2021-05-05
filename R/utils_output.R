#' Post-processing output from MSstats summarization
#' 
#' @param input `data.table` in MSstats format
#' @param summarized output of the `MSstatsSummarize` function
#' @param summary_method name of the summarization method
#' (`summaryMethod` parameter to `dataProcess`)
#' 
#' @return list that consists of the following elements:
#' \itemize{
#' \item{ProcessedData}{ - feature-level data after processing} 
#' \item{RunlevelData}{ - run-level (summarized) data}
#' \item{SummaryMethod}{ (string) - name of summarization method that was used}
#' }
#' 
#' @export
#' 
MSstatsSummarizationOutput = function(input, summarized, summary_method) {
    GROUP = Protein = RUN = NULL
    
    if (inherits(summarized, "try-error")) {
        msg = paste("*** error : can't summarize per subplot with ", 
                    summary_method, ".")
        getOption("MSstatsLog")("ERROR", msg)
        getOption("MSstatsMsg")("ERROR", msg)
        rqall = NULL
        rqmodelqc = NULL
        workpred = NULL
    } else {
        input[LABEL == "L", TotalGroupMeasurements := uniqueN(.SD),
              by = c("PROTEIN", "GROUP"), .SDcols = c("FEATURE", "originalRUN")]
        cols = intersect(c("PROTEIN", "originalRUN", "RUN", "GROUP", "GROUP_ORIGINAL", 
                           "SUBJECT_ORIGINAL", "TotalGroupMeasurements",
                           "NumMeasuredFeature", "MissingPercentage", 
                           "more50missing", "NumImputedFeature"),
                         colnames(input))
        merge_col = ifelse(is.element("RUN", colnames(summarized)), "RUN", "SUBJECT_ORIGINAL")
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
        workpred <- NULL
    }
    
    if (is.element("RUN", colnames(rqall)) & !is.null(rqall)) {
        rqall = rqall[order(Protein, as.numeric(as.character(RUN))), ]
        rownames(rqall) = NULL
    }
    
    output_cols = intersect(c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE",
                    "LABEL", "GROUP", "RUN", "SUBJECT", "FRACTION",
                    "originalRUN", "censored", "INTENSITY", "ABUNDANCE",
                    "newABUNDANCE", "predicted", "feature_quality", "is_outlier"), colnames(input))
    list(ProcessedData = as.data.frame(input)[, output_cols], 
         RunlevelData = as.data.frame(rqall), 
         SummaryMethod = summary_method, 
         ModelQC = NULL, 
         PredictBySurvival = NULL)
}
