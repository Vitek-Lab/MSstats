#' Post-processing output from MSstats summarization
#' 
#' @param input `data.table` in MSstats format
#' @param summarized output of the `MSstatsSummarize` function
#' @param summary_method name of the summarization method
#' (`summaryMethod` parameter to `dataProcess`)
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
        cols = intersect(c("PROTEIN", "originalRUN", "RUN", "GROUP_ORIGINAL", 
                           "SUBJECT_ORIGINAL", "NumMeasuredFeature", 
                           "MissingPercentage", "more50missing", 
                           "NumImputedFeature"),
                         colnames(input))
        merge_col = ifelse(is.element("RUN", colnames(summarized)), "RUN", "SUBJECT_ORIGINAL")
        lab = unique(input[, cols, with = FALSE])
        if (nlevels(input$LABEL) > 1) {
            lab = lab[GROUP != 0]
        }
        rqall = merge(summarized, lab, by.x = c(merge_col, "Protein"),
                      by.y = c(merge_col, "PROTEIN"))
        data.table::setnames(rqall, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL"),
                             c("GROUP", "SUBJECT"))
        
        rqall$GROUP = factor(rqall$GROUP)
        rqall$Protein = factor(rqall$Protein)
        rqmodelqc = summarized$ModelQC
        workpred <- NULL
    }
    
    if (is.element("RUN", colnames(rqall))) {
        rqall = rqall[order(Protein, as.numeric(as.character(RUN))), ]
        rownames(rqall) = NULL
    }
    
    output_cols = intersect(c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE",
                    "LABEL", "GROUP", "RUN", "SUBJECT", "FRACTION",
                    "originalRUN", "censored", "INTENSITY", "ABUNDANCE",
                    "newABUNDANCE", "predicted"), colnames(input))
    list(ProcessedData = as.data.frame(input)[, output_cols], 
         RunlevelData = as.data.frame(rqall), 
         SummaryMethod = summary_method, 
         ModelQC = NULL, 
         PredictBySurvival = NULL)
}
