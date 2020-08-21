.processFinalOutput = function(input, rqresult, summary_method) {
    if (inherits(rqresult, "try-error")) {
        msg = paste("*** error : can't summarize per subplot with ", 
                    summary_method, ".")
        getOption("MSstatsLog")("ERROR", msg)
        getOption("MSstatsMsg")("ERROR", msg)
        rqall = NULL
        rqmodelqc = NULL
        workpred = NULL
    } else {
        label <- nlevels(input$LABEL) == 2
        if (!is.element("RUN", colnames(rqresult))) {
            lab = unique(input[, c("GROUP", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT")])
            if (label) {
                lab = lab[GROUP != 0, ]
            }
            rqall = merge(rqresult, lab, by = "SUBJECT_ORIGINAL")
        } else {
            lab = unique(input[, c("RUN", "originalRUN", "GROUP", "GROUP_ORIGINAL", 
                                   "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT")])
            if (label) {
                lab = lab[lab$GROUP != 0, ]
            }
            rqall = merge(rqresult, lab, by = "RUN")
        }
        
        rqall$GROUP = factor(rqall$GROUP)
        rqall$Protein = factor(rqall$Protein)
        rqmodelqc = rqresult$ModelQC
        #MC : can't use this predicted value.
        #workpred <- rqresult$PredictedBySurvival
        workpred <- NULL
        msg = " == the summarization per subplot is done."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    if (is.element("RUN", colnames(rqall))) {
        rqall = rqall[order(Protein, as.numeric(as.character(RUN))), ]
        rownames(rqall) = NULL
    }
    
    list(ProcessedData = as.data.frame(input), 
         RunlevelData = as.data.frame(rqall), 
         SummaryMethod = summary_method, 
         ModelQC = NULL, 
         PredictBySurvival = NULL)
}