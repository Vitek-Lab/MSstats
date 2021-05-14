#' Check if groupComparison input was processed by the dataProcess function
#' @param input data.table
#' @keywords internal
.checkGroupComparisonInput = function(input) {
    cols = c("RUN", "Protein", "LogIntensities", "originalRUN",
             "GROUP", "SUBJECT", "more50missing", "NumMeasuredFeature")
    
    if (length(setdiff(toupper(cols), toupper(colnames(input)))) > 0) {
        msg = paste("The `data` input was not processed by the dataProcess function.",
                    "Please use the dataProcess function first.")
        getOption("MSstatsLog")("INFO", msg)
        stop(msg)
    }
    input
}


#' Check if contrast matrix includes all conditions
#' @param contrast_matrix contrast matrix
#' @param input data.table of summarized data
#' @keywords internal
.checkContrastMatrix = function(contrast_matrix, input) {
    if (ncol(contrast_matrix) != nlevels(input$GROUP)) {
        msg = paste("Number of columns of the contrast.matrix parameter must be",
                    "equal to the number of groups.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    
    if (any(is.null(row.names(contrast_matrix)))) {
        msg = paste("All rows of the contrast.matrix parameter",
                    "must be named")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    contrast_matrix
}


#' Check if there is only single subject
#' @param input data.table
#' @keywords internal
.checkSingleSubject = function(input) {
    SUBJECT = GROUP = NULL
    
    unique_annot = unique(input[, list(GROUP, SUBJECT)])
    subject_counts = unique_annot[, list(NumSubjects = data.table::uniqueN(SUBJECT)),
                                  by = "GROUP"]
    all(subject_counts$NumSubject == 1)
}


#' Check if there are technical replicates
#' @param input data.table
#' @keywords internal
.checkTechReplicate = function(input) {
    GROUP = RUN = SUBJECT = NULL
    
    unique_annot = unique(input[, list(RUN,
                                       SUBJECT_NESTED = paste(GROUP,
                                                              SUBJECT,
                                                              sep = "."))])
    run_counts = unique_annot[, list(NumRuns = data.table::uniqueN(RUN)),
                              by = "SUBJECT_NESTED"]
    all(run_counts$NumRuns != 1)
}
