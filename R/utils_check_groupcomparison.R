.checkGroupComparisonInput = function(input) {
    cols = c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE",
             "LABEL", "GROUP", "RUN", "SUBJECT", "FRACTION", "originalRUN",
             "INTENSITY", "ABUNDANCE")
    
    if (length(setdiff(toupper(cols), toupper(colnames(input$ProcessedData)))) > 0) {
        msg = paste("The `data` input was not processed by the dataProcess function.",
                    "Please use the dataProcess function first.")
        getOption("MSstatsLog")("INFO", msg)
        stop(msg)
    }
    input
}


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


.checkRepeated = function(input) {
    subject_by_group = table(input[LABEL == "L", list(SUBJECT, GROUP)])
    subject_appearances = apply(subject_by_group, 1, function(x) sum(x > 0))
    repeated = any(subject_appearances > 1)
    if (repeated) {
        msg = "Time course design of experiment - okay"
    } else {
        msg = "Case control design of experiment - okay"
    }
    getOption("MSstatsLog")("INFO", msg)
    repeated
}


.checkSingleSubject = function(input) {
    SUBJECT = NULL
    
    unique_annot = unique(input[, list(GROUP, SUBJECT)])
    subject_counts = unique_annot[, list(NumSubjects = data.table::uniqueN(SUBJECT)),
                                  by = "GROUP"]
    all(subject_counts$NumSubject == 1)
}


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


.getLogBaseName = function(processed) {
    tmp = head(processed[!is.na(ABUNDANCE) & !is.na(INTENSITY) & ABUNDANCE > 2], 
               1)
    log_2_diff = abs(log(tmp$INTENSITY, 2) - tmp$ABUNDANCE)
    log_10_diff = abs(log(tmp$INTENSITY, 10) - tmp$ABUNDANCE)
    if (log_2_diff < log_10_diff) {
        "log2FC"
    } else {
        "log10FC"
    }
}
