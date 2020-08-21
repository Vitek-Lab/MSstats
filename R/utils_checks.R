#' Save information about R session to sessionInfo.txt file.
#' @importFrom utils sessionInfo
#' @keywords internal
.saveSessionInfo = function() {
    session_info = utils::sessionInfo()
    sink("sessionInfo.txt")
    print(session_info)
    sink()
}

#' Check validity of parameters to dataProcess function
#' @param log_base of logarithmic transformation
#' @param normalization_method string: "quantile", "equalizemedians", "FALSE",
#' "NONE" or "globalStandards"
#' @param address string
#' @param fill_rows logical, if TRUE, missing run observations for each feature
#' will be added with INTENSITY = NA
#' @param feature_selection list with elements: remove_uninformative
#' @param summarization list with elements: method.
#' @param imputation list with elements: cutoff, symbol.
#' @param n_clusters integer
#' @keywords internal
.checkDataProcessParams = function(log_base, normalization_method,
                                   standards_names, address, fill_rows,
                                   feature_selection, summarization,
                                   imputation, n_clusters) {
    checkmate::assertChoice(log_base, c(2, 10), .var.name = "logTrans")
    checkmate::assertLogical(fill_rows, .var.name = "fillIncompleteRows")
    checkmate::assertChoice(summarization$method, c("linear", "TMP"),
                            .var.name = "summaryMethod") 
    getOption("MSstatsLog")("INFO", paste("Summary method:", 
                                          summarization$method))
    checkmate::assertChoice(imputation$cutoff, c("minFeature", "minRun", 
                                                 "minFeatureNRun"),
                            .var.name = "cutoffCensored")
    getOption("MSstatsLog")("INFO", paste("cutOffCensored:", imputation$cutoff))
    checkmate::assertChoice(imputation$symbol, c("0", "NA"), 
                            null.ok = TRUE, .var.name = "censoredInt")
    getOption("MSstatsLog")("INFO", paste("censoredInt:", imputation$symbol))
    ## check input for censoredInt and MBimpute
    if (summarization$method == "TMP" & imputation$MB & is.null(imputation$symbol)) {
        msg = paste("The combination of required input",
                    "MBimpute == TRUE and censoredInt = NULL",
                    "has no censore missing values.",
                    "Imputation will not be performed.- stop")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    checkmate::assertChoice(toupper(as.character(normalization_method)),
                            c("NONE", "FALSE", "EQUALIZEMEDIANS", "QUANTILE", "GLOBALSTANDARDS"),
                            .var.name = "normalization")
    if (toupper(as.character(normalization_method)) == "GLOBALSTANDARDS" &
        is.null(standards_names)) {
        msg = paste("For normalization with global standards,",
                    "the names of global standards are needed.",
                    "Please add 'nameStandards' input.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    # TODO: more checks
}


#' Check if a given column exists in the data
#' @param input data.table
#' @param column_name chr, name of a column to check
#' @keywords internal
.checkExperimentDesign = function(input, column_name) {
    if (any(is.na(input[[column_name]]))) {
        msg = paste("Missing information in the", column_name, "column. ",
                    "Please check the", column_name, "column.", collapse = " ")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
}


.updateColnames = utils::getFromNamespace(".updateColnames", "MSstatsConvert")
#' Check validity of data that were not processed by MSstats converter
#' @param input data.table
#' @importFrom data.table uniqueN as.data.table
#' @keywords internal
.checkUnProcessedDataValidity = function(input) {
    input = data.table::as.data.table(input)
    colnames(input) = toupper(colnames(input))
    cols = c("ProteinName", "PeptideSequence", "PeptideModifiedSequence",
             "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", 
             "Condition", "BioReplicate", "Run", "Intensity")
    cols = toupper(cols)
    provided_cols = intersect(cols, colnames(input))
    
    if (length(provided_cols) < 10) {
        msg = paste("Missing columns in the input:", 
                    paste(setdiff(cols, colnames(input)), collapse = " "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }    
    
    if (is.element("PEPTIDEMODIFIEDSEQUENCE", provided_cols)) {
        colnames(input) = .updateColnames(
            input, "PEPTIDEMODIFIEDSEQUENCE", "PEPTIDESEQUENCE")
    }
    
    if (!is.numeric(input$INTENSITY)) {	
        suppressWarnings({
            input$INTENSITY = as.numeric(as.character(input$INTENSITY))
        })
    }
    
    .checkExperimentDesign(input, "RUN")
    .checkExperimentDesign(input, "BIOREPLICATE")
    .checkExperimentDesign(input, "CONDITION")
    
    cols = intersect(c(cols, "FRACTION", "TECHREPLICATE"),
                     colnames(input))
    input = input[, cols, with = FALSE]
    input$PEPTIDE = paste(input$PEPTIDESEQUENCE, input$PRECURSORCHARGE, sep = "_")
    input$TRANSITION = paste(input$FRAGMENTION, input$PRODUCTCHARGE, sep = "_")
    # TODO: := ?
    
    if (data.table::uniqueN(input$ISOTOPELABELTYPE) > 2) {
        getOption("MSstatsLog")("ERROR",  paste("There are more than two levels of labeling.",
                                                "So far, only label-free or reference-labeled experiment are supported. - stop"))
        stop("Statistical tools in MSstats are only proper for label-free or with reference peptide experiments.")
    }
    
    input$ISOTOPELABELTYPE = factor(input$ISOTOPELABELTYPE)
    if (data.table::uniqueN(input$ISOTOPELABELTYPE) == 2) {
        levels(input$ISOTOPELABELTYPE) = c("H", "L")
    } else {
        levels(input$ISOTOPELABELTYPE) = "L"
    }
    input
}


#' Check validity of data already processed by MSstats converter
#' @param input data.frame of class `MSstatsValidated`
#' @keywords internal
.prepareForDataProcess = function(input) {
    colnames(input) = toupper(colnames(input))
    
    if (is.element("PEPTIDEMODIFIEDSEQUENCE", colnames(input))) {
        colnames(input) = MSstatsConvert:::.updateColnames(
            input, "PEPTIDEMODIFIEDSEQUENCE", "PEPTIDESEQUENCE")
    }
    
    input$PEPTIDE = paste(input$PEPTIDESEQUENCE, input$PRECURSORCHARGE, sep = "_")
    input$TRANSITION = paste(input$FRAGMENTION, input$PRODUCTCHARGE, sep = "_")
    # TODO: := ?
    
    input$ISOTOPELABELTYPE = factor(input$ISOTOPELABELTYPE)
    if (data.table::uniqueN(input$ISOTOPELABELTYPE) == 2) {
        levels(input$ISOTOPELABELTYPE) = c("H", "L")
    } else {
        levels(input$ISOTOPELABELTYPE) = "L"
    }
    input
}

setGeneric(".checkDataValidity", 
           function(input) standardGeneric(".checkDataValidity"))
setMethod(".checkDataValidity", "data.frame", .checkUnProcessedDataValidity)
setMethod(".checkDataValidity", "MSstatsValidated", .prepareForDataProcess) # TODO: or change column names


#' Create ABUNDANCE column and log-transform intensities
#' @param input data.table
#' @param log_base base of the logarithm
#' @keywords internal
.preProcessIntensities = function(input, log_base) {
    n_smaller_than_1 = sum(input$INTENSITY < 1, na.rm = TRUE)
    if (n_smaller_than_1 > 0) {
        input[INTENSITY < 1, INTENSITY := 1]
        msg = paste("** There are", n_smaller_than_1, 
                    "intensities which are zero or less than 1.",
                    "These intensities are replaced with 1",
                    collapse = " ")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    } 
    input[, ABUNDANCE := log(INTENSITY, log_base)]
    getOption("MSstatsLog")("INFO",
                            paste("Logarithm transformation with base",
                                  log_base,
                                  "is done - okay",
                                  collapse = " "))
}


#' Create columns for data processing
#' @param input data.table
#' @keywords internal
.updateColumnsForProcessing = function(input) {
    input[, FEATURE := paste(PEPTIDE, TRANSITION, sep = "_")]
    input[, GROUP := ifelse(LABEL == "L", GROUP_ORIGINAL, 0)]
    input[, SUBJECT := ifelse(LABEL == "L", SUBJECT_ORIGINAL, 0)]
    input[, SUBJECT_NESTED := paste(GROUP, SUBJECT, sep = ".")]
    # as.factor() / factor() everywhere?
    
    cols = c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE", "LABEL", 
             "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "RUN", "GROUP", "SUBJECT", 
             "SUBJECT_NESTED", "INTENSITY")
    input[!is.na(PROTEIN) & PROTEIN != "", cols, with = FALSE]
    # processout <- rbind(processout, c("New input format : made new columns for analysis - okay"))
    # write.table(processout, file=finalfile, row.names=FALSE)
}


#' Make factor columns where needed
#' @param input data.table
#' @keywords internal
.makeFactorColumns = function(input) {
    input$PROTEIN = factor(input$PROTEIN)
    input$PEPTIDE = factor(input$PEPTIDE)
    input$TRANSITION = factor(input$TRANSITION)
    input = input[order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL, 
                        RUN, PROTEIN, PEPTIDE, TRANSITION), ]
    input$GROUP = factor(input$GROUP)
    input$SUBJECT = factor(input$SUBJECT)
    input$SUBJECT_NESTED = factor(input$SUBJECT_NESTED, levels = unique(input$SUBJECT_NESTED))
    input$FEATURE = factor(input$FEATURE, levels = unique(input$FEATURE))
    input$originalRUN = input$RUN
    input$RUN = factor(input$RUN, levels = unique(input$RUN), labels = seq(1, length(unique(input$RUN))))
    msg = paste("Factorize in columns(GROUP, SUBJECT, GROUP_ORIGINAL,",
                "SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN) - okay")
    getOption("MSstatsLog")("INFO", msg)
    input
}
