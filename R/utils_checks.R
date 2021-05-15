#' Prepare data for processing by `dataProcess` function
#' 
#' @param input `data.table` in MSstats format
#' @param log_base base of the logarithm to transform intensities
#' @inheritParams MSstatsConvert::MSstatsBalancedDesign
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
#' head(input)
#' 
MSstatsPrepareForDataProcess = function(input, log_base, fix_missing) {
    input = .checkDataValidity(input, fix_missing = fix_missing)
    input = .updateColumnsForProcessing(input)
    .preProcessIntensities(input, log_base)
    input = .makeFactorColumns(input)
    input
}


#' Save information about R session to sessionInfo.txt file.
#' @importFrom utils sessionInfo
#' @keywords internal
.saveSessionInfo = function() {
    file_name = paste0("msstats_sessionInfo_", 
                       gsub(" " , "T", as.character(Sys.time())), 
                       ".txt")
    file_name = gsub(":", "_", file_name, fixed = TRUE)
    session_info = utils::sessionInfo()
    sink(file_name)
    print(session_info)
    sink()
}

#' Check validity of parameters to dataProcess function
#' @param log_base of logarithmic transformation
#' @param normalization_method string: "quantile", "equalizemedians", "FALSE",
#' "NONE" or "globalStandards"
#' @param feature_selection list with elements: remove_uninformative
#' @param summarization list with elements: method.
#' @param imputation list with elements: cutoff, symbol.
#' @keywords internal
.checkDataProcessParams = function(log_base, normalization_method,
                                   standards_names, feature_selection, 
                                   summarization, imputation) {
    checkmate::assertChoice(log_base, c(2, 10), .var.name = "logTrans")
    checkmate::assertChoice(summarization$method, c("linear", "TMP"),
                            .var.name = "summaryMethod") 
    getOption("MSstatsLog")("INFO", paste("Summary method:", 
                                          summarization$method))
    checkmate::assertChoice(imputation$symbol, c("0", "NA"), 
                            null.ok = TRUE, .var.name = "censoredInt")
    getOption("MSstatsLog")("INFO", paste("censoredInt:", imputation$symbol))
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


#' Check validity of data that were not processed by MSstats converter
#' @param input data.table
#' @inheritParams MSstatsPrepareForDataProcess
#' @importFrom data.table uniqueN as.data.table
#' @keywords internal
.checkUnProcessedDataValidity = function(input, fix_missing, fill_incomplete) {
    input = data.table::as.data.table(unclass(input))
    cols = c("ProteinName", "PeptideSequence", "PeptideModifiedSequence",
             "PrecursorCharge", "FragmentIon", "ProductCharge", 
             "IsotopeLabelType", "Condition", "BioReplicate", "Run", "Intensity")
    provided_cols = intersect(cols, colnames(input))
    
    if (length(provided_cols) < 10) {
        msg = paste("Missing columns in the input:", 
                    paste(setdiff(cols, colnames(input)), collapse = " "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    data.table::setnames(input, "PeptideModifiedSequence", "PeptideSequence", 
                         skip_absent = TRUE)
    
    balanced_cols = c("PeptideSequence", "PrecursorCharge", 
                      "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsBalancedDesign(
        input, balanced_cols, TRUE, TRUE, fix_missing)
    input = data.table::as.data.table(unclass(input))
    data.table::setnames(input, colnames(input), toupper(colnames(input)))
    
    
    if (!is.numeric(input$INTENSITY)) {	
        suppressWarnings({
            input$INTENSITY = as.numeric(as.character(input$INTENSITY))
        })
    }
    
    .checkExperimentDesign(input, "RUN")
    .checkExperimentDesign(input, "BIOREPLICATE")
    .checkExperimentDesign(input, "CONDITION")
    
    cols = toupper(cols)
    cols = intersect(c(cols, "FRACTION", "TECHREPLICATE"),
                     colnames(input))
    input = input[, cols, with = FALSE]
    
    input$PEPTIDE = paste(input$PEPTIDESEQUENCE, input$PRECURSORCHARGE, sep = "_")
    input$TRANSITION = paste(input$FRAGMENTION, input$PRODUCTCHARGE, sep = "_")
    
    if (data.table::uniqueN(input$ISOTOPELABELTYPE) > 2) {
        getOption("MSstatsLog")("ERROR",  
                                paste("There are more than two levels of labeling.",
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
#' @param .. additional parameters, currently ignored
#' @keywords internal
.prepareForDataProcess = function(input, ...) {
    input = as.data.table(unclass(input))
    colnames(input) = toupper(colnames(input))
    if (is.element("PEPTIDEMODIFIEDSEQUENCE", colnames(input))) {
        data.table::setnames(
            input, "PEPTIDEMODIFIEDSEQUENCE", "PEPTIDESEQUENCE")
    }
    input$PEPTIDE = paste(input$PEPTIDESEQUENCE, input$PRECURSORCHARGE, sep = "_")
    input$TRANSITION = paste(input$FRAGMENTION, input$PRODUCTCHARGE, sep = "_")
    input$ISOTOPELABELTYPE = factor(input$ISOTOPELABELTYPE)
    if (data.table::uniqueN(input$ISOTOPELABELTYPE) == 2) {
        levels(input$ISOTOPELABELTYPE) = c("H", "L")
    } else {
        levels(input$ISOTOPELABELTYPE) = "L"
    }
    input
}

setGeneric(".checkDataValidity", 
           function(input, ...) standardGeneric(".checkDataValidity"))
setMethod(".checkDataValidity", "data.frame", .checkUnProcessedDataValidity)
setMethod(".checkDataValidity", "MSstatsValidated", .prepareForDataProcess)


#' Create ABUNDANCE column and log-transform intensities
#' @param input data.table
#' @param log_base base of the logarithm
#' @keywords internal
.preProcessIntensities = function(input, log_base) {
    INTENSITY = ABUNDANCE = NULL
    
    if (any(!is.na(input$INTENSITY) & input$INTENSITY < 1, na.rm = TRUE)) {
        n_smaller_than_1 = sum(!is.na(input$INTENSITY) & input$INTENSITY < 1, 
                               na.rm = TRUE)
        input[, INTENSITY := ifelse(!is.na(INTENSITY) & INTENSITY < 1, 
                                    1, INTENSITY)]
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
                                  "is done",
                                  collapse = " "))
}


#' Create columns for data processing
#' @param input data.table
#' @keywords internal
.updateColumnsForProcessing = function(input) {
    FEATURE = PEPTIDE = TRANSITION = GROUP = LABEL = GROUP_ORIGINAL = NULL
    SUBJECT = SUBJECT_ORIGINAL = PROTEIN = NULL
    
    data.table::setnames(
        input, c("PROTEINNAME", "ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE"), 
        c("PROTEIN", "LABEL", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL"),
        skip_absent = TRUE)
    
    input[, FEATURE := paste(PEPTIDE, TRANSITION, sep = "_")]
    input[, GROUP := ifelse(LABEL == "L", GROUP_ORIGINAL, "0")]
    input[, SUBJECT := ifelse(LABEL == "L", SUBJECT_ORIGINAL, "0")]

    cols = c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE", "LABEL", 
             "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "RUN", "GROUP", 
             "SUBJECT", "FRACTION", "INTENSITY")
    input[!is.na(PROTEIN) & PROTEIN != "", cols, with = FALSE]
}


#' Make factor columns where needed
#' @param input data.table
#' @keywords internal
.makeFactorColumns = function(input) {
    PROTEIN = PEPTIDE = TRANSITION = LABEL = GROUP_ORIGINAL = RUN = GROUP = NULL
    SUBJECT_ORIGINAL = FEATURE = originalRUN = SUBJECT = NULL
    
    input[, PROTEIN := factor(PROTEIN)]
    input[, PEPTIDE := factor(PEPTIDE)]
    input[, TRANSITION := factor(TRANSITION)]
    input = input[order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL,
                        RUN, PROTEIN, PEPTIDE, TRANSITION), ]
    input[, GROUP := factor(GROUP)]
    input[, SUBJECT := factor(SUBJECT)]
    input[, FEATURE := factor(FEATURE)]
    input[, originalRUN := factor(as.character(RUN))]
    input[, RUN := factor(RUN, levels = unique(RUN),
                          labels = seq_along(unique(RUN)))]
    
    msg = paste("Factorize in columns(GROUP, SUBJECT, GROUP_ORIGINAL,",
                "SUBJECT_ORIGINAL, FEATURE, RUN)")
    getOption("MSstatsLog")("INFO", msg)
    input
}


#' Prepare a peptides dictionary for global standards normalization
#' 
#' @param input `data.table` in MSstats standard format
#' @param normalization normalization method
#' 
#' @details This function extracts information required to perform normalization
#' with global standards. It is useful for running the summarization workflow
#' outside of the dataProcess function.
#' 
#' @export
#' 
#' @examples 
#' input = data.table::as.data.table(DDARawData)
#' peptides_dict = makePeptidesDictionary(input, "GLOBALSTANDARDS")
#' head(peptides_dict) # ready to be passed to the MSstatsNormalize function
#' 
makePeptidesDictionary = function(input, normalization) {
    PEPTIDE = PeptideSequence = PrecursorCharge = NULL
    
    if (toupper(normalization) == "GLOBALSTANDARDS") {
        cols = intersect(c("PeptideSequence", "PeptideModifiedSequence", 
                           "PrecursorCharge"), colnames(input))
        peptides_dict = unique(input[, cols, with = FALSE])
        colnames(peptides_dict)[1] = "PeptideSequence"
        peptides_dict[, PEPTIDE := paste(PeptideSequence, PrecursorCharge, sep = "_")]
    } else {
        NULL
    }
}
