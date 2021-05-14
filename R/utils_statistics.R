#' Log information about feature-level data
#' @param input data.table
#' @return TRUE invisibly after successful logging
#' @keywords internal
.logDatasetInformation = function(input) {
    .logSummaryStatistics(input)
    .checkSingleLabelProteins(input)
    .logMissingness(input)
    invisible(TRUE)
}


#' Print summary statistics to the log file
#' @param input data.table
#' @return TRUE invisibly
#' @keywords internal
.logSummaryStatistics = function(input) {
    RUN = SUBJECT_ORIGINAL = FRACTION = GROUP_ORIGINAL = NumRuns = NULL
    NumBioReplicates = NumFractions = NULL
    
    PEPTIDE = FEATURE = PROTEIN = feature_count = NULL
    
    num_proteins = data.table::uniqueN(input$PROTEIN)
    peptides_per_protein = input[, list(peptide_count = data.table::uniqueN(PEPTIDE)),
                                 by = "PROTEIN"]
    features_per_peptide = input[, list(feature_count = data.table::uniqueN(FEATURE)),
                                 by = "PEPTIDE"]
    features_per_protein = input[, list(feature_count = data.table::uniqueN(FEATURE)),
                                 by = "PROTEIN"]
    features_per_protein = features_per_protein[feature_count == 1L, PROTEIN]
    
    pep_per_prot = range(peptides_per_protein$peptide_count)
    feat_per_pept = range(features_per_peptide$feature_count)
    counts_msg = paste("", paste("# proteins:", num_proteins),
                       paste("# peptides per protein:", 
                             paste(pep_per_prot, sep = "-", collapse = "-")),
                       paste("# features per peptide:",
                             paste(feat_per_pept, sep = "-", collapse = "-")),
                       sep = "\n ", collapse = "\n ")
    getOption("MSstatsMsg")("INFO", counts_msg)
    getOption("MSstatsLog")("INFO", counts_msg)
    
    if (length(features_per_protein) > 0) {
        single_features = unique(as.character(features_per_protein))
        n_feat = min(length(single_features), 5)
        msg = paste("Some proteins have only one feature:", "\n",
                    paste(unique(as.character(features_per_protein))[1:n_feat],
                          sep = ",\n ", collapse = ",\n "),
                    "...")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    samples_info = input[, list(NumRuns = data.table::uniqueN(RUN),
                                NumBioReplicates = data.table::uniqueN(SUBJECT_ORIGINAL),
                                NumFractions = data.table::uniqueN(FRACTION)), 
                         by = "GROUP_ORIGINAL"]
    samples_info = samples_info[
        , 
        list(GROUP_ORIGINAL, NumRuns, NumBioReplicates, 
             NumTechReplicates = as.integer(
                 round(NumRuns / (NumBioReplicates * NumFractions)
                 )))]
    samples_info = data.table::dcast(data.table::melt(samples_info, 
                                                      id.vars = "GROUP_ORIGINAL"), 
                                     variable ~ GROUP_ORIGINAL)
    colnames(samples_info)[1] = ""
    samples_info[, 1] = c("# runs", "# bioreplicates", "# tech. replicates")
    
    samples_info = rbind(data.table::data.table(t(colnames(samples_info))), 
                         samples_info, use.names = FALSE)
    samples_msg = apply(samples_info, 2, .nicePrint)
    samples_msg = apply(samples_msg, 1, function(x) paste0(x, collapse = ""))
    samples_msg = paste("", samples_msg, sep = "\n")
    getOption("MSstatsLog")("INFO", samples_msg)
    getOption("MSstatsMsg")("INFO", samples_msg)    
    invisible(TRUE)
}


#' Print a table nicely 
#' @param string_vector character
#' @return character
#' @keywords internal
.nicePrint = function(string_vector) {
    max_chars = max(nchar(string_vector))   
    string_vector = sapply(string_vector, 
                           function(x) paste(paste0(rep(" ", max_chars - nchar(x)), collapse = ""), x),
                           USE.NAMES = FALSE)
} 


#' Print proteins with a single label to the log file
#' @param input data.table
#' @param label label ("L" or "H")
#' @return TRUE invisibly
#' @keywords internal
.logSingleLabeledProteins = function(input, label) {
    LABEL = PROTEIN = NULL
    
    name = ifelse(label == "L", "endogeneous", "reference")
    proteins = unique(input[LABEL == label, as.character(PROTEIN)])
    if (length(proteins) > 0) {
        n_prot = min(length(proteins), 5)
        msg = paste(paste("Some proteins only have", name,
                          "intensities in label-based experiment",
                          "Please check or remove these proteins:"),
                    paste(proteins[1:n_prot], sep = ", \n ", collapse = ", \n "),
                    "... (see the log file for a full list)")
        getOption("MSstatsMsg")("WARN", msg)
        getOption("MsstatsLog")("WARN", msg)
    }
    invisible(TRUE)
}


#' Check if there are proteins with a single label in a labeled dataset
#' @param input data.table
#' @return TRUE invisibly
#' @keywords internal
.checkSingleLabelProteins = function(input) {
    LABEL = label_count = NULL
    
    if (data.table::uniqueN(input$LABEL) == 2) {
        labels_by_protein = unique(input[, 
                                         list(LABEL, 
                                              label_count = data.table::uniqueN(LABEL)),
                                         by = "PROTEIN"])
        labels_by_protein = labels_by_protein[label_count == 1L, ]
        
        .logSingleLabeledProteins(labels_by_protein, "L")
        .logSingleLabeledProteins(labels_by_protein, "H")
    }
    invisible(TRUE)
}


#' Log information about missing data
#' @param input data.table
#' @return TRUE invisibly
#' @keywords internal
.logMissingness = function(input) {
    ABUNDANCE = AllMissing = NumMissing = NumTotal = AnyAllMissing = NULL
    FEATURE = FractionMissing = RUN = NULL
    
    missing = input[, list(NumMissing = sum(is.na(ABUNDANCE), na.rm = TRUE),
                           NumTotal = .N), 
                    by = c("LABEL", "GROUP", "FEATURE")]
    missing[, AllMissing := NumMissing == NumTotal]
    missing[, AnyAllMissing := any(AllMissing), by = c("LABEL", "FEATURE")]
    missing_in_any = as.character(missing[(AnyAllMissing), FEATURE])
    # TODO: LOG, which features have AnyAllMissing
    
    missing_by_run = input[, list(NumMissing = sum(is.na(ABUNDANCE), na.rm = TRUE),
                                  NumTotal = .N), by = "RUN"]
    missing_by_run[, FractionMissing := NumMissing / NumTotal]
    missing_by_run = as.character(missing_by_run[FractionMissing > 0.75, 
                                                 unique(RUN)])
    
    if (length(missing_in_any) > 0) {
        num_features = min(length(missing_in_any), 5)
        msg = paste("Some features are completely",
                    "missing in at least one condition: ", "\n",
                    paste(unique(as.character(missing_in_any))[1:num_features], 
                          sep = ",\n ", 
                          collapse = ",\n "), "...")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    if (length(missing_by_run) > 0) {
        msg = paste("The following runs have more than 75% missing values:",
                    paste(missing_by_run, sep = ",\n ", collapse = ",\n "))
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    invisible(TRUE)
}
