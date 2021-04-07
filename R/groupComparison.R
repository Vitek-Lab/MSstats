#' Whole plot testing
#'
#' @param contrast.matrix comparison between conditions of interests.
#' @param data name of the (output of dataProcess function) data set.
#' @param save_fitted_models logical, if TRUE, fitted models will be added to
#' the output.
#'
#' @details
#' contrast.matrix : comparison of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically. Command levels(QuantData$ProcessedData$GROUP_ORIGINAL) can illustrate the actual order of the levels of conditions.
#' The underlying model fitting functions are lm and lmer for the fixed effects model and mixed effects model, respectively.
#' The input of this function is the quantitative data from function (dataProcess).
#'
#' @return list that consists of the following elements: (TODO)
#' @export groupComparison
#' @import lme4
#' @import limma
#' @importFrom data.table rbindlist
#'

groupComparison <- function(contrast.matrix, data, save_fitted_models = TRUE) {
    Protein = issue = pvalue = PROTEIN = NULL
    # Save session information here
    # Initialize groupComparison log here
    getOption("MSstatsLog")("INFO", "MSstats - groupComparison function")
    
    data = .checkGroupComparisonInput(data)
    summarized = data.table::as.data.table(data$RunlevelData)
    processed = data.table::as.data.table(data$ProcessedData)
    contrast_matrix = .checkContrastMatrix(contrast.matrix, processed)
    scopeOfBioReplication = "expanded"
    labeled = nlevels(data$ProcessedData$LABEL) > 1
    repeated = .checkRepeated(processed)
    has_imputed = is.element("NumImputedFeature", colnames(data$RunlevelData))
    
    getOption("MSstatsLog")("INFO", paste0("labeled = ", labeled))
    getOption("MSstatsLog")("INFO", paste0("scopeOfBioReplication = ",
                                           scopeOfBioReplication))
    getOption("MSstatsLog")("INFO",
                            "** Start to test and get inference in whole plot")
    
    groups = unique(summarized$GROUP)
    all_proteins = unique(summarized$Protein)
    group_comparison = vector("list", length(all_proteins))
    model_qc_data = vector("list", length(all_proteins))
    fitted_models = vector("list", length(all_proteins))
    pb = txtProgressBar(max = nlevels(all_proteins), style = 3)
    for (i in 1:nlevels(all_proteins)) {
        single_protein = summarized[Protein == all_proteins[i]]
        processed_single = processed[PROTEIN == all_proteins[i]]
        comparison_outputs = .groupComparisonSingleProtein(
            single_protein, contrast_matrix,
            repeated, groups, save_fitted_models,
            processed_single, has_imputed
        )
        model_qc_data[[i]] = comparison_outputs[[1]]
        group_comparison[[i]] = comparison_outputs[[2]]
        fitted_models[[i]] = comparison_outputs[[3]]
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    getOption("MSstatsLog")("INFO",
                            "Comparisons for all proteins are done.- okay")
    comparisons = data.table::rbindlist(group_comparison, fill = TRUE)
    comparisons[, adj.pvalue := p.adjust(pvalue, method = "BH"),
                by = "Label"]
    logFC_colname = .getLogBaseName(processed)
    comparisons[, adj.pvalue := ifelse(!is.na(issue) & 
                                           issue == "oneConditionMissing", 
                                       0, adj.pvalue)]
    data.table::setnames(comparisons, "logFC", logFC_colname)
    qc = rbindlist(model_qc_data, fill = TRUE)
    cols = c("Protein", "Label", logFC_colname, "SE", "Tvalue", "DF", 
             "pvalue", "adj.pvalue", "issue", "MissingPercentage",
             "ImputationPercentage")
    if (!has_imputed) {
        cols = cols[1:10]
    }
    getOption("MSstatsLog")("INFO", "Group comparison is done. - okay")
    list(ComparisonResult = as.data.frame(comparisons)[, cols],
         ModelQC = as.data.frame(qc),
         fittedmodel = fitted_models)   
}


#' Group comparison for a single protein
#' @param single_protein data.table with summarized data for a single protein
#' @param contrast_matrix contrast matrix
#' @param repeated if TRUE, repeated measurements will be modeled
#' @param groups unique labels of experimental conditions
#' @param save_fitted_models if TRUE, fitted model will be saved.
#' If not, it will be replaced with NULL
#' @param processed data.table with processed data for a single protein
#' @param has_imputed TRUE if missing values have been imputed
#' @keywords internal
.groupComparisonSingleProtein = function(single_protein, contrast_matrix,
                                         repeated, groups, save_fitted_models,
                                         processed, has_imputed) {
    single_protein = .prepareSingleProteinForGC(single_protein)
    is_single_subject = .checkSingleSubject(single_protein)
    has_tech_reps = .checkTechReplicate(single_protein)
    
    fitted_model = try(.fitModelSingleProtein(single_protein, contrast_matrix,
                                              has_tech_reps, is_single_subject,
                                              repeated, groups, 
                                              save_fitted_models, processed,
                                              has_imputed),
                       silent = TRUE)
    if (inherits(fitted_model, "try-error")) {
        result = list(list(Protein = unique(single_protein$Protein),
                           Label = row.names(contrast_matrix),
                           logFC = NA, SE = NA, Tvalue = NA,
                           DF = NA, pvalue = NA, issue = NA), NULL)
    } else {
        result = fitted_model
    }
    list(single_protein, result[[1]], result[[2]])
}
