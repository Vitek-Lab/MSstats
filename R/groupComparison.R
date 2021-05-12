#' Whole plot testing
#'
#' @param contrast.matrix comparison between conditions of interests.
#' @param data name of the (output of dataProcess function) data set.
#' @param save_fitted_models logical, if TRUE, fitted models will be added to
#' the output.
#' @param log_base base of the logarithm used in dataProcess.
#' @inheritParams .documentFunction
#'
#' @details
#' contrast.matrix : comparison of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically. Command levels(QuantData$FeatureLevelData$GROUP_ORIGINAL) can illustrate the actual order of the levels of conditions.
#' The underlying model fitting functions are lm and lmer for the fixed effects model and mixed effects model, respectively.
#' The input of this function is the quantitative data from function (dataProcess).
#'
#' @return list that consists of the following elements: (TODO)
#' @export 
#' @import lme4
#' @import limma
#' @importFrom data.table rbindlist
#'
groupComparison = function(contrast.matrix, data, 
                           save_fitted_models = TRUE, log_base = 2,
                           use_log_file = TRUE, append = FALSE, 
                           verbose = TRUE, log_file_path = NULL
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        "MSstats_groupComparison_log_")
    getOption("MSstatsLog")("INFO", "MSstats - groupComparison function")
    labeled = data.table::uniqueN(data$FeatureLevelData$Label) > 1
    split_summarized = MSstatsPrepareForGroupComparison(data)
    repeated = checkRepeatedDesign(data)
    samples_info = getSamplesInfo(data)
    groups = unique(data$ProteinLevelData$GROUP)
    contrast_matrix = MSstatsContrastMatrix(contrast.matrix, groups)
    getOption("MSstatsLog")("INFO",
                            "== Start to test and get inference in whole plot")
    getOption("MSstatsMsg")("INFO",
                            " == Start to test and get inference in whole plot ...")
    testing_results = MSstatsGroupComparison(split_summarized, contrast_matrix,
                                             save_fitted_models, repeated, samples_info)
    getOption("MSstatsLog")("INFO",
                            "== Comparisons for all proteins are done.")
    getOption("MSstatsMsg")("INFO",
                            " == Comparisons for all proteins are done.")
    MSstatsGroupComparisonOutput(testing_results, data, log_base)
}


#' Prepare output for dataProcess for group comparison
#' 
#' @param summarization_output output of dataProcess
#' 
#' @export
#' 
MSstatsPrepareForGroupComparison = function(summarization_output) {
    has_imputed = is.element("NumImputedFeature", colnames(summarization_output$ProteinLevelData))
    summarized = data.table::as.data.table(summarization_output$ProteinLevelData)
    summarized = .checkGroupComparisonInput(summarized)
    labeled = nlevels(summarization_output$FeatureLevelData$LABEL) > 1
    
    getOption("MSstatsLog")("INFO", paste0("labeled = ", labeled))
    getOption("MSstatsLog")("INFO", "scopeOfBioReplication = expanded")
    output = split(summarized, summarized$Protein)
    attr(output, "has_imputed") = has_imputed
    output
}


#' Group comparison
#' 
#' @param summarized_list output of MSstatsPrepareForGroupComparison
#' @param contrast_matrix contrast matrix
#' @param save_fitted_models if TRUE, fitted models will be included in the output
#' @param repeated logical, output of checkRepeatedDesign function
#' @param samples_info data.table, output of getSamplesInfo function
#' 
#' @export
#' 
MSstatsGroupComparison = function(summarized_list, contrast_matrix,
                                  save_fitted_models, repeated, samples_info) {
    groups = sort(colnames(contrast_matrix))
    has_imputed = attr(summarized_list, "has_imputed")
    all_proteins_id = seq_along(summarized_list)
    test_results = vector("list", length(all_proteins_id))
    pb = txtProgressBar(max = length(all_proteins_id), style = 3)
    for (i in all_proteins_id) {
        comparison_outputs = MSstatsGroupComparisonSingleProtein(
            summarized_list[[i]], contrast_matrix, repeated, 
            groups, samples_info, save_fitted_models, has_imputed
        )
        test_results[[i]] = comparison_outputs
        setTxtProgressBar(pb, i)
    }
    close(pb)
    test_results
}


#' Create output of group comparison based on results for individual proteins
#' 
#' @param input output of MSstatsGroupComparison function
#' @param summarization_output output of dataProcess function
#' @param log_base base of the logarithm used in fold-change calculation
#' 
#' @export
#' 
MSstatsGroupComparisonOutput = function(input, summarization_output, log_base = 2) {
    has_imputed = is.element("NumImputedFeature", colnames(summarization_output$ProteinLevelData))
    model_qc_data = lapply(input, function(x) x[[1]])
    comparisons = lapply(input, function(x) x[[2]])
    fitted_models = lapply(input, function(x) x[[3]])
    comparisons = data.table::rbindlist(comparisons, fill = TRUE)
    comparisons[, adj.pvalue := p.adjust(pvalue, method = "BH"),
                by = "Label"]
    logFC_colname = paste0("log", log_base, "FC")
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
    getOption("MSstatsLog")("INFO", "The output for groupComparison is ready.")
    list(ComparisonResult = as.data.frame(comparisons)[, cols],
         ModelQC = as.data.frame(qc),
         FittedModel = fitted_models)   
}


#' Group comparison for a single protein
#' 
#' @param single_protein data.table with summarized data for a single protein
#' @param contrast_matrix contrast matrix
#' @param repeated if TRUE, repeated measurements will be modeled
#' @param groups unique labels of experimental conditions
#' @param samples_info number of runs per group
#' @param save_fitted_models if TRUE, fitted model will be saved.
#' If not, it will be replaced with NULL
#' @param has_imputed TRUE if missing values have been imputed
#' 
#' @export
#' 
MSstatsGroupComparisonSingleProtein = function(single_protein, contrast_matrix,
                                               repeated, groups, samples_info,
                                               save_fitted_models,
                                               has_imputed) {
    single_protein = .prepareSingleProteinForGC(single_protein)
    is_single_subject = .checkSingleSubject(single_protein)
    has_tech_reps = .checkTechReplicate(single_protein)
    
    fitted_model = try(.fitModelSingleProtein(single_protein, contrast_matrix,
                                              has_tech_reps, is_single_subject,
                                              repeated, groups, samples_info,
                                              save_fitted_models, has_imputed),
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
