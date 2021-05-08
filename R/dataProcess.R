#' Process MS data: clean, normalize and summarize before differential analysis
#' 
#' @param raw name of the raw (input) data set.
#' @param logTrans base of logarithm transformation: 2 (default) or 10.
#' @param normalization normalization to remove systematic bias between MS runs. 
#' There are three different normalizations supported:
#' 'equalizeMedians' (default) represents constant normalization (equalizing the medians) 
#' based on reference signals is performed. 
#' 'quantile' represents quantile normalization based on reference signals 
#' 'globalStandards' represents normalization with global standards proteins. 
#' If FALSE, no normalization is performed.
#' @param nameStandards optional vector of global standard peptide names. 
#' Required only for normalization with global standard peptides.
#' @param featureSubset "all" (default) uses all features that the data set has. 
#' "top3" uses top 3 features which have highest average of log-intensity across runs. 
#' "topN" uses top N features which has highest average of log-intensity across runs. 
#' It needs the input for n_top_feature option. 
#' "highQuality" flags uninformative feature and outliers.
#' @param remove_uninformative_feature_outlier optional. Only required if 
#' featureSubset = "highQuality". TRUE allows to remove 
#' 1) noisy features (flagged in the column feature_quality with "Uninformative"),
#' 2) outliers (flagged in the column, is_outlier with TRUE, 
#' before run-level summarization. FALSE (default) uses all features and intensities 
#' for run-level summarization.
#' @param n_top_feature optional. Only required if featureSubset = 'topN'.  
#' It that case, it specifies number of top features that will be used.
#' Default is 3, which means to use top 3 features.
#' @param summaryMethod "TMP" (default) means Tukey's median polish, 
#' which is robust estimation method. "linear" uses linear mixed model.
#' @param equalFeatureVar only for summaryMethod = "linear". default is TRUE. 
#' Logical variable for whether the model should account for heterogeneous variation 
#' among intensities from different features. Default is TRUE, which assume equal 
#' variance among intensities from features. FALSE means that we cannot assume equal 
#' variance among intensities from features, then we will account for heterogeneous 
#' variation from different features.
#' @param censoredInt Missing values are censored or at random. 
#' 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. 
#' '0' uses zero intensities as censored intensity. 
#' In this case, NA intensities are missing at random. 
#' The output from Skyline should use '0'. 
#' Null assumes that all NA intensites are randomly missing.
#' @param MBimpute only for summaryMethod = "TMP" and censoredInt = 'NA' or '0'. 
#' TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) 
#' by Accelated failure model. FALSE uses the values assigned by cutoffCensored.
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the runs 
#' which have more than 50\% missing values. FALSE is default.
#' @param maxQuantileforCensored Maximum quantile for deciding censored missing values, default is 0.999
#' @inheritParams .documentFunction
#' 
#' @importFrom utils sessionInfo
#' @importFrom data.table as.data.table 
#' 
#' @export
#' 

dataProcess = function(
    raw, logTrans = 2, normalization = "equalizeMedians", nameStandards = NULL,
    featureSubset = "all", remove_uninformative_feature_outlier = FALSE, 
    n_top_feature = 3, summaryMethod = "TMP", equalFeatureVar = TRUE, 
    censoredInt = "NA", MBimpute = TRUE, remove50missing = FALSE,
    fix_missing = NULL, maxQuantileforCensored = 0.999, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path,
                                        base = "MSstats_dataProcess_log_")
    getOption("MSstatsLog")("INFO", "MSstats - dataProcess function")
    .checkDataProcessParams(
        logTrans, normalization, nameStandards,
        list(method = featureSubset, n_top = n_top_feature,
             remove_uninformative = remove_uninformative_feature_outlier),
        list(method = summaryMethod, equal_var = equalFeatureVar),
        list(symbol = censoredInt, MB = MBimpute))
    
    peptides_dict = makePeptidesDictionary(as.data.table(unclass(raw)), normalization)
    input = MSstatsPrepareForDataProcess(raw, logTrans, fix_missing)
    # Normalization, Imputation and feature selection ----
    input = MSstatsNormalize(input, normalization, peptides_dict, nameStandards) # MSstatsNormalize
    input = MSstatsMergeFractions(input)
    input = MSstatsHandleMissing(input, summaryMethod, MBimpute,
                                 censoredInt, maxQuantileforCensored)
    input = MSstatsSelectFeatures(input, featureSubset, n_top_feature,
                                  min_feature_count = 2)
    # Record statistics about the dataset
    .logDatasetInformation(input)
    # Summarization per subplot (per RUN) ----
    getOption("MSstatsMsg")("INFO",
                            " == Start the summarization per subplot...")
    processed = getProcessed(input)
    input = MSstatsPrepareForSummarization(input, summaryMethod, MBimpute, censoredInt,
                                           remove_uninformative_feature_outlier)
    input_split = split(input, input$PROTEIN)
    summarized = tryCatch(MSstatsSummarize(input_split, summaryMethod,
                                           MBimpute, censoredInt, 
                                           remove50missing, equalFeatureVar),
                          error = function(e) {
                              print(e)
                              NULL
                          })
    output = MSstatsSummarizationOutput(input, summarized, processed,
                                        summaryMethod, MBimpute, censoredInt)
    output
}
