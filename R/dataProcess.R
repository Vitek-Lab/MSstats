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
#' @param address the name of folder that will store the results. Default folder 
#' is the current working directory. The other assigned folder has to be existed 
#' under the current working directory. An output csv file is automatically created 
#' with the default name of "BetweenRunInterferenceFile.csv". 
#' The command address can help to specify where to store the file as well as how 
#' to modify the beginning of the file name.
#' @param fillIncompleteRows If the input dataset has incomplete rows, 
#' TRUE (default) adds the rows with intensity value = NA for missing peaks. 
#' FALSE reports error message with list of features which have incomplete rows.
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
#' @param cutoffCensored Cutoff value for censoring with censoredInt = 'NA' or '0'. 
#' Default is 'minFeature', which uses minimum value for each feature.
#' 'minFeatureNRun' uses the smallest between minimum value of corresponding feature 
#' and minimum value of corresponding run. 'minRun' uses minumum value for each run
#' @param MBimpute only for summaryMethod = "TMP" and censoredInt = 'NA' or '0'. 
#' TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) 
#' by Accelated failure model. FALSE uses the values assigned by cutoffCensored.
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the runs 
#' which have more than 50\% missing values. FALSE is default.
#' @param maxQuantileforCensored Maximum quantile for deciding censored missing values, default is 0.999
#' @param clusters a user specified number of clusters for parallel computing. 
#' Default is NULL, which does not use cluster.
#' 
#' @importFrom utils sessionInfo
#' @importFrom data.table as.data.table 
#' 
#' @export
#' 

dataProcess = function(
    raw, logTrans = 2, normalization = "equalizeMedians", nameStandards = NULL,
    address = "", fillIncompleteRows = TRUE, featureSubset = "all", 
    remove_uninformative_feature_outlier = FALSE, n_top_feature = 3, 
    summaryMethod = "TMP", equalFeatureVar = TRUE, censoredInt = "NA", 
    cutoffCensored = "minFeature", MBimpute = TRUE, remove50missing = FALSE,
    maxQuantileforCensored = 0.999, clusters = NULL
) {
    # Setup: check validity of data and parameters ----    
    # Logging is taken care of by MSstatsLogging()
    .saveSessionInfo()
    getOption("MSstatsLog")("INFO", "MSstats - dataProcess function")
    .checkDataProcessParams(
        logTrans, normalization, nameStandards, address, fillIncompleteRows, 
        list(method = featureSubset, n_top = n_top_feature,
             remove_uninformative = remove_uninformative_feature_outlier), 
        list(method = summaryMethod, equal_var = equalFeatureVar),
        list(symbol = censoredInt,
             cutoff = cutoffCensored,
             MB = MBimpute), 
        clusters)
    input = data.table::as.data.table(as(raw, "data.frame"))
    peptides_dict = unique(input[, list(PeptideSequence, PrecursorCharge)])
    peptides_dict[, PEPTIDE := paste(PeptideSequence, PrecursorCharge, sep = "_")]
    input = .checkDataValidity(raw)
    input = .updateColumnsForProcessing(input)
    .preProcessIntensities(input, logTrans)    # rm(raw) # here?
    # Handle fractions and missing run values ----
    check_multi_run = .checkMultiRun(input)
    input = .handleFractions(input, check_multi_run)
    input = .makeBalancedDesign(input, fillIncompleteRows)
    .checkDuplicatedMeasurements(input)
    input = .makeFactorColumns(input)
    # Normalization, Imputation and feature selection ----
    input = .normalize(input, normalization, peptides_dict, nameStandards)
    input = .flagCensored(input, summaryMethod, MBimpute, 
                          censoredInt, maxQuantileforCensored)
    input = .selectFeatures(input, featureSubset, n_top_feature)
    # Record statistics about the dataset
    input = .removeSingleLabelFeatures(input)
    .logDatasetInformation(input)    
    # Summarization per subplot (per RUN) ----
    getOption("MSstatsMsg")("INFO", 
                            "\n == Start the summarization per subplot...")
    summarization = try(.runQuantification(
        input, summaryMethod, equalFeatureVar, cutoffCensored, censoredInt, 
        remove50missing, MBimpute, original_scale = FALSE, logsum = FALSE, 
        featureSubset, remove_uninformative_feature_outlier, 
        message.show = FALSE, clusters = clusters), silent = TRUE)
    .processFinalOutput(input, summarization, summaryMethod)
}
