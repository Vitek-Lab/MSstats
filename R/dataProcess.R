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
#' @param min_feature_count optional. Only required if featureSubset = "highQuality".
#' Defines a minimum number of informative features a protein needs to be considered
#' in the feature selection algorithm.
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
#' @param fix_missing Optional, same as the `fix_missing` parameter in MSstatsConvert::MSstatsBalancedDesign function
#' @inheritParams .documentFunction
#' 
#' @importFrom utils sessionInfo
#' @importFrom data.table as.data.table 
#' 
#' @export
#' 
#' @examples 
#' # Consider a raw data (i.e. SRMRawData) for a label-based SRM experiment from a yeast study
#' # with ten time points (T1-T10) of interests and three biological replicates.
#' # It is a time course experiment. The goal is to detect protein abundance changes
#' # across time points.
#' head(SRMRawData)
#' # Log2 transformation and normalization are applied (default)
#' QuantData<-dataProcess(SRMRawData, use_log_file = FALSE)
#' head(QuantData$FeatureLevelData)
#' # Log10 transformation and normalization are applied
#' QuantData1<-dataProcess(SRMRawData, logTrans=10, use_log_file = FALSE)
#' head(QuantData1$FeatureLevelData)
#' # Log2 transformation and no normalization are applied
#' QuantData2<-dataProcess(SRMRawData,normalization=FALSE, use_log_file = FALSE)
#' head(QuantData2$FeatureLevelData)
#' 
dataProcess = function(
    raw, logTrans = 2, normalization = "equalizeMedians", nameStandards = NULL,
    featureSubset = "all", remove_uninformative_feature_outlier = FALSE, 
    min_feature_count = 2, n_top_feature = 3, summaryMethod = "TMP", 
    equalFeatureVar = TRUE, censoredInt = "NA", MBimpute = TRUE, 
    remove50missing = FALSE, fix_missing = NULL, maxQuantileforCensored = 0.999, 
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
    input = MSstatsNormalize(input, normalization, peptides_dict, nameStandards)
    input = MSstatsMergeFractions(input)
    input = MSstatsHandleMissing(input, summaryMethod, MBimpute,
                                 censoredInt, maxQuantileforCensored)
    input = MSstatsSelectFeatures(input, featureSubset, n_top_feature,
                                  min_feature_count)
    .logDatasetInformation(input)
    getOption("MSstatsLog")("INFO",
                            "== Start the summarization per subplot...")
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
    getOption("MSstatsLog")("INFO",
                            "== Summarization is done.")
    getOption("MSstatsMsg")("INFO",
                            " == Summarization is done.")
    output = MSstatsSummarizationOutput(input, summarized, processed,
                                        summaryMethod, MBimpute, censoredInt)
    output
}


#' Feature-level data summarization
#' 
#' @param proteins_list list of processed feature-level data
#' @param method summarization method: "linear" or "TMP" 
#' @param equal_variance only for summaryMethod = "linear". Default is TRUE. 
#' Logical variable for whether the model should account for heterogeneous variation 
#' among intensities from different features. Default is TRUE, which assume equal
#' variance among intensities from features. FALSE means that we cannot assume 
#' equal variance among intensities from features, then we will account for
#' heterogeneous variation from different features.
#' @param censored_symbol Missing values are censored or at random. 
#' 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. 
#' '0' uses zero intensities as censored intensity. 
#' In this case, NA intensities are missing at random. 
#' The output from Skyline should use '0'. 
#' Null assumes that all NA intensites are randomly missing.
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the runs 
#' which have more than 50\% missing values. FALSE is default.
#' @param impute only for summaryMethod = "TMP" and censoredInt = 'NA' or '0'. 
#' TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by Accelated failure model. 
#' FALSE uses the values assigned by cutoffCensored
#' 
#' @importFrom data.table uniqueN
#' @importFrom utils setTxtProgressBar
#' 
#' @return list of length one with run-level data.
#' 
#' @export
#' 
#' @examples
#' raw = DDARawData 
#' method = "TMP"
#' cens = "NA"
#' impute = TRUE
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input = MSstatsPrepareForDataProcess(raw, 2, NULL)
#' input = MSstatsNormalize(input, "EQUALIZEMEDIANS")
#' input = MSstatsMergeFractions(input)
#' input = MSstatsHandleMissing(input, "TMP", TRUE, "NA", 0.999)
#' input = MSstatsSelectFeatures(input, "all")
#' processed = getProcessed(input)
#' input = MSstatsPrepareForSummarization(input, method, impute, cens, FALSE)
#' input_split = split(input, input$PROTEIN)
#' summarized = MSstatsSummarize(input_split, method, impute, cens, FALSE, TRUE)
#' length(summarized) # list of summarization outputs for each protein
#' head(summarized[[1]][[1]]) # run-level summary
#' 
MSstatsSummarize = function(proteins_list, method, impute, censored_symbol,
                            remove50missing, equal_variance) {
    
    num_proteins = length(proteins_list)
    summarized_results = vector("list", num_proteins)
    if (method == "TMP") {
        pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
        for (protein_id in seq_len(num_proteins)) {
            single_protein = proteins_list[[protein_id]]
            summarized_results[[protein_id]] = MSstatsSummarizeSingleTMP(
                single_protein, impute, censored_symbol, remove50missing)
            setTxtProgressBar(pb, protein_id)
        }
        close(pb)
    } else {
        pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
        for (protein_id in seq_len(num_proteins)) {
            single_protein = proteins_list[[protein_id]]
            summarized_result = MSstatsSummarizeSingleLinear(single_protein,
                                                             equal_variance)
            summarized_results[[protein_id]] = summarized_result
            setTxtProgressBar(pb, protein_id)
        }
        close(pb)
    }
    summarized_results
}


#' Linear model-based summarization for a single protein
#' 
#' @param single_protein feature-level data for a single protein
#' @param equal_variances if TRUE, observation are assumed to be homoskedastic
#' 
#' @return list with protein-level data
#' 
#' @importFrom stats xtabs
#' 
#' @export
#' 
#' @examples
#' raw = DDARawData 
#' method = "linear"
#' cens = NULL
#' impute = FALSE 
#' # currently, MSstats only supports MBimpute = FALSE for linear summarization
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input = MSstatsPrepareForDataProcess(raw, 2, NULL)
#' input = MSstatsNormalize(input, "EQUALIZEMEDIANS")
#' input = MSstatsMergeFractions(input)
#' input = MSstatsHandleMissing(input, "TMP", TRUE, "NA", 0.999)
#' input = MSstatsSelectFeatures(input, "all")
#' input = MSstatsPrepareForSummarization(input, method, impute, cens, FALSE)
#' input_split = split(input, input$PROTEIN)
#' single_protein_summary = MSstatsSummarizeSingleLinear(input_split[[1]])
#' head(single_protein_summary[[1]])
#' 
MSstatsSummarizeSingleLinear = function(single_protein, equal_variances = TRUE) {
    ABUNDANCE = RUN = FEATURE = PROTEIN = LogIntensities = NULL
    
    label = data.table::uniqueN(single_protein$LABEL) > 1
    single_protein = single_protein[!is.na(ABUNDANCE)]
    single_protein[, RUN := factor(RUN)]
    single_protein[, FEATURE := factor(FEATURE)]
    
    counts = xtabs(~ RUN + FEATURE, 
                   data = unique(single_protein[, list(FEATURE, RUN)]))
    counts = as.matrix(counts)
    is_single_feature = .checkSingleFeature(single_protein)
    
    fit = try(.fitLinearModel(single_protein, is_single_feature, is_labeled = label, 
                              equal_variances), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
        msg = paste("*** error : can't fit the model for ", unique(single_protein$PROTEIN))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        result = NULL
    } else {
        cf = summary(fit)$coefficients[, 1]
        result = unique(single_protein[, list(Protein = PROTEIN, RUN = RUN)])
        log_intensities = get_linear_summary(single_protein, cf,
                                             counts, label)
        result[, LogIntensities := log_intensities]
    }
    list(result)
}


#' Tukey Median Polish summarization for a single protein
#' 
#' @param single_protein feature-level data for a single protein
#' @inheritParams MSstatsSummarize
#' 
#' @return list of two data.tables: one with fitted survival model,
#' the other with protein-level data
#' 
#' @importFrom stats predict
#' 
#' @export
#' 
#' @examples
#' raw = DDARawData 
#' method = "TMP"
#' cens = "NA"
#' impute = TRUE 
#' # currently, MSstats only supports MBimpute = FALSE for linear summarization
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input = MSstatsPrepareForDataProcess(raw, 2, NULL)
#' input = MSstatsNormalize(input, "EQUALIZEMEDIANS")
#' input = MSstatsMergeFractions(input)
#' input = MSstatsHandleMissing(input, "TMP", TRUE, "NA", 0.999)
#' input = MSstatsSelectFeatures(input, "all")
#' input = MSstatsPrepareForSummarization(input, method, impute, cens, FALSE)
#' input_split = split(input, input$PROTEIN)
#' single_protein_summary = MSstatsSummarizeSingleTMP(input_split[[1]],
#'                                                    impute, cens, FALSE)
#' head(single_protein_summary[[1]])
#' 
MSstatsSummarizeSingleTMP = function(single_protein, impute, censored_symbol, 
                                     remove50missing) {
    newABUNDANCE = n_obs = n_obs_run = RUN = FEATURE = LABEL = NULL
    predicted = censored = NULL
    cols = intersect(colnames(single_protein), c("newABUNDANCE", "cen", "RUN",
                                                 "FEATURE", "ref"))
    single_protein = single_protein[(n_obs > 1 & !is.na(n_obs)) &
                                        (n_obs_run > 0 & !is.na(n_obs_run))]
    if (nrow(single_protein) == 0) {
        return(list(NULL, NULL))
    }
    single_protein[, RUN := factor(RUN)]
    single_protein[, FEATURE := factor(FEATURE)]
    if (impute & any(single_protein[["censored"]])) {
        survival_fit = .fitSurvival(single_protein[LABEL == "L", cols,
                                                   with = FALSE])
        single_protein[, predicted := predict(survival_fit,
                                              newdata = .SD)]
        single_protein[, predicted := ifelse(censored & (LABEL == "L"), predicted, NA)]
        single_protein[, newABUNDANCE := ifelse(censored & LABEL == "L",
                                                predicted, newABUNDANCE)]
        survival = single_protein[, c(cols, "predicted"), with = FALSE]
    } else {
        survival = single_protein[, cols, with = FALSE]
        survival[, predicted := NA]
    }
    
    single_protein = .isSummarizable(single_protein, remove50missing)
    if (is.null(single_protein)) {
        return(list(NULL, NULL))
    } else {
        single_protein = single_protein[!is.na(newABUNDANCE), ]
        is_labeled = nlevels(single_protein$LABEL) > 1
        result = .runTukey(single_protein, is_labeled, censored_symbol, 
                           remove50missing)
    }
    list(result, survival)
}
