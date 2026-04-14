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
#' If FALSE, no normalization is performed.  See MSstats vignettes for 
#' recommendations on which normalization option to use.
#' @param nameStandards optional vector of global standard peptide names. 
#' Required only for normalization with global standard peptides.
#' @param featureSubset "topN" (default) uses top N features which has highest average of log-intensity across runs. 
#' "top3" uses top 3 features which have highest average of log-intensity across runs. 
#' "all" uses all features that the data set has (not recommended in DIA experiments).
#' It needs the input for n_top_feature option. 
#' "highQuality" flags uninformative feature and outliers. See MSstats vignettes for 
#' recommendations on which feature selection option to use.
#' @param remove_uninformative_feature_outlier optional. Only required if 
#' featureSubset = "highQuality". TRUE allows to remove 
#' 1) noisy features (flagged in the column feature_quality with "Uninformative"),
#' 2) outliers (flagged in the column, is_outlier with TRUE, 
#' before run-level summarization. FALSE (default) uses all features and intensities 
#' for run-level summarization.
#' @param min_feature_count optional. Only required if featureSubset = "highQuality".
#' Defines a minimum number of informative features a protein needs to be considered
#' in the feature selection algorithm.
#' @param n_top_feature Specifies the number of top features to use in summarization (100 default). 
#' Only required if featureSubset = 'topN'.  
#' Default is 100, which means to use top 100 features. 
#' Smaller numbers can be set to improve processing times. This option is by default on 
#' at a high number (100) to improve processing times without affecting differential analysis.
#' @param summaryMethod "TMP" (default) means Tukey's median polish, 
#' which is robust estimation method. "linear" uses linear mixed model. If 
#' anomaly detection algorithm is performed, "linear" must be used.
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
#' TRUE (default) imputes missing values with 'NA' or '0' (depending on censoredInt option) 
#' by Accelerated failure model. If set to FALSE, no missing values are imputed. 
#' FALSE is appropriate only when missingness is assumed to be at random.
#' See MSstats vignettes for recommendations on which imputation option to use.
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the proteins 
#' where every run has at least 50\% missing values for each peptide. FALSE is default.
#' @param maxQuantileforCensored Maximum quantile for deciding censored missing values, default is 0.999
#' @param fix_missing Optional, same as the `fix_missing` parameter in MSstatsConvert::MSstatsBalancedDesign function
#' @param numberOfCores Number of cores for parallel processing. When > 1, 
#' a logfile named `MSstats_dataProcess_log_progress.log` is created to 
#' track progress. Only works for Linux & Mac OS. Default is 1.
#' @param aft_iterations Number of iterations for AFT model fitting. Default is 90.
#' @inheritParams .documentFunction
#' 
#' @importFrom utils sessionInfo
#' @importFrom data.table as.data.table
#' 
#' @return A list containing:
#' \describe{
#'   \item{FeatureLevelData}{A data frame with feature-level information after processing. Columns include:
#'     \describe{
#'       \item{PROTEIN}{Identifier for the protein associated with the feature.}
#'       \item{PEPTIDE}{Identifier for the peptide sequence.}
#'       \item{TRANSITION}{Identifier for the transition, typically representing a specific ion pair.}
#'       \item{FEATURE}{Unique identifier for the feature, which could be a combination of peptide and transition.}
#'       \item{LABEL}{Specifies the isotopic labeling of peptides, notably for SRM-based experiments. "L" indicates light-labeled peptides while "H" denotes heavy-labeled peptides.}
#'       \item{GROUP}{Experimental group identifier.}
#'       \item{RUN}{Identifier for the specific MS run.}
#'       \item{SUBJECT}{Subject identifier within the experimental group.}
#'       \item{FRACTION}{Fraction identifier if fractionation was performed.}
#'       \item{originalRUN}{Original run identifier before any processing.}
#'       \item{censored}{Logical indicator of whether the intensity value is considered missing or below limit of detection.}
#'       \item{INTENSITY}{Original intensity measurement of the feature in the given run.}
#'       \item{ABUNDANCE}{Processed abundance or intensity value after log-transformation and normalization.}
#'       \item{newABUNDANCE}{The ABUNDANCE column but includes imputed missing values. It is the column that is used for protein summarization.}
#'       \item{predicted}{Predicted intensity values for censored data, typically derived from a statistical model.}
#'     }
#'   }
#'   \item{ProteinLevelData}{A data frame with run-level summarized information for each protein. Columns include:
#'     \describe{
#'       \item{RUN}{Identifier for the specific MS run.}
#'       \item{Protein}{Identifier for the protein.}
#'       \item{LogIntensities}{Log-transformed intensities for the protein in each run.}
#'       \item{originalRUN}{Original run identifier before any processing.}
#'       \item{GROUP}{Experimental group identifier.}
#'       \item{SUBJECT}{Subject identifier within the experimental group.}
#'       \item{TotalGroupMeasurements}{Total number of feature measurements for the protein in the given group.}
#'       \item{NumMeasuredFeatures}{Number of features measured for the protein in the given run.}
#'       \item{MissingPercentage}{Percentage of missing feature values for the protein in the given run.}
#'       \item{more50missing}{Logical indicator of whether more than 50 percent of the features values are missing for the protein in the given run.}
#'       \item{NumImputedFeature}{Number of features for which values were imputed due to missing or censored data for the protein in the given run.}
#'     }
#'   }
#' }
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
    featureSubset = "topN", remove_uninformative_feature_outlier = FALSE, 
    min_feature_count = 2, n_top_feature = 100, summaryMethod = "TMP", 
    equalFeatureVar = TRUE, censoredInt = "NA", MBimpute = TRUE, 
    remove50missing = FALSE, fix_missing = NULL, maxQuantileforCensored = 0.999, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    numberOfCores = 1, aft_iterations=90
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
        list(symbol = censoredInt, MB = MBimpute),
        colnames(raw))
    
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
    summarized = tryCatch(MSstatsSummarizeWithMultipleCores(input, summaryMethod,
                                           MBimpute, censoredInt, 
                                           remove50missing, equalFeatureVar, 
                                           numberOfCores, aft_iterations),
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

#' Feature-level data summarization with multiple cores
#' 
#' @param input feature-level data processed by dataProcess subfunctions
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
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the proteins 
#' where every run has at least 50\% missing values for each peptide. FALSE is default.
#' @param impute only for summaryMethod = "TMP" and censoredInt = 'NA' or '0'. 
#' TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by Accelated failure model. 
#' FALSE uses the values assigned by cutoffCensored
#' @param numberOfCores Number of cores for parallel processing. When > 1, 
#' a logfile named `MSstats_dataProcess_log_progress.log` is created to 
#' track progress. Only works for Linux & Mac OS. Default is 1.
#' @param aft_iterations Number of iterations for AFT model fitting. Default is 90.
#' 
#' @importFrom parallel makeCluster parLapply stopCluster clusterExport
#' 
#' @return list of length one with run-level data.
#' 
MSstatsSummarizeWithMultipleCores = function(input, method, impute, censored_symbol,
                              remove50missing, equal_variance, numberOfCores = 1,
                              aft_iterations = 90) {
    if (numberOfCores > 1) {
        is_labeled_reference = "is_labeled_ref" %in% colnames(input) && any(input$is_labeled_ref, na.rm = TRUE)
        if (is_labeled_reference) {
            protein_indices = split(seq_len(nrow(input)), list(input$PROTEIN))
        } else {
            protein_indices = split(seq_len(nrow(input)), list(input$PROTEIN, input$LABEL))
        }
        num_proteins = length(protein_indices)
        function_environment = environment()
        cl = parallel::makeCluster(numberOfCores)
        getOption("MSstatsLog")("INFO",
                                "Starting the cluster setup for summarization")
        parallel::clusterExport(cl, c("MSstatsSummarizeSingleTMP", 
                                      "MSstatsSummarizeSingleLinear",
                                      "input", "impute", "censored_symbol",
                                      "remove50missing", "protein_indices", 
                                      "equal_variance", "aft_iterations"), 
                                envir = function_environment)
        cat(paste0("Number of proteins to process: ", num_proteins), 
            sep = "\n", file = "MSstats_dataProcess_log_progress.log")
        if (method == "TMP") {
            summarized_results = parallel::parLapply(cl, seq_len(num_proteins), function(i) {
                if (i %% 100 == 0) {
                    cat("Finished processing an additional 100 proteins", 
                        sep = "\n", file = "MSstats_dataProcess_log_progress.log", append = TRUE)
                }
                single_protein = input[protein_indices[[i]],]
                MSstatsSummarizeSingleTMP(
                    single_protein, impute, censored_symbol, remove50missing,
                    aft_iterations)
            })
        } else {
            summarized_results = parallel::parLapply(cl, seq_len(num_proteins), function(i) {
                if (i %% 100 == 0) {
                    cat("Finished processing an additional 100 proteins", 
                        sep = "\n", file = "MSstats_dataProcess_log_progress.log", append = TRUE)
                }
                single_protein = input[protein_indices[[i]],]
                MSstatsSummarizeSingleLinear(
                    single_protein,
                    impute, 
                    censored_symbol, 
                    remove50missing,
                    aft_iterations)
            })
        }
        parallel::stopCluster(cl)
        return(summarized_results)
    } else {
        return(MSstatsSummarizeWithSingleCore(input, method, impute, 
                                              censored_symbol, 
                                              remove50missing, 
                                              equal_variance,
                                              aft_iterations))
    }
}

#' Feature-level data summarization with 1 core
#' 
#' @inheritParams MSstatsSummarizeWithMultipleCores
#' @param aft_iterations Number of iterations for AFT model fitting. Default is 90.
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
#' summarized = MSstatsSummarizeWithSingleCore(input, method, impute, cens, FALSE, TRUE, 100)
#' length(summarized) # list of summarization outputs for each protein
#' head(summarized[[1]][[1]]) # run-level summary
#' 
MSstatsSummarizeWithSingleCore = function(input, method, impute, censored_symbol,
                            remove50missing, equal_variance, aft_iterations = 90) {


    is_labeled_reference = "is_labeled_ref" %in% colnames(input) && any(input$is_labeled_ref, na.rm = TRUE)
    if (is_labeled_reference) {
        protein_indices = split(seq_len(nrow(input)), list(input$PROTEIN))
    } else {
        protein_indices = split(seq_len(nrow(input)), list(input$PROTEIN, input$LABEL))
    }
    num_proteins = length(protein_indices)
    summarized_results = vector("list", num_proteins)
    if (method == "TMP") {
        pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
        for (protein_id in seq_len(num_proteins)) {
            single_protein = input[protein_indices[[protein_id]],]
            summarized_results[[protein_id]] = MSstatsSummarizeSingleTMP(
                single_protein, impute, censored_symbol, remove50missing, 
                aft_iterations)
            setTxtProgressBar(pb, protein_id)
        }
        close(pb)
    } else {
        pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
        for (protein_id in seq_len(num_proteins)) {
            single_protein = input[protein_indices[[protein_id]],]
            summarized_result = MSstatsSummarizeSingleLinear(
                single_protein, impute, censored_symbol, 
              remove50missing, aft_iterations)

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
#' @param impute boolean for whether imputation should be performed
#' @param censored_symbol Character string indicating how censored values are represented 
#' @param remove50missing if TRUE, proteins with more than 50\% missing values in each run are removed
#' @param aft_iterations number of iterations for AFT model fitting
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
#' single_protein_summary = MSstatsSummarizeSingleLinear(input_split[[1]], impute, cens, TRUE, 100)
#' head(single_protein_summary[[1]])
#' 
#' 
MSstatsSummarizeSingleLinear = function(single_protein,
                                        impute,
                                        censored_symbol,
                                        remove50missing,
                                        aft_iterations = 90,
                                        equal_variances = TRUE) {
    ABUNDANCE = RUN = FEATURE = PROTEIN = LogIntensities = NULL

    cols = intersect(
      colnames(single_protein),
      c("newABUNDANCE", "cen", "RUN", "FEATURE", "ref_covariate")
    )

    is_labeled_reference = "is_labeled_ref" %in% colnames(single_protein) &&
        any(single_protein$is_labeled_ref, na.rm = TRUE)

    single_protein = single_protein[
      (n_obs > 1 & !is.na(n_obs)) &
        (n_obs_run > 0 & !is.na(n_obs_run))
    ]

    if (nrow(single_protein) == 0) {
        return(list(NULL, NULL))
    }

    single_protein[, RUN := factor(RUN)]
    single_protein[, FEATURE := factor(FEATURE)]

    if (impute & any(single_protein[["censored"]])) {
        fit_data = if (is_labeled_reference) {
            single_protein[is_labeled_ref == FALSE, cols, with = FALSE]
        } else {
            single_protein[, cols, with = FALSE]
        }
        survival_fit = .fitSurvival(fit_data, aft_iterations)
        sigma2 = survival_fit$scale^2

        single_protein[, c("predicted", "imputation_var") := {
            pred = predict(survival_fit, newdata = .SD, se.fit = TRUE)
            list(pred$fit, pred$se.fit^2 + sigma2)
        }]

        if (is_labeled_reference) {
            single_protein[, predicted := ifelse(censored & is_labeled_ref == FALSE, predicted, NA)]
            single_protein[, newABUNDANCE := ifelse(censored & is_labeled_ref == FALSE, predicted, newABUNDANCE)]
        } else {
            single_protein[, predicted := ifelse(censored, predicted, NA)]
            single_protein[, newABUNDANCE := ifelse(censored, predicted, newABUNDANCE)]
        }

        survival = single_protein[, intersect(c(cols, "LABEL", "predicted"), colnames(single_protein)), with = FALSE]
    } else {
        survival = single_protein[, intersect(c(cols, "LABEL"), colnames(single_protein)), with = FALSE]
        survival[, predicted := NA]
    }
    
    if (all(!is.na(single_protein$ANOMALYSCORES))) {
        single_protein[, weights :=
            anomaly_weights_z_vec(ANOMALYSCORES),
          by = .(PROTEIN, PEPTIDE)]
    } else {
        single_protein[, weights := NA_real_]
    }
    
    label = data.table::uniqueN(single_protein$LABEL) > 1
    
    single_protein = single_protein[!is.na(newABUNDANCE)]
    
    if (nrow(single_protein) == 0) {
        return(list(NULL, survival))
    }
    
    single_protein[, RUN := factor(RUN)]
    single_protein[, FEATURE := factor(FEATURE)]
    
    is_single_feature = .checkSingleFeature(single_protein)
    
    if (is_single_feature) {
        result = single_protein[, .(LogIntensities = mean(newABUNDANCE)), by = RUN]
        result[, Protein := unique(single_protein$PROTEIN)]
        result[, LABEL := unique(single_protein$LABEL)]
        result[, Variance := NA_real_]
        setcolorder(result, c("Protein", "RUN", "LogIntensities", "Variance"))

        return(list(result, survival))
    } else {
        counts = xtabs(
          ~ RUN + FEATURE,
          data = unique(single_protein[, .(FEATURE, RUN)])
        )
        counts = as.matrix(counts)

        fit = try(
          .fitLinearModel(single_protein, is_single_feature,
                          is_labeled = label,
                          equal_variances),
          silent = TRUE
        )

        if (inherits(fit, "try-error")) {
            msg = paste("*** error : can't fit the model for",
                        unique(single_protein$PROTEIN))
            getOption("MSstatsLog")("WARN", msg)
            getOption("MSstatsMsg")("WARN", msg)
            result = NULL
        } else {
            cf = summary(fit)$coefficients[, 1]
            cov_mat = vcov(fit)

            result = unique(single_protein[, .(Protein = PROTEIN, RUN = RUN)])
            extracted_values = get_linear_summary(single_protein, cf,
                                                  counts, label, cov_mat)
            result = cbind(result, extracted_values)
            result[, LABEL := unique(single_protein$LABEL)]
        }
        
        return(list(result, survival))
    }
}


#' Tukey Median Polish summarization for a single protein
#' 
#' @param single_protein feature-level data for a single protein
#' @param aft_iterations number of iterations for AFT model fitting
#' @inheritParams MSstatsSummarizeWithSingleCore
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
#'                                                    impute, cens, FALSE, 100)
#' head(single_protein_summary[[1]])
#' 
MSstatsSummarizeSingleTMP = function(single_protein, impute, censored_symbol,
                                     remove50missing, aft_iterations = 90) {
    newABUNDANCE = n_obs = n_obs_run = RUN = FEATURE = LABEL = NULL
    predicted = censored = NULL
    cols = intersect(colnames(single_protein), c("newABUNDANCE", "cen", "RUN",
                                                 "FEATURE", "ref_covariate"))
    is_labeled_reference = "is_labeled_ref" %in% colnames(single_protein) &&
        any(single_protein$is_labeled_ref, na.rm = TRUE)
    single_protein = single_protein[(n_obs > 1 & !is.na(n_obs)) &
                                        (n_obs_run > 0 & !is.na(n_obs_run))]
    if (nrow(single_protein) == 0) {
        return(list(NULL, NULL))
    }
    single_protein[, RUN := factor(RUN)]
    single_protein[, FEATURE := factor(FEATURE)]
    if (impute & any(single_protein[["censored"]])) {

        # Flag to track convergence warning
        converged = TRUE

        fit_data = if (is_labeled_reference) {
            single_protein[is_labeled_ref == FALSE, cols, with = FALSE]
        } else {
            single_protein[, cols, with = FALSE]
        }

        # Try to fit survival model and catch convergence warnings
        survival_fit = withCallingHandlers({
            .fitSurvival(fit_data, aft_iterations)
        }, warning = function(w) {
            if (grepl("converge", conditionMessage(w), ignore.case = TRUE)) {
                message("Convergence warning caught: ", conditionMessage(w))
                converged <<- FALSE
            }
        })

        if (converged) {
            single_protein[, predicted := predict(survival_fit, newdata = .SD)]
        } else {
            single_protein[, predicted := NA_real_]
        }

        if (is_labeled_reference) {
            single_protein[, predicted := ifelse(censored & is_labeled_ref == FALSE, predicted, NA)]
            single_protein[, newABUNDANCE := ifelse(censored & is_labeled_ref == FALSE, predicted, newABUNDANCE)]
        } else {
            single_protein[, predicted := ifelse(censored, predicted, NA)]
            single_protein[, newABUNDANCE := ifelse(censored, predicted, newABUNDANCE)]
        }
        survival = single_protein[, intersect(c(cols, "LABEL", "predicted"), colnames(single_protein)), with = FALSE]
    } else {
        survival = single_protein[, intersect(c(cols, "LABEL"), colnames(single_protein)), with = FALSE]
        survival[, predicted := NA]
    }

    single_protein = .isSummarizable(single_protein, remove50missing)
    if (is.null(single_protein)) {
        return(list(NULL, NULL))
    } else {
        single_protein = single_protein[!is.na(newABUNDANCE), ]
        result = .runTukey(single_protein, is_labeled_reference, censored_symbol,
                           remove50missing)
        if (!is.null(result) && !is.element("LABEL", colnames(result))) {
            result[, LABEL := "L"]
        }
    }
    list(result, survival)
}
