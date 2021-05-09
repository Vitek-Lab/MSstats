#' Feature-level data summarization
#' 
#' @param input list of processed feature-level data
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
#' @param remove_uninformative_feature_outlier It only works after users used featureSubset = "highQuality" 
#' in dataProcess. TRUE allows to remove 1) the features are flagged in the column, 
#' feature_quality = "Uninformative" which are features with bad quality, 
#' 2) outliers that are flagged in the column, is_outlier = TRUE, 
#' for run-level summarization. FALSE (default) uses all features and intensities 
#' for run-level summarization.
#' 
#' @importFrom data.table uniqueN
#' 
#' @export
#' 
MSstatsSummarize = function(proteins_list, method, impute, censored_symbol,
                            remove50missing, equal_variance) {
    num_proteins = length(proteins_list)
    summarized_results = vector("list", num_proteins)
    if (method == "TMP") {
        pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
        for (protein_id in 1:num_proteins) {
            single_protein = proteins_list[[protein_id]]
            summarized_results[[protein_id]] = MSstatsSummarizeSingleTMP(
                single_protein, impute, censored_symbol, remove50missing)
            setTxtProgressBar(pb, protein_id)
        }
        close(pb)
    } else {
        pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
        for (protein_id in 1:num_proteins) {
            single_protein = proteins_list[[protein_id]]
            summarized_result = MSstatsSummarizeSingleLinear(single_protein,
                                                             equal_variances)
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
#' @export
#' 
MSstatsSummarizeSingleLinear = function(single_protein, equal_variances = TRUE) {
    ABUNDANCE = RUN = FEATURE = NULL
    
    label = data.table::uniqueN(single_protein$LABEL) > 1
    single_protein = single_protein[!is.na(ABUNDANCE)]
    single_protein[, RUN := factor(RUN)]
    single_protein[, FEATURE := factor(FEATURE)]
    
    counts = xtabs(~ RUN + FEATURE, 
                   data = unique(single_protein[, .(FEATURE, RUN)]))
    counts = as.matrix(counts)
    is_single_feature = .checkSingleFeature(single_protein)
    
    # TODO: message about i of n proteins, after adding n parameter
    fit = try(.fitLinearModel(single_protein, is_single_feature, is_labeled = label, 
                              equal_variances), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
        msg = paste("*** error : can't fit the model for ", unique(single_protein$PROTEIN))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        result = NULL
    } else {
        if (class(fit) == "lm") {
            cf = summary(fit)$coefficients[, 1]
        } else{
            cf = fixef(fit)
        }
        
        result = unique(single_protein[, .(Protein = PROTEIN, RUN = RUN)])
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
#' @export
#' 
MSstatsSummarizeSingleTMP = function(single_protein, impute, censored_symbol, 
                                     remove50missing) {
    newABUNDANCE = NULL
    cols = intersect(colnames(single_protein), c("newABUNDANCE", "cen", "RUN",
                                                 "FEATURE", "ref"))
    single_protein = single_protein[(n_obs > 1 & !is.na(n_obs)) &
                                        (n_obs_run > 0 & !is.na(n_obs_run))]
    if (nrow(single_protein) == 0) {
        return(list(NULL, NULL))
    }
    single_protein[, RUN := factor(RUN)]
    single_protein[, FEATURE := factor(FEATURE)]
    if (impute) {
        survival_fit = .fitSurvival(single_protein[LABEL == "L", cols,
                                                   with = FALSE])
        single_protein[, predicted := predict(survival_fit,
                                              newdata = .SD)]
        single_protein[, predicted := ifelse(censored & (LABEL == "L"), predicted, NA)]
        single_protein[, newABUNDANCE := ifelse(censored & LABEL == "L",
                                                predicted, newABUNDANCE)]
        survival = single_protein[, c(cols, "predicted"), with = FALSE]
    } else {
        survival = NULL
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
