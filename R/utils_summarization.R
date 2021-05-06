#' Feature-level data summarization
#' 
#' @param input `data.table` in MSstats standard format
#' @param summaryMethod summarization method: "linear" or "TMP" 
#' @param equalFeatureVar only for summaryMethod = "linear". Default is TRUE. 
#' Logical variable for whether the model should account for heterogeneous variation 
#' among intensities from different features. Default is TRUE, which assume equal
#' variance among intensities from features. FALSE means that we cannot assume 
#' equal variance among intensities from features, then we will account for
#' heterogeneous variation from different features.
#' @param cutoffCensored Cutoff value for censoring. only with censoredInt = 'NA' or '0'. 
#' Default is 'minFeature', which uses minimum value for each feature. 
#' 'minFeatureNRun' uses the smallest between minimum value of corresponding 
#' feature and minimum value of corresponding run. 'minRun' uses minimum value for each run.
#' @param censoredInt Missing values are censored or at random. 
#' 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. 
#' '0' uses zero intensities as censored intensity. 
#' In this case, NA intensities are missing at random. 
#' The output from Skyline should use '0'. 
#' Null assumes that all NA intensites are randomly missing.
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the runs 
#' which have more than 50\% missing values. FALSE is default.
#' @param MBimpute only for summaryMethod = "TMP" and censoredInt = 'NA' or '0'. 
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
MSstatsSummarize = function(
    input, summaryMethod = "TMP", equalFeatureVar = TRUE, 
    cutoffCensored = "minFeature", censoredInt = "NA", 
    remove50missing = FALSE, MBimpute = TRUE, 
    remove_uninformative_feature_outlier = FALSE) {
    ABUNDANCE = feature_quality = is_outlier = NULL
    
    input = .removeSingleLabelFeatures(input)
    label = data.table::uniqueN(input$LABEL) == 2
    if (label) {
        input$ref = factor(ifelse(input$LABEL == "L", 
                                  input$RUN[input$LABEL == "L"], 0))
    }
    
    if (is.element("remove", colnames(input))) {
        input = input[!(remove), ]
    }
    
    if (remove_uninformative_feature_outlier & 
        is.element("feature_quality", colnames(input))) {
        input[, ABUNDANCE := ifelse(feature_quality == "Uninformative", 
                                    NA, ABUNDANCE)]
        input[, ABUNDANCE := ifelse(is_outlier, NA, ABUNDANCE)]
        msg = "** Filtered out uninformative features and outliers."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    if (summaryMethod == "linear") {
        result = .summarizeLinear(input, MBimpute, cutoffCensored, 
                                  censoredInt, equalFeatureVar)
    } else if (summaryMethod == "TMP") {
        result = .summarizeTukey(input, MBimpute, cutoffCensored, 
                                 censoredInt, remove50missing)
    }
    msg = " == the summarization per subplot is done."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("\nINFO ", msg)
    result
}
