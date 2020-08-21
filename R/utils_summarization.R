.runQuantification = function(
    input, summaryMethod, equalFeatureVar, cutoffCensored, censoredInt, 
    remove50missing, MBimpute, original_scale, logsum, featureSubset,
    remove_uninformative_feature_outlier, message.show, clusters) {
    
    label = data.table::uniqueN(input$LABEL) == 2
    
    if (label) {
        input$ref = factor(ifelse(input$LABEL == "L", 
                                  input$RUN[input$LABEL == "L"], 0))
    }
    
    if (is.element("remove", colnames(input))) { # todo: always have "remove" column
        input = input[!(remove), ]
    }
    
    if (remove_uninformative_feature_outlier & 
        is.element("feature_quality", colnames(input))) {
        input[feature_quality == "Uninformative", ABUNDANCE] = NA
        input[(is_outlier), ABUNDANCE] = NA
        msg = "** Filtered out uninformative feature and outliers."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    if (summaryMethod == "linear") {
        result = .summarizeLinear(input, cutoffCensored, censoredInt, 
                                  remove50missing)
    } else if (summaryMethod == "TMP") {
        result = .summarizeTukey(input, MBimpute, cutoffCensored, censoredInt, 
                                 remove50missing, original_scale, clusters)
    }
    result
}
