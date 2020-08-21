
.flagCensored = function(input, summary_method, impute, missing_symbol,
                         censored_cutoff) {
    if (summary_method == "TMP" & impute) {
        is_labeled = nlevels(input$LABEL) == 2 # label <- nlevels(input$LABEL) == 2
        input$censored <- FALSE
        ## if intensity = 1, but abundance > cutoff after normalization, it also should be censored.
        if (!is.null(censored_cutoff)) {
            quantiles = input[!is.na(INTENSITY) & INTENSITY > 1 & LABEL == "L", 
                              quantile(ABUNDANCE, 
                                       prob = c(0.01, 0.25, 0.5, 0.76, 
                                                censored_cutoff), 
                                       na.rm = TRUE)]
            iqr = quantiles[4] - quantiles[2]
            multiplier = (quantiles[5] - quantiles[4]) / iqr
            cutoff_lower = (quantiles[2] - multiplier * iqr) 
            input$censored = !is.na(input$INTENSITY) & 
                input$LABEL == "L" &
                input$ABUNDANCE < cutoff_lower
            if (cutoff_lower <= 0 & !is.null(missing_symbol) & missing_symbol == "0") {
                zero_one_filter = (!is.na(input$INTENSITY) & input$INTENSITY == 1) |
                    (!is.na(input$ABUNDANCE & input$ABUNDANCE) <= 0)
                input$censored = ifelse(zero_one_filter, TRUE, input$censored)
            }
            if (!is.null(missing_symbol) & missing_symbol == "NA") {
                input$censored = ifelse(is.na(input$INTENSITY), TRUE, 
                                        input$censored)
            }
            
            if (!is_labeled) {
                msg = paste('** Log2 intensities under cutoff =', 
                            format(cutoff_lower, digits = 5), 
                            ' were considered as censored missing values.')
                msg_2 = paste("** Log2 intensities =", missing_symbol, "were considered as censored missing values.")
            } else {
                msg = paste('** Log2 endogenous intensities under cutoff =', 
                      format(cutoff_lower, digits = 5), 
                      ' were considered as censored missing values.')
                getOption("MSstatsMsg")("INFO", msg)
                getOption("MSstatsMsg")("INFO", msg)
            }
            getOption("MSstatsMsg")("INFO", msg)
            getOption("MSstatsMsg")("INFO", msg_2)
        } else {
            if (missing_symbol == '0') {
                input$censored = input$LABEL == "L" & 
                    !is.na(input$INTENSITY) &
                    (input$INTENSITY == 1 | input$ABUNDANCE <= 0)
            } else if (missing_symbol == 'NA') {
                input$censored = input$LABEL == "L" & is.na(input$ABUNDANCE)
            }
        }
    }
    input
}
