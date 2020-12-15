#' @importFrom data.table uniqueN
#' @importFrom stats na.omit
.removeSingleLabelFeatures = function(input) {
    n_obs = NULL
    
    if (data.table::uniqueN(input$LABEL) == 2) {
        counts = na.omit(input)[, list(n_obs = .N), 
                                by = c("LABEL", "FEATURE", 
                                       "GROUP_ORIGINAL")]
        counts[, AnyMissing := any(n_obs == 0), by = c("FEATURE", "LABEL")]
        counts = counts[(AnyMissing), ]
        input$SuggestToFilter = ifelse(input$FEATURE %in% unique(counts$FEATURE),
                                       1, 0)
    }
    input 
}
