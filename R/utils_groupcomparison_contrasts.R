#' Create a contrast matrix for groupComparison function
#' 
#' @param contrast list of lists. Each sub-list consists of two vectors that 
#' name conditions that will be compared. See the details section for more 
#' information
#' @param conditions unique condition labels  
#' @param labels labels for contrasts (row.names of the contrast matrix)
#' 
#' @export
#' 
MSstatsContrastMatrix = function(contrasts, conditions, labels = NULL) {
    num_conditions = length(conditions)
    contrast_matrix = matrix(0, nrow = length(contrasts),
                             ncol = num_conditions)
    for (contrast_id in seq_along(contrasts)) {
        contrast = contrasts[[contrast_id]]
        contrast_vector = rep(0, num_conditions)
        positive = conditions %in% contrast[[1]]
        negative = conditions %in% contrast[[2]]
        contrast_vector[positive] = 1 / sum(positive)
        contrast_vector[negative] = -1 / sum(negative)
        contrast_matrix[contrast_id, ] = contrast_vector
    }
    if (is.null(labels)) {
        row.names(contrast_matrix) = .getContrastLabels(contrasts)
    } else {
        row.names(contrast_matrix) = labels
    }
    colnames(contrast_matrix) = conditions
    contrast_matrix
}


#' Get labels for contrasts 
#' @param contrasts list of lists of condition labels
#' @keywords internal
.getContrastLabels = function(contrasts) {
    prelabels = lapply(contrasts,
                       function(x) sapply(x,
                                          function(y)
                                              paste(y, sep = ",",
                                                    collapse = ",")))
    labels = sapply(prelabels, function(x) paste(x, sep = " vs ",
                                                 collapse = " vs "))
    labels
}


#' Create a contrast for a model with only group as a fixed effect
#' @param input summarized data for a single protein
#' @param contrast_matrix row of a contrast_matrix
#' @param coefs coefficients of a linear model (named vector)
#' @keywords internal
.getContrast = function(input, contrast_matrix, coefs) {
    coef_names = names(coefs)
    intercept = grep("Intercept", coef_names, value = TRUE)
    if (length(intercept) > 0) {
        intercept_term = rep(0, length(intercept))
        names(intercept_term) = intercept
    } else {
        intercept_term = NULL
    }
    group = grep("GROUP", coef_names, value = TRUE)
    interaction = grep(":", coef_names, value = TRUE)
    group = setdiff(group, interaction)
    if (length(group) > 0) {
        # THE LINE BELOW WILL REQUIRE CHANGE WHEN SWITCHING TO V4
        group_term = contrast_matrix[, as.numeric(levels(input$GROUP))]
        group_term = group_term[-1]
        names(group_term) = group
    } else {
        group_term = NULL
    }
    contrast = c(intercept_term, group_term)
    contrast[!is.na(coefs)]
}
