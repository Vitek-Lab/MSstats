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
    UseMethod("MSstatsContrastMatrix", contrasts)
}

MSstatsContrastMatrix.list = function(contrasts, conditions, labels = NULL) {
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


MSstatsContrastMatrix.character = function(contrasts, conditions, labels = NULL) {
    if (contrasts == "pairwise") {
        contrast_combinations = combn(conditions, 2)
        contrasts_list = lapply(1:ncol(contrast_combinations),
                                function(x) list(contrast_combinations[1, x],
                                                 contrast_combinations[2, x]))
        MSstatsContrastMatrix.list(contrasts_list, conditions)
    } else {
        stop(paste("Contrast matrix of type", contrasts, "not implemented"))
    }
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


