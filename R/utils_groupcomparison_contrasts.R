#' @export
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


#' @importFrom utils combn
#' @export
MSstatsContrastMatrix.character = function(contrasts, conditions, labels = NULL) {
    if (contrasts == "pairwise") {
        contrast_combinations = combn(conditions, 2)
        num_combinations = ncol(contrast_combinations)
        contrasts_list = lapply(seq_len(num_combinations),
                                function(x) list(contrast_combinations[1, x],
                                                 contrast_combinations[2, x]))
        MSstatsContrastMatrix.list(contrasts_list, conditions)
    } else {
        stop(paste("Contrast matrix of type", contrasts, "not implemented"))
    }
}


#' @export
MSstatsContrastMatrix.matrix = function(contrasts, conditions, labels = NULL) {
    groups_matrix = colnames(contrasts)
    checkmate::assertSetEqual(groups_matrix, as.character(conditions), 
                              .var.name = "colnames of contrast matrix")
    if (is.null(row.names(contrasts))) {
        stop("Row names of the contrast matrix must be the contrast labels")
    }
    contrasts
}


#' @export
MSstatsContrastMatrix.data.frame = function(contrasts, conditions, labels) {
    groups_matrix = colnames(contrasts)
    checkmate::assertSetEqual(groups_matrix, as.character(conditions), 
                              .var.name = "colnames of contrast matrix")
    if (is.null(row.names(contrasts))) {
        stop("Row names of the contrast matrix must be the contrast labels")
    }
    as.matrix(contrasts)
}


#' Create a contrast matrix for groupComparison function
#' 
#' @param contrasts One of the following:
#' i) list of lists. Each sub-list consists of two vectors that name 
#' conditions that will be compared. See the details section for more information
#' ii) matrix. In this case, it's correctness will be checked
#' iii) "pairwise". In this case, pairwise comparison matrix will be generated
#' iv) data.frame. In this case, input will be converted to matrix
#' @param conditions unique condition labels  
#' @param labels labels for contrasts (row.names of the contrast matrix)
#' 
#' @export
#' 
MSstatsContrastMatrix = function(contrasts, conditions, labels = NULL) {
    UseMethod("MSstatsContrastMatrix", contrasts)
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
