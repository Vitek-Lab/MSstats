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

.getContrastLabels = function(constrasts) {
    prelabels = lapply(contrasts,
                       function(x) sapply(x,
                                          function(y)
                                              paste(y, sep = ",",
                                                    collapse = ",")))
    labels = sapply(prelabels, function(x) paste(x, sep = " vs ",
                                                 collapse = " vs "))
    labels
}


.getSubjectTable = function(patients, contrast_matrix) {
    patients$Value = contrast_matrix[patients$GROUP]
    patients$GROUP = paste0("GROUP", patients$GROUP)
    patients$SUBJECT = paste0("SUBJECT", patients$SUBJECT)
    patients
}

.makeContrastFreeSingle = function(input, contrast_matrix, coefs) {
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
        group_term = contrast_matrix[, as.numeric(levels(input$GROUP))] # GROUP_ORIGINAL?
        group_term = group_term[-1]
        names(group_term) = group
    } else {
        group_term = NULL
    }

    contrast = c(intercept_term, group_term)
    contrast[!is.na(coefs)]
}


.getContrast = function(input, fit, contrast_matrix, coefs, model_data) {
    if (is.element("SUBJECT", colnames(model_data))) {
        model_data = unique(model_data[, c("GROUP", "SUBJECT")])
        patients = .getSubjectTable(model_data, contrast_matrix)
        patient_count = tapply(patients$SUBJECT, patients$GROUP,
                               data.table::uniqueN)
        get_contrast_free(input, patients, patient_count, contrast_matrix, coefs)
    } else {
        .makeContrastFreeSingle(input, contrast_matrix, coefs)
    }
}
