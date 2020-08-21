.summarizeLinear = function(input, cutoff_base, censored_symbol, remove50missing) {
    if (is.null(censored_symbol)) {
        .summarizeLinearNoImputation(input)
    } else {
        .summarizeLinearImputation(input, cutoff_base, censored_symbol, remove50missing)
    }
}

.summarizeLinearNoImputation = function(input) {
    input = input[!is.na(input$ABUNDANCE),]
    proteins = unique(input$PROTEIN)
    n_proteins = data.table::uniqueN(input$PROTEIN)
    
    summarized_result = vector("list", n_proteins)
    for (protein_id in seq_along(unique(input$PROTEIN))) {
        protein = proteins[protein_id]
        single_protein = input[PROTEIN == protein, ]
        single_protein_output = .summarizeLinearSingleProtein(single_protein)
        summarized_result[[protein_id]] = single_protein_output
    }
    
    list(
        input,    
        summarized = data.table::rbindlist(summarized_result)
    )
}

.summarizeLinearSingleProtein = function(input) {
    input$RUN = factor(input$RUN)
    input$FEATURE = factor(input$FEATURE)
    label = data.table::uniqueN(input$LABEL) > 1
    result = .getResultTableTemplate(input, label)
    
    singleFeature <- .checkSingleFeature(input)
    singleSubject <- .checkSingleSubject(input)
    TechReplicate <- .checkTechReplicate(input) ## use for label-free model
    
    # TODO: message about i of n proteins, after adding n parameter
    fit <- try(.fit.quantification.run(input, singleFeature, singleSubject, 
                                       TechReplicate, labeled = label, 
                                       equalFeatureVar), silent=TRUE)
    
    if (inherits(fit, "try-error")) {
        msg = paste("*** error : can't fit the model for ", unique(input$PROTEIN))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        
        # TODO: residuals and fitted not needed anymore - no need to return two data.frames
        if (nrow(input) != 0) {
            input$residuals = NA
            input$fitted = NA
            result = NULL
        }
    } else {
        if (class(fit) == "lm") {
            cf = summary(fit)$coefficients
        } else{
            cf = fixef(fit)
        }
        
        input$LogIntensities = .getQuantificationLinear(input, fit, cf, label)
        input$residuals = resid(fit)
        input$fitted = fitted(fit)
        # TODO: make sure fitted and residuals are in the output. Labeled is tricky
    }
    
    result
}


.getResultTableTemplate = function(input, label) {
    if (label) {
        merge(input[, list(Protein = unique(PROTEIN),
                           NumFeature = data.table::uniqueN(FEATURE),
                           Run = unique(RUN))],
              input[, list(NumPeaks = data.table::.N), by = "RUN"])
    } else {
        input$ref = factor(input$ref)
        counts = input[, list(n_obs = data.table::.N), by = "ref"]
        data.table::data.table(Protein = unique(input$PROTEIN),
                               RUN = c(levels(input$ref)[-1], "Ref"),
                               LogIntensities = NA, 
                               NumFeature = length(unique(input$FEATURE)),
                               NumPeaks = c(counts[-1, n_obs], 
                                            counts[1, n_obs]))
        
    }
    input$LogIntensities = NA
    input
}

.getQuantificationLinear = function(input, fit, cf, label) {
    result = vector("numeric", data.table::uniqueN(input$RUN))
    for (i in seq_along(setdiff(unique(input$RUN), "Ref"))) {
        contrast_matrix = rep(0, nlevels(input$RUN))
        contrast_matrix[i] = 0
        contrast = .make.contrast.run.quantification(fit, contrast_matrix,
                                                     input, label)
        if (inherits(fit, "lm")) {
            result[i] = .estimableFixedQuantification(cf, contrast)
        } else {
            result[i] = .estimableRandomQuantification(cf, contrast)
        }
    }
    
    ## for label-based case, need reference quantification
    if (label) {
        contrast <- .make.contrast.run.quantification.reference(fit, 
                                                                contrast_matrix,
                                                                input)
        if (class(fit) == "lm") {
            result[length(result)] <- .estimableFixedQuantification(cf,contrast)
        } else{
            result[length(result)] <- .estimableRandomQuantification(cf,contrast)
        }
    }
    result
}

.summarizeLinearImputation = function(input, cutoff_base, censored_symbol, 
                                      remove50missing) {
    input$PROTEIN <- factor(input$PROTEIN)
    input$RUN <- factor(input$RUN)
    proteins = unique(input$PROTEIN)
    n_proteins = length(proteins)
    
    summarized_result = vector("list", n_proteins)
    
    for (protein_id in seq_along(unique(input$PROTEIN))) {
        protein = proteins[protein_id]
        single_protein = input[PROTEIN == protein, ]
        single_protein_output = .summarizeLinearSingleProteinImpute(
            single_protein, cutoff_base, censored_symbol, remove50missing)
        summarized_result[[protein_id]] = single_protein_output
    }
    # dataafterfit <- NULL	
    # TODO: survival here?
    
    list(
        input,    
        summarized = data.table::rbindlist(summarized_result)
    )
}


.summarizeLinearSingleProteinImpute = function(input, cutoff_base, censored_symbol,
                                               remove50missing) {
    label = data.table::uniqueN(input$LABEL) > 1
    input$FEATURE <- factor(input$FEATURE)
    # TODO: add message protein i of n after adding n_proteins parameter
    if (all(is.na(input$ABUNDANCE))) {
        msg = paste("Can't summarize for ", unique(input$PROTEIN), 
                    "because all measurements are NAs.")
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        return(NULL)
    }
    
    count_measurements = input[LABEL == "RUN" & !is.na(input$INTENSITY),
                               list(n_obs = .N),
                               by = "RUN"]
    count_measurements = count_measurements[n_obs == 0, ]
    if (nrow(count_measurements) > 0) {
        input = input[!(RUN %in% count_measurements$RUN), ]
    }
    input$RUN = factor(input$RUN)
    if (data.table::uniqueN(input$RUN) == 1L) {
        msg = paste("* Only 1 MS run in", unique(input$PROTEIN), 
                    "has measurement. Can't summarize with censored intensities.")
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        return(NULL)
    }
    count_by_feature = input[LABEL == "L" & !is.na(INTENSITY) & INTENSITY != 0, 
                             list(n_obs = .N), by = "FEATURE"]
    count_by_feature = count_by_feature[n_obs == 0L, ]
    if (nrow(count_by_feature) > 0) {
        input = input[!(FEATURE %in% unique(count_by_feature$FEATURE)), ]
        input$FEATURE = factor(input$FEATURE)
    }		
    if (nrow(input) == 0) {
        message(paste("* All measurements are NAs or only one measurement per feature in",
                      unique(input$PROTEIN), ". Can't summarize with censored intensities."))
        
        return(NULL)
    }	

    if (censored_symbol == "0") {
        input$cen = ifelse(!is.na(input$INTENSITY) & input$INTENSITY == 0, 0, 1)
    } else if (censored_symbol == "NA") {
        input$cen = ifelse(is.na(input$INTENSITY), 0, 1)
    }
    input = .setCensoredByThreshold(input, cutoff_base, censored_symbol,
                                    remove50missing)
    survival_fit = .fitSurvival(input)
    surv_coef = summary(survival_fit)$coefficients
    # TODO: IMPUTATION SEEMS TO BE MISSING!!
    
    result = data.table::data.table(Protein = unique(input$PROTEIN),
                                    RUN = unique(input$RUN))
    
    log_intensities = vector("numeric", data.table::uniqueN(input$RUN))
    for (run_id in seq_along(unique(input$RUN))) {
        contrast_matrix = rep(0, nlevels(input$RUN))
        contrast_matrix[run_id] = 1
        contrast = .make.contrast.run.quantification.Survival(survival_fit, 
                                                              contrast_matrix,
                                                              input, label)
        log_intensities[run_id] = .estimableFixedQuantificationSurvival(surv_coef, contrast)
    }
    result$LogIntensities = log_intensities
    # What was the reshape-melt code for?
    result
    # Add survival predictions to the "input" output [if not imputed, put NA, put NA on H]
    # 
}


.setCensoredByThreshold = function(input, cutoff_base, censored_symbol,
                                   remove50missing) {
    if (censored_symbol == "NA") {
        nonmissing_filter = !is.na(input$ABUNDANCE)
    } else if (censored_symbol == "0") {
        nonmissing_filter = !is.na(input$INTENSITY) & input$INTENSITY != 0
    }
    
    if (cutoff_base == "minFeature") {
        grouping_vars = c("FEATURE", "LABEL")
    } else if (cutoff_base == "minRun") {
        grouping_vars = c("RUN", "LABEL")
    } 
    
    if (cutoff_base %in% c("minFeature", "minRun")) {
        cutoffs = input[nonmissing_filter, 
                        list(ABUNDANCE_cut = 0.99*min(ABUNDANCE)), 
                        by = grouping_vars]
        input = merge(input, cutoffs, by = grouping_vars)
        input$ABUNDANCE = ifelse(nonmissing_filter, input$ABUNDANCE,
                                 input$ABUNDANCE_cut)
        input = input[, !(colnames(input) == "ABUNDANCE_cut"), with = FALSE]
    } else {
        feature_cutoffs = input[nonmissing_filter, 
                                list(ABUNDANCE_cut_fea = 0.99*min(ABUNDANCE)),
                                by = c("FEATURE", "LABEL")]
        
        if (remove50missing & censored_symbol == "0") {
            n_features = data.table::uniqueN(input$FEATURE)
            missing_runs = input[, list(perc_nm = .N / n_features), 
                                 by = "RUN"]
            missing_runs = missing_runs[perc_nm < 0.5, ]
            if (nrow(missing_runs) > 0) {
                input = input[!(RUN %in% unique(missing_runs$RUN)), ]
                input$RUN = factor(input$RUN)
            }
        }
        
        run_cutoffs = input[nonmissing_filter, 
                            list(ABUNDANCE_cut_run = 0.99*min(ABUNDANCE)),
                            by = c("RUN", "LABEL")]
        cutoffs = merge(feature_cutoffs, run_cutoffs, by = c("LABEL"),
                        allow.cartesian = TRUE)
        cutoffs$final_cutoff = sapply(1:nrow(cutoffs), 
                                      function(i) min(cutoffs$ABUNDANCE_cut_fea[i],
                                                      cutoffs$ABUNDANCE_cut_run[i]))
        if (data.table::uniqueN(input[, list(FEATURE, RUN)]) > 1) {
            input$ABUNDANCE = ifelse(nonmissing_filter,
                                     input$ABUNDANCE, cutoffs$final_cutoff)
        } else {
            if (censored_symbol == "0") {
                input$ABUNDANCE = ifelse(nonmissing_filter, input$ABUNDANCE, 
                                         cutoffs$ABUNDANCE_cut_fea)
            }
        }
    }
    input
}

#' @importFrom data.table uniqueN
.fitSurvival = function(input) {
    missingness_filter = !is.na(input$ABUNDANCE)
    n_total = nrow(input[missingness_filter, ])
    n_features = data.table::uniqueN(input[missingness_filter, FEATURE])
    n_runs = data.table::uniqueN(input[missingness_filter, RUN])
    is_labeled = data.table::uniqueN(input$LABEL) > 1
    countdf = n_total  < n_features + n_runs - 1
    
    # TODO: set.seed here?
    set.seed(100)
    if (is_labeled) {
        if (length(unique(input$FEATURE)) == 1) {
            # with single feature, not converge, wrong intercept
            # need to check
            fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN + ref,
                                         data = input, dist = "gaussian")
        } else {
            if (countdf) {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN + ref,
                                             data = input, dist = "gaussian")
            } else {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ FEATURE + RUN + ref,
                                             data = input, dist = "gaussian")
            }
        }
    } else {
        if (n_features == 1L) {
            fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type = "left") ~ RUN,
                                        data = input, dist = "gaussian")    
        } else {
            if (countdf) {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type = "left") ~ RUN,
                                            data = input, dist = "gaussian")
            } else {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type = "left") ~ FEATURE + RUN,
                                            data = input, dist = "gaussian")
            }
        }  
    }
    fit
}
