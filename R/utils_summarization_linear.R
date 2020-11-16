#' Linear summarization
#' @inheritParams .summarizeTukey
#' @return data.table
#' @keywords internal
.summarizeLinear = function(input, impute, cutoff_base, censored_symbol) {
    censored = ABUNDANCE = FEATURE = more50missing = PROTEIN = NULL
    prop_features = RUN = NumImputedFeature = LABEL = INTENSITY = NULL
    
    input[, newABUNDANCE := ABUNDANCE]
    input[, nonmissing := .getNonMissingFilter(.SD, impute, censored_symbol)]
    input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE")]
    # remove feature with 1 measurement
    input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)] 
    input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]
    
    input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
    input[, prop_features := sum(nonmissing) / total_features,
          by = c("PROTEIN", "RUN")] 
    
    proteins = unique(input$PROTEIN)
    n_proteins = length(proteins)
    summarized_results = vector("list", n_proteins)
    pb = utils::txtProgressBar(min = 0, max = n_proteins, style = 3)
    for (protein_id in seq_along(proteins)) {
        single_protein = input[PROTEIN == proteins[protein_id]]
        single_protein = single_protein[(n_obs > 1 & !is.na(n_obs)) &
                                            (n_obs_run > 0 & !is.na(n_obs_run))]
        single_protein[, RUN := factor(RUN)]
        single_protein[, FEATURE := factor(FEATURE)]
        summarized_result = .summarizeLinearSingleProtein(single_protein)
        summarized_results[[protein_id]] = summarized_result
        setTxtProgressBar(pb, protein_id)
    }
    input[, NonMissingStats := .getNonMissingFilterStats(.SD, censored_symbol)]
    input[, NumMeasuredFeature := sum(NonMissingStats), 
          by = c("PROTEIN", "RUN")]
    input[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
    input[, more50missing := MissingPercentage >= 0.5]
    if (!is.null(censored_symbol)) {
        if (is.element("censored", colnames(input))) {
            input[, nonmissing_orig := LABEL == "L" & !censored]
        } else {
            input[, nonmissing_orig := LABEL == "L" & !is.na(INTENSITY)]
        }
        input[, nonmissing_orig := ifelse(is.na(newABUNDANCE), TRUE, nonmissing_orig)]
        input[, NumImputedFeature := 0]
    }
    
    summarized_results = data.table::rbindlist(summarized_results)
    list(input, summarized_results)
}


.summarizeLinearSingleProtein = function(input, equalFeatureVar = TRUE) {
    label = data.table::uniqueN(input$LABEL) > 1
    input = input[!is.na(ABUNDANCE)]
    counts = xtabs(~ RUN + FEATURE, 
                   data = unique(input[, .(FEATURE, RUN)]))
    counts = as.matrix(counts)
    
    is_single_feature = .checkSingleFeature(input)
    is_single_subject = .checkSingleSubject(input)
    has_techreps = .checkTechReplicate(input) ## use for label-free model
    
    # TODO: message about i of n proteins, after adding n parameter
    fit = try(.fit.quantification.run(input, is_single_feature, 
                                      is_single_subject, 
                                      has_techreps, labeled = label, 
                                      equalFeatureVar), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
        msg = paste("*** error : can't fit the model for ", unique(input$PROTEIN))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        
        # # TODO: residuals and fitted not needed anymore - no need to return two data.frames
        # if (nrow(input) != 0) {
        #     input$residuals = NA
        #     input$fitted = NA
        #     result = NULL
        # }
    } else {
        if (class(fit) == "lm") {
            cf = summary(fit)$coefficients[, 1]
        } else{
            cf = fixef(fit)
        }
        
        result = unique(input[, .(Protein = PROTEIN, RUN = RUN)])
        log_intensities = get_linear_summary(input, cf,
                                             counts, label)
        result[, LogIntensities := log_intensities]
        # input$residuals = resid(fit)
        # input$fitted = fitted(fit)
        # TODO: make sure fitted and residuals are in the output. Labeled is tricky
    }
    
    result
}


## check single subject for both case-control and time-course?
.checkSingleSubject = function(input) {
    SUBJECT_ORIGINAL = num_subjects = NULL
    
    all(input[, list(num_subjects = data.table::uniqueN(SUBJECT_ORIGINAL)), 
              by = "GROUP_ORIGINAL"][, num_subjects] == 1)
}

.checkSingleFeature = function(input) {
    data.table::uniqueN(input$FEATURE) < 2
}

.checkTechReplicate = function(input) {
    num_replicates = RUN = NULL
    
    all(input[,
              list(num_replicates = data.table::uniqueN(RUN)),
              by = "SUBJECT_NESTED"][, num_replicates] != 1)
}
