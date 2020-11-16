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
    # TODO: add #peaks
    
    summarized_results = data.table::rbindlist(summarized_results)
    list(input, summarized_results)
}


.summarizeLinearSingleProtein = function(input, equal_variances = TRUE) {
    ABUNDANCE = NULL
    
    label = data.table::uniqueN(input$LABEL) > 1
    input = input[!is.na(ABUNDANCE)]
    counts = xtabs(~ RUN + FEATURE, 
                   data = unique(input[, .(FEATURE, RUN)]))
    counts = as.matrix(counts)
    
    is_single_feature = .checkSingleFeature(input)
    
    # TODO: message about i of n proteins, after adding n parameter
    fit = try(.fitLinearModel(input, is_single_feature, 
                              is_single_subject, 
                              has_techreps, is_labeled = label, 
                              equal_variances), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
        msg = paste("*** error : can't fit the model for ", unique(input$PROTEIN))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        result = NULL
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
    }
    
    result
}

.fitLinearModel = function(input, is_single_feature, is_labeled,
                           equal_variances) {
    if (!is_labeled) {
        if (is_single_feature) {
            linear_model = lm(ABUNDANCE ~ RUN , data = input)
        } else {
            linear_model = lm(ABUNDANCE ~ FEATURE + RUN , data = input)
        }
    } else {
        if (is_single_feature) {
            linear_model = lm(ABUNDANCE ~ RUN + ref , data = input)
        } else {
            linear_model = lm(ABUNDANCE ~ FEATURE + RUN + ref, data = input)
        }
    }
    
    ## make equal variance for feature : need to be updated
    if (!equal_variances) {
        linear_model = .updateUnequalVariances(input = input, 
                                               fit = linear_model,
                                               num_iter = 1)
    }
    
    linear_model
}



.updateUnequalVariances = function(input, fit, num_iter) {
    for (i in 1:num_iter) {
        if (i == 1) {
            if (class(fit) == "lm") {
                abs.resids = data.frame(abs.resids = abs(fit$residuals))
                fitted = data.frame(fitted = fit$fitted.values)
            } else {
                abs.resids = data.frame(abs.resids = abs(resid(fit)))
                fitted = data.frame(fitted = fitted(fit))
            }
            input = data.frame(input, 
                               "abs.resids" = abs.resids, 
                               "fitted" = fitted)
        }
        fit.loess = loess(abs.resids ~ fitted, data = input)
        loess.fitted = data.frame(loess.fitted = fitted(fit.loess))
        input = data.frame(input, "loess.fitted" = loess.fitted)
        
        ## loess fitted valuaes are predicted sd
        input$weight = 1 / (input$loess.fitted ^ 2)
        input = input[, -which(colnames(input) %in% "abs.resids")]
        
        ## re-fit using weight
        if (class(fit) == "lm") {
            wls.fit = lm(formula(fit), data = input, weights = weight)
        } else {
            wls.fit = lmer(formula(fit), data = input, weights = weight)
        }
        
        ## lm or lmer
        if (class(wls.fit) == "lm") {
            abs.resids = data.frame(abs.resids = abs(wls.fit$residuals))
        } else {
            abs.resids = data.frame(abs.resids = abs(resid(wls.fit)))
        }
        input = data.frame(input, "abs.resids" = abs.resids)
        
        input = input[, -which(colnames(input) %in% c("loess.fitted", "weight"))]
    }
    
    return(wls.fit)
}


.checkSingleFeature = function(input) {
    data.table::uniqueN(input$FEATURE) < 2
}
