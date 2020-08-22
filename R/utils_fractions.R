.countCommonFeatures = function(features_1, features_2) {
    data.table::uniqueN(intersect(as.character(features_1), as.character(features_2)))
}

.checkMultiRun = function(input) {
    if (is.element("FRACTION", colnames(input))) {
        return(list(has_fractions = TRUE, is_risky = FALSE))
    } else {
        count_techreps = input[, list(n_techreps = data.table::uniqueN(RUN)),
                               by = c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")]
        if (!any(count_techreps$n_techreps > 1)) {
            has_fractions = FALSE
            is_risky = FALSE
        } else {
            info = unique(input[, list(GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN)])
            info[, condition := paste(GROUP_ORIGINAL, SUBJECT_ORIGINAL, sep = "_")]
            single_sample = unique(info[condition == unique(condition)[1],
                                 list(GROUP_ORIGINAL, SUBJECT_ORIGINAL)])
            single_sample_data = input[!is.na(ABUNDANCE) & 
                                           GROUP_ORIGINAL == single_sample$GROUP_ORIGINAL &
                                           SUBJECT_ORIGINAL == single_sample$SUBJECT_ORIGINAL]
            single_run_features = unique(single_sample_data[RUN == unique(RUN)[1], 
                                                            as.character(FEATURE)])
            common = single_sample_data[, list(n_common = .countCommonFeatures(FEATURE, single_run_features)), 
                                        by = "RUN"]
            common$fraction = common$n_common / max(common$n_common)
            overlap = common$fraction[-1]
            
            if (all(overlap > 0.5)) {
                has_fractions = FALSE
                is_risky = FALSE
            } else if (all(overlap < 0.5)) {
                has_fractions = TRUE
                is_risky = FALSE
            } else {
                has_fractions = FALSE
                is_risky = TRUE
            }
            
        }
    }
    list(has_fractions = has_fractions,
         is_risky = is_risky)
}


.handleFractions = function(input, multi_run_check) {
    if (multi_run_check$is_risky) {
        msg = paste("** MSstats suspects that there are fractionations and",
                    "potentially technical replicates too.",
                    "Please add Fraction column to the input.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    } else if (multi_run_check$has_fractions) {
        if (!is.element("FRACTION", colnames(input))) {
            input = .addFractions(input)
            has_missing_fractions = any(is.na(input$FRACTION))
            if (has_missing_fractions) {
                msg = paste("** It is hard to find the same fractionation across sample,",
                            "due to lots of overlapped features between fractionations.",
                            "Please add Fraction column in input.")
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            }
        }
        n_fractions = data.table::uniqueN(input$FRACTION)
        msg = paste("Multiple fractionations exist:",
                    n_fractions, "fractionations per MS replicate.")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        
        input = .removeOverlappingFeatures(input)
        input = .checkOverlappedFeatures(input)
        # TODO: check if there are still overlapped fractions here, stop if so?
    } else {
        input$FRACTION = 1L
    }
    input
}

.addFractions = function(input) {
    input$FRACTION <- NA
    run_info = unique(input[, list(GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN,
                                   CONDITION = paste(GROUP_ORIGINAL, SUBJECT_ORIGINAL, sep = "_"))])
    single_condition = run_info[CONDITION == unique(CONDITION)[1], ]
    single_condition$FRACTION = 1:nrow(single_condition)
    for (run_id in seq_along(single_condition$RUN)) {
        features_in_run = as.character(unique(input[!is.na(ABUNDANCE) & RUN == single_condition$RUN[run_id], FEATURE]))
        same_features = input[!is.na(ABUNDANCE) & FEATURE %in% features_in_run, 
                              list(GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, ABUNDANCE)]
        count_features = data.table::dcast(same_features, RUN ~ GROUP_ORIGINAL + SUBJECT_ORIGINAL,
                                           fun.aggregate = length, value.var = "ABUNDANCE")
        same_fraction = apply(count_features[, -1], 2, function(x) as.character(count_features[which.max(x)," RUN"]))
        input[RUN %in% same_fraction, "FRACTION"] = single_condition[RUN == single_condition$RUN[run_id], FRACTION]
    }
    input
}

.removeOverlappingFeatures = function(input) {
    if (data.table::uniqueN(input$Fraction) > 1) {
        input[, fraction_keep := .getCorrectFraction(.SD), 
              by = "FEATURE", 
              .SDcols = c("FEATURE", "FRACTION", "RUN", "ABUNDANCE")] # by = c("LABEL", "PROTEINNAME", "FEATURE")?
        # input = input[FRACTION == fraction_keep, ] # ?
        input$ABUNDANCE[input$FRACTION != input$fraction_keep] = NA
        input[, !(colnames(input) == "fraction_keep"), with = FALSE]
    }
    input
}


.getCorrectFraction = function(input) {
    measurement_count = input[!is.na(ABUNDANCE) & ABUNDANCE > 0, 
                              list(n_obs = data.table::uniqueN(RUN)),
                              by = "FRACTION"]
    which_max_measurements = which.max(measurement_count$n_obs)
    if (length(which_max_measurements) == 1L) {
        return(unique(measurement_count$n_obs[which_max_measurements]))
    } else {
        average_abundance = input[!is.na(ABUNDANCE) & ABUNDANCE > 0, 
                                  list(mean_abundance = mean(ABUNDANCE)),
                                  by = "FRACTION"]
        which_max_abundance = which.max(average_abundance$mean_abundance)
        unique(average_abundance$mean_abundance[which_max_abundance])
    }
}

.checkOverlappedFeatures = function(input) {
    count_fractions = input[, 
                            list(n_fractions = data.table::uniqueN(FRACTION)),
                            by = "FEATURE"]
    count_fractions = count_fractions[n_fractions > 1, ]
    if (nrow(count_fractions) > 1) {
        overlapped_features = unique(as.character(count_fractions$FEATURE))
        msg = paste("Those features are measured across all fractionations.",
                    "Please keep only one intensity of listed features", 
                    "among fractionations from one sample.\n",
                    paste(overlapped_features, sep = ", ", collapse = ", "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    } else {
        input    
    }
}