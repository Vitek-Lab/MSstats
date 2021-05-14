#' Normalize MS data
#' 
#' @param input data.table in MSstats format
#' @param normalization_method name of a chosen normalization method: "NONE" or
#' "FALSE" for no normalization, "EQUALIZEMEDIANS" for median normalization,
#' "QUANTILE" normalization for quantile normalization from `preprocessCore` package,
#' "GLOBALSTANDARDS" for normalization based on selected peptides or proteins.
#' @param peptides_dict `data.table` of names of peptides and their corresponding 
#' features. 
#' @param standards character vector with names of standards, required if 
#' "GLOBALSTANDARDS" method was selected.
#' 
#' @export
#' 
#' @return data.table
#' 
#' @examples
#' raw = DDARawData 
#' method = "TMP"
#' cens = "NA"
#' impute = TRUE
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input = MSstatsPrepareForDataProcess(raw, 2, NULL)
#' input = MSstatsNormalize(input, "EQUALIZEMEDIANS") # median normalization
#' head(input)
#' 
MSstatsNormalize = function(input, normalization_method, peptides_dict = NULL, standards = NULL) {
    normalization_method = toupper(normalization_method)
    if (normalization_method == "NONE" | normalization_method == "FALSE") {
        return(input)
    } else if (normalization_method == "EQUALIZEMEDIANS") {
        input = .normalizeMedian(input)
    } else if (normalization_method == "QUANTILE") {
        input = .normalizeQuantile(input)
    } else if (normalization_method == "GLOBALSTANDARDS") {
        input = .normalizeGlobalStandards(input, peptides_dict, standards)
    }
    input
}


#' Median normalization
#' @param input `data.table` in standard MSstats format
#' @importFrom stats median
#' @keywords internal
.normalizeMedian = function(input) {
    ABUNDANCE_RUN = ABUNDANCE_FRACTION = ABUNDANCE = NULL
    
    if (length(unique(input$LABEL)) == 1L) {
        label = "L"
    } else {
        label = "H"
    }
    input[, ABUNDANCE_RUN := .getMedian(.SD, label),
          by = c("RUN", "FRACTION"), .SDcols = c("ABUNDANCE", "LABEL")]
    input[, ABUNDANCE_FRACTION := median(ABUNDANCE_RUN, na.rm = TRUE),
          by = "FRACTION"]
    input[, ABUNDANCE := ABUNDANCE - ABUNDANCE_RUN + ABUNDANCE_FRACTION]
    input = input[, !(colnames(input) %in% c("ABUNDANCE_RUN", "ABUNDANCE_FRACTION")),
                  with = FALSE]
    getOption("MSstatsLog")("Normalization based on median: OK")
    input
}


#' Get median of protein abundances for a given label
#' @param df `data.table`
#' @param label "L" for light isotopes, "H" for heavy isotopes.
#' @importFrom stats median
#' @keywords internal
.getMedian = function(df, label) {
    median(df$ABUNDANCE[df$LABEL == label], na.rm = TRUE)
}


#' Quantile normalization based on the `preprocessCore` package
#' @param input `data.table` in MSstats standard format
#' @keywords internal
.normalizeQuantile = function(input) {
    ABUNDANCE = FRACTION = RUN = LABEL = GROUP_ORIGINAL = NULL
    SUBJECT_ORIGINAL = PROTEIN = PEPTIDE = TRANSITION = INTENSITY = NULL
    
    input[ABUNDANCE == 0, "ABUNDANCE"] = 1
    fractions = unique(input$FRACTION)
    is_labeled = data.table::uniqueN(input$LABEL) > 1
    per_fraction_normalized = vector("list", length(fractions))
    
    if (!is_labeled) {
        grouping_cols = c("FEATURE", "RUN")
        for (fraction_id in fractions) {
            fraction_runs = as.character(unique(input[FRACTION == fractions[fraction_id], RUN]))
            wide = .getWideTable(input, fraction_runs)
            normalized = .quantileNormalizationSingleLabel(wide, fraction_runs)
            normalized = data.table::melt(normalized, id.vars = "FEATURE",
                                          variable.name = "RUN", value.name = "ABUNDANCE")
            per_fraction_normalized[[fraction_id]] = normalized
        }
    } else {
        grouping_cols = c("LABEL", "FEATURE", "RUN")
        for (fraction_id in fractions) {
            fraction_runs = as.character(unique(input[FRACTION == fractions[fraction_id], RUN]))
            wide_h = .getWideTable(input, fraction_runs, "H", FALSE)
            normalized_h = .quantileNormalizationSingleLabel(wide_h,
                                                             fraction_runs, "H")
            normalized_h$LABEL = "H"
            wide_l = .getWideTable(input, fraction_runs)
            for (run_col in fraction_runs) {
                wide_l[[run_col]] = wide_l[[run_col]] - wide_h[[run_col]] + normalized_h[[run_col]]
            }
            wide_l$LABEL = "L"
            per_fraction_normalized[[fraction_id]] = data.table::melt(
                rbind(normalized_h, wide_l), id.vars = c("LABEL", "FEATURE"),
                variable.name = "RUN", value.name = "ABUNDANCE")
        }
    }
    
    per_fraction_normalized = data.table::rbindlist(per_fraction_normalized)
    input = merge(input[, colnames(input) != "ABUNDANCE", with = FALSE],
                  per_fraction_normalized, by = grouping_cols)
    input = input[order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL,
                        RUN, PROTEIN, PEPTIDE, TRANSITION), ]
    input[!is.na(INTENSITY) & INTENSITY == 1, "ABUNDANCE"] = 0 # Skyline
    
    if (length(fractions) == 1L) {
        msg = "Normalization : Quantile normalization - okay"
    } else {
        msg = "Normalization : Quantile normalization per fraction - okay"
    }
    getOption("MSstatsLog")("INFO", msg)
    
    input
}


#' Utility function for quantile normalization - get table in wide format
#' @param input `data.table` in MSstats standard format
#' @param vector of run labels
#' @param label "L" for light isotopes, "H" for heavy isotopes
#' @param remove_missing if TRUE, only non-missing values will be considered
#' @keywords internal
.getWideTable = function(input, runs, label = "L", remove_missing = TRUE) {
    RUN = NULL
    
    if (remove_missing) {
        nonmissing_filter = !is.na(input$INTENSITY)
    } else {
        nonmissing_filter = rep(TRUE, nrow(input))
    }
    label_filter = input$LABEL == label
    
    wide = data.table::dcast(input[nonmissing_filter & label_filter & (RUN %in% runs)],
                             FEATURE ~ RUN, value.var = "ABUNDANCE")
    wide = wide[, lapply(.SD, .replaceZerosWithNA)]
    colnames(wide)[-1] = runs
    wide
}


#' Quantile normalization for a single label
#' @param input `data.table` in MSstats standard format
#' @param runs run labels
#' @param label "L" for light isotopes, "H" for heavy isotopes
#' @importFrom preprocessCore normalize.quantiles
#' @keywords internal
.quantileNormalizationSingleLabel = function(input, runs, label = "L") {
    FEATURE = NULL
    
    normalized = input[, list(FEATURE, preprocessCore::normalize.quantiles(as.matrix(.SD))),
                       .SDcols = runs]
    colnames(normalized)[-1] = runs
    normalized
}


#' Utility function for normalization: replace 0s by NA
#' @param vec vector
#' @keywords internal
.replaceZerosWithNA = function(vec) {
    vec = unlist(vec, FALSE, FALSE)
    if (is.character(vec) | is.factor(vec)) {
        vec
    } else {
        ifelse(vec == 0, NA, vec)
    }
}


#' Normalization based on standards
#' @inheritParams MSstatsNormalize
#' @importFrom data.table melt uniqueN
#' @keywords internal
.normalizeGlobalStandards = function(input, peptides_dict, standards) {
    PeptideSequence = PEPTIDE = PROTEIN = median_by_fraction = NULL
    Standard = FRACTION = LABEL = ABUNDANCE = RUN = GROUP = NULL
    
    proteins = as.character(unique(input$PROTEIN))
    means_by_standard = unique(input[, list(RUN)])
    for (standard_id in seq_along(standards)) {
        peptide_name = unlist(peptides_dict[PeptideSequence == standards[standard_id],
                                            as.character(PEPTIDE)], FALSE, FALSE)
        if (length(peptide_name) > 0) {
            standard = input[PEPTIDE == peptide_name, ]
        } else {
            if (standards[standard_id] %in% proteins) {
                standard = input[PROTEIN == standards[standard_id], ]
            } else {
                msg = paste("global standard peptides or proteins, ",
                            standards[standard_id],", is not in dataset.",
                            "Please check whether 'nameStandards' input is correct or not.")
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            }
        }
        mean_by_run = standard[GROUP != "0" & !is.na(ABUNDANCE),
                               list(mean_abundance = mean(ABUNDANCE, na.rm = TRUE)),
                               by = "RUN"]
        colnames(mean_by_run)[2] = paste0("meanStandard", standard_id)
        means_by_standard = merge(means_by_standard, mean_by_run,
                                  by = "RUN", all.x = TRUE)
    }
    means_by_standard = data.table::melt(means_by_standard, id.vars = "RUN",
                                         variable.name = "Standard", value.name = "ABUNDANCE")
    means_by_standard[, mean_by_run := mean(ABUNDANCE, na.rm = TRUE), by = "RUN"]
    means_by_standard = merge(means_by_standard, unique(input[, list(RUN, FRACTION)]),
                              by = "RUN")
    means_by_standard[, median_by_fraction := median(mean_by_run, na.rm = TRUE),
                      by = "FRACTION"]
    means_by_standard[, ABUNDANCE := NULL]
    means_by_standard[, Standard := NULL]
    means_by_standard = unique(means_by_standard)
    
    input = merge(input, means_by_standard, all.x = TRUE, by = c("RUN", "FRACTION"))
    input[, ABUNDANCE := ifelse(LABEL == "L", ABUNDANCE - mean_by_run + median_by_fraction, ABUNDANCE)]
    
    if (data.table::uniqueN(input$FRACTION) == 1L) {
        msg = "Normalization : normalization with global standards protein - okay"
    } else {
        msg = "Normalization : normalization with global standards protein - okay"
    }
    getOption("MSstatsLog")("INFO", msg)
    input[ , !(colnames(input) %in% c("mean_by_run", "median_by_fraction")), with = FALSE]
}


#' Re-format the data before feature selection
#' 
#' @param input `data.table` in MSstats format
#' 
#' @importFrom data.table uniqueN
#' 
#' @export
#' 
#' @return data.table
#' 
#' @examples 
#' raw = DDARawData 
#' method = "TMP"
#' cens = "NA"
#' impute = TRUE
#' MSstatsConvert::MSstatsLogsSettings(FALSE)
#' input = MSstatsPrepareForDataProcess(raw, 2, NULL)
#' input = MSstatsNormalize(input, "EQUALIZEMEDIANS")
#' input = MSstatsMergeFractions(input)
#' head(input)
#' 
MSstatsMergeFractions = function(input) {
    ABUNDANCE = INTENSITY = GROUP_ORIGINAL = SUBJECT_ORIGINAL = RUN = NULL
    originalRUN = FRACTION = TECHREPLICATE = tmp = merged = newRun = NULL
    ncount = FEATURE = NULL
    
    input[!is.na(ABUNDANCE) & ABUNDANCE < 0, "ABUNDANCE"] = 0
    input[!is.na(INTENSITY) & INTENSITY == 1, "ABUNDANCE"] = 0
    if (data.table::uniqueN(input$FRACTION) == 1) {
        return(input)
    } else {
        if (is.element("TECHREPLICATE", colnames(input))) {
            run_info = unique(input[, 
                                    list(GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, 
                                         originalRUN, FRACTION, TECHREPLICATE)])
            match_runs = try(
                data.table::dcast(run_info, 
                                  GROUP_ORIGINAL + SUBJECT_ORIGINAL + TECHREPLICATE ~ FRACTION,
                                  value.var = "originalRUN"), silent = TRUE
            )
            if (inherits(match_runs, "try-error")) {
                msg = "*** error : can't figure out which multiple runs come from the same sample."
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            } else {
                input$newRun = NA
                input$newRun = as.character(input$newRun)
                run_info[, GROUP_ORIGINAL := as.character(GROUP_ORIGINAL)]
                run_info[, SUBJECT_ORIGINAL := as.character(SUBJECT_ORIGINAL)]
                for (k in seq_len(nrow(run_info))) {
                    input[originalRUN %in% run_info[k, 4:ncol(run_info)], "newRun"] = paste(paste(run_info[k, 1:4], collapse = "_"), 'merged', sep = "_")   
                }
                
                select_fraction = input[!is.na(ABUNDANCE) & ABUNDANCE > 0, 
                                        list(ncount = .N),
                                        by = c("FEATURE", "FRACTION")]
                select_fraction = select_fraction[ncount == data.table::uniqueN(input$newRun)]
                select_fraction[, tmp := paste(FEATURE, FRACTION, sep = "_")]
                input$tmp = paste(input$FEATURE, input$FRACTION, sep="_")
                input = input[!(tmp %in% select_fraction$tmp), ]
                input$originalRUN = input$newRun
                input$RUN = input$originalRUN
                input$RUN = factor(input$RUN, levels=unique(input$RUN), labels=seq(1, length(unique(input$RUN))))
                input = input[, !(colnames(input) %in% c('tmp','newRun')), with = FALSE]
            }
        } else {
            run_info = unique(input[, 
                                    list(GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, 
                                         originalRUN, FRACTION)])
            match_runs = try(
                data.table::dcast(run_info, 
                                  GROUP_ORIGINAL + SUBJECT_ORIGINAL ~ FRACTION,
                                  value.var = "originalRUN"), silent = TRUE
            )
            if (inherits(match_runs, "try-error")) {
                msg = "*** error : can't figure out which multiple runs come from the same sample."
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            } else {
                match_runs[, merged := "merged"]
                match_runs[, newRun := do.call(paste, c(.SD, sep = "_")), 
                           .SDcols = c(1:3, ncol(match_runs))]
                match_runs = unique(match_runs[, list(GROUP_ORIGINAL,
                                                      SUBJECT_ORIGINAL,
                                                      newRun)])
                
                input = merge(input, match_runs,
                              by = c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL"),
                              all.x = TRUE)
                select_fraction = input[!is.na(ABUNDANCE) & input$ABUNDANCE > 0,
                                        list(ncount = .N),
                                        by = c("FEATURE", "FRACTION")]
                select_fraction = select_fraction[ncount != 0]
                select_fraction[, tmp := paste(FEATURE, FRACTION, sep = "_")]
                input$tmp = paste(input$FEATURE, input$FRACTION, sep = "_")
                input = input[tmp %in% select_fraction$tmp, ]
                input$originalRUN = input$newRun
                input$RUN = input$originalRUN
                input$RUN = factor(input$RUN, levels = unique(input$RUN), 
                                   labels = seq_along(unique(input$RUN)))
                input = input[, !(colnames(input) %in% c('tmp','newRun')), 
                              with = FALSE]
            }
        }
    }
    input
}
