.makeBalancedDesign = function(input, fill_missing) {
    n_rows = nrow(input)
    n_full = data.table::uniqueN(input$LABEL) * data.table::uniqueN(input$FEATURE) * data.table::uniqueN(input$RUN)
    if (n_rows < n_full) {
        if (fill_missing) {
            all_possibilities = data.table::as.data.table(
                expand.grid(LABEL = unique(input[, LABEL]),
                            FEATURE = unique(input[, FEATURE]),
                            RUN = unique(input[, RUN])))
            all_possibilities = merge(all_possibilities,
                                      unique(input[, list(FEATURE, LABEL)]),
                                      all.x = TRUE, by = c("FEATURE", "LABEL"))
            all_possibilities = merge(
                all_possibilities, 
                unique(input[, list(PROTEIN, PEPTIDE, FEATURE, TRANSITION, FRACTION)]), 
                all.x = TRUE, by = c("FEATURE"))
            all_possibilities = merge(
                all_possibilities, 
                unique(input[, list(RUN, GROUP, GROUP_ORIGINAL, SUBJECT,
                                    SUBJECT_ORIGINAL, SUBJECT_NESTED)]),
                all.x = TRUE, by = "RUN")
            input = merge(all_possibilities, input[, list(FEATURE, RUN, LABEL, INTENSITY, ABUNDANCE)],
                  all.x = TRUE, by = c("FEATURE", "RUN", "LABEL"))
            # TODO: log, whether any changes were made here
        } else {
            any_missing = as.character(unique(.getMissingRunsPerFeature(input)[, FEATURE]))
            msg = paste("The following features have missing values in at least one run.",
                        paste(any_missing, sep = ",\n ", collapse = ",\n "))
            getOption("MSstatsLog")("WARN", msg)
            getOption("MSstatsMsg")("WARN", msg)
        }
    }
    input
}

.getMissingRunsPerFeature = function(input) {
    n_runs = data.table::uniqueN(input$RUN)
    any_missing = input[, list(n_measurements = data.table::uniqueN(RUN)),
                        by = c("FRACTION", "LABEL", "FEATURE")]
    any_missing = any_missing[n_measurements < n_runs]
    any_missing
}


.checkDuplicatedMeasurements = function(input) {
    n_measurements = FEATURE = NULL
    counts = input[, list(n_measurements = .N),
                   by = c("FRACTION", "PROTEIN", "LABEL", "FEATURE", "RUN")]
    counts = unique(counts[n_measurements > 1, as.character(FEATURE)])
    if (length(counts) > 0) {
        msg = paste("The following features have duplicated measurements in some runs:",
                    "please remove the duplicates.", "\n",
                    paste(counts, sep = ",\n ", collapse = ",\n "))
        # TODO: report separately for L and H
        getOption("MSstatsLog")("WARN", msg)
        stop(msg)
    }
}

