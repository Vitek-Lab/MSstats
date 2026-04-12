make_srm_prep_input <- function() {
    data.table::data.table(
        PROTEIN   = rep("P1", 8),
        FEATURE   = rep(c("F1", "F2"), each = 4),
        RUN       = c("R1","R1","R2","R2","R1","R1","R2","R2"),
        LABEL     = c("H","L","H","L","H","L","H","L"),
        ABUNDANCE = c(10, 14, 11, 15, 10.2, 14.2, 11.2, 15.2)
    )
}

make_labelfree_prep_input <- function() {
    data.table::data.table(
        PROTEIN   = rep("P1", 8),
        FEATURE   = rep(c("F1", "F2"), each = 4),
        RUN       = rep(c("R1","R2","R3","R4"), 2),
        LABEL     = rep("L", 8),
        ABUNDANCE = c(10, 11, 12, 13, 10.2, 11.2, 12.2, 13.2)
    )
}

MSstatsConvert::MSstatsLogsSettings(FALSE)

prep_labeled <- MSstatsPrepareForSummarization(
    make_srm_prep_input(),
    method = "TMP", impute = FALSE,
    censored_symbol = NULL,
    remove_uninformative_feature_outlier = FALSE
)

expect_true(
    "ref_covariate" %in% colnames(prep_labeled),
    info = "ref_covariate: MSstatsPrepareForSummarization must create column for labeled data"
)
expect_true(
    is.factor(prep_labeled$ref_covariate),
    info = "ref_covariate: must be stored as a factor"
)
h_vals <- unique(as.character(prep_labeled$ref_covariate[prep_labeled$LABEL == "H"]))
expect_identical(
    h_vals, "0",
    info = "ref_covariate: heavy (H) rows must map to '0'"
)
l_rows <- prep_labeled[prep_labeled$LABEL == "L", ]
expect_true(
    all(as.character(l_rows$ref_covariate) == l_rows$RUN),
    info = "ref_covariate: light (L) rows must map to their RUN value"
)
l_runs    <- unique(prep_labeled$RUN[prep_labeled$LABEL == "L"])
rc_levels <- levels(prep_labeled$ref_covariate)
expect_true(
    "0" %in% rc_levels,
    info = "ref_covariate: '0' must be a factor level (placeholder for H rows)"
)
expect_true(
    all(l_runs %in% rc_levels),
    info = "ref_covariate: all light-run IDs must appear as factor levels"
)
expect_equal(
    sort(rc_levels),
    sort(c("0", l_runs)),
    info = "ref_covariate: factor levels must be exactly '0' plus the L-run IDs"
)
prep_labelfree <- MSstatsPrepareForSummarization(
    make_labelfree_prep_input(),
    method = "TMP", impute = FALSE,
    censored_symbol = NULL,
    remove_uninformative_feature_outlier = FALSE
)
expect_false(
    "ref_covariate" %in% colnames(prep_labelfree),
    info = "ref_covariate: must NOT be created for label-free (single LABEL) data"
)
