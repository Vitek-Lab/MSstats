make_srm_prep_input <- function() {
    data.table::data.table(
        PROTEIN   = rep("P1", 8),
        FEATURE   = rep(c("F1", "F2"), each = 4),
        RUN       = c("R1","R1","R2","R2","R1","R1","R2","R2"),
        LABEL     = c("H","L","H","L","H","L","H","L"),
        ABUNDANCE = c(10, 14, 11, 15, 10.2, 14.2, 11.2, 15.2),
        is_labeled_ref = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
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

# --- .prepareLinear: n_obs computed per PROTEIN+FEATURE+LABEL ----------------

make_two_label_input <- function() {
    data.table::data.table(
        PROTEIN  = rep("P1", 8),
        FEATURE  = factor(rep(c("F1", "F2"), each = 4)),
        LABEL    = rep(c("L","L","H","H"), 2),
        RUN      = factor(rep(c("R1","R2","R1","R2"), 2)),
        ABUNDANCE = c(10, 11, 8, 9,   10.5, 11.5, 8.5, 9.5),
        INTENSITY = rep(100L, 8),
        ANOMALYSCORES = rep(NA_real_, 8)
    )
}

prep_input <- make_two_label_input()
result_prep <- MSstats:::.prepareLinear(prep_input, impute = FALSE, censored_symbol = NULL)

# Each FEATURE+LABEL combination has 2 runs → n_obs should be 2
expect_equal(
    unique(result_prep[LABEL == "L" & FEATURE == "F1", n_obs]),
    2L,
    info = ".prepareLinear: n_obs must be 2 for L rows of F1 (2 L runs)"
)
expect_equal(
    unique(result_prep[LABEL == "H" & FEATURE == "F1", n_obs]),
    2L,
    info = ".prepareLinear: n_obs must be 2 for H rows of F1 (2 H runs)"
)
# Without LABEL grouping the old code would have returned 4 (all runs pooled)
expect_false(
    any(result_prep$n_obs == 4L),
    info = ".prepareLinear: n_obs must not be 4 (old un-grouped value)"
)

# total_features per PROTEIN+LABEL: 2 features per label
expect_equal(
    unique(result_prep[LABEL == "L", total_features]),
    2L,
    info = ".prepareLinear: total_features for L must be 2"
)
expect_equal(
    unique(result_prep[LABEL == "H", total_features]),
    2L,
    info = ".prepareLinear: total_features for H must be 2"
)

# --- .prepareLinear: is_labeled_reference=TRUE groups WITHOUT LABEL -----------
# For SRM, H is the normalization reference, not an independent label.
# n_obs must combine L and H observations for each feature (no per-label split).

result_prep_srm <- MSstats:::.prepareLinear(
    make_two_label_input(),
    impute = FALSE, censored_symbol = NULL,
    is_labeled_reference = TRUE
)

# 2 L runs + 2 H runs per feature → combined n_obs = 4
expect_equal(
    unique(result_prep_srm[LABEL == "L" & FEATURE == "F1", n_obs]),
    4L,
    info = ".prepareLinear(is_labeled_reference=TRUE): n_obs must combine L and H observations"
)
expect_equal(
    unique(result_prep_srm[LABEL == "H" & FEATURE == "F1", n_obs]),
    4L,
    info = ".prepareLinear(is_labeled_reference=TRUE): H rows must share the same n_obs as L rows"
)

# --- .prepareTMP: is_labeled_reference=FALSE groups by LABEL -----------------

result_tmp_unlabeled <- MSstats:::.prepareTMP(
    make_two_label_input(),
    impute = FALSE, censored_symbol = NULL,
    is_labeled_reference = FALSE
)

# Each FEATURE+LABEL combination has 2 runs → n_obs per label = 2
expect_equal(
    unique(result_tmp_unlabeled[LABEL == "L" & FEATURE == "F1", n_obs]),
    2L,
    info = ".prepareTMP(is_labeled_reference=FALSE): n_obs must be per-label (2 L runs)"
)
expect_equal(
    unique(result_tmp_unlabeled[LABEL == "H" & FEATURE == "F1", n_obs]),
    2L,
    info = ".prepareTMP(is_labeled_reference=FALSE): H n_obs must be counted independently"
)
expect_false(
    any(result_tmp_unlabeled$n_obs == 4L),
    info = ".prepareTMP(is_labeled_reference=FALSE): n_obs must not be 4 (pooled across labels)"
)

# total_features per PROTEIN+LABEL
expect_equal(
    unique(result_tmp_unlabeled[LABEL == "L", total_features]),
    2L,
    info = ".prepareTMP(is_labeled_reference=FALSE): total_features for L must be 2"
)
expect_equal(
    unique(result_tmp_unlabeled[LABEL == "H", total_features]),
    2L,
    info = ".prepareTMP(is_labeled_reference=FALSE): total_features for H must be 2"
)

# --- .prepareTMP: is_labeled_reference=TRUE groups WITHOUT LABEL -------------
# For SRM, .getNonMissingFilter marks H rows as non-informative (is_labeled_ref=TRUE
# → use_for_analysis=FALSE). Without LABEL grouping, each feature's n_obs is the
# count of nonmissing L rows, but that count is assigned to ALL rows of the feature
# (L and H alike). This prevents H rows from getting n_obs=0 and being filtered
# out before they can serve as the normalization reference in .adjustLRuns.

result_tmp_srm <- MSstats:::.prepareTMP(
    make_srm_prep_input(),
    impute = FALSE, censored_symbol = NULL,
    is_labeled_reference = TRUE
)

# H rows share the L-derived n_obs for their feature (= 2); they must not get 0
expect_equal(
    unique(result_tmp_srm[LABEL == "H", n_obs]),
    2L,
    info = ".prepareTMP(is_labeled_reference=TRUE): H rows must share n_obs with L rows (not 0)"
)
# n_obs_run: similarly, H rows must have n_obs_run > 0 so they survive .isSummarizable
expect_true(
    all(result_tmp_srm[LABEL == "H", n_obs_run] > 0),
    info = ".prepareTMP(is_labeled_reference=TRUE): H rows must have n_obs_run > 0"
)

