# Tests for PR 4: Per-label statistics and multi-label summarization
#
# Key changes:
#   - .prepareLinear / .prepareTMP: n_obs, n_obs_run, total_features,
#     prop_features now grouped by PROTEIN+FEATURE+LABEL (not just PROTEIN+FEATURE)
#     so H and L observations are counted independently.
#   - .runTukey / .fitTukey: when is_labeled_reference=FALSE (turnover), results
#     include a LABEL column so both H and L summaries are returned.
#   - .getNonMissingFilterStats: LABEL == "L" guard removed so filter is
#     applied to all rows, not just light rows.
#   - rbindlist(fill=TRUE) allows mixing H-only and L-only protein result lists.

# Helpers ----------------------------------------------------------------------

# Data with 2 features, 2 runs per label, both H and L labels
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

# --- .prepareLinear: n_obs computed per PROTEIN+FEATURE+LABEL ----------------

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

# --- .runTukey(is_labeled_reference=FALSE): result includes LABEL column -----

# Single-feature input — exercises the non-.fitTukey branch
tukey_sf <- data.table::data.table(
    PROTEIN      = rep("P1", 8),
    FEATURE      = factor(rep("F1", 8)),
    LABEL        = rep(c("L","H"), each = 4),
    RUN          = factor(rep(c("R1","R2","R3","R4"), 2)),
    newABUNDANCE = c(10, 11, 12, 13, 9, 10, 11, 12)
)

result_turnover <- MSstats:::.runTukey(
    tukey_sf, is_labeled_reference = FALSE,
    censored_symbol = NULL, remove50missing = FALSE
)

expect_true(
    "LABEL" %in% colnames(result_turnover),
    info = ".runTukey(is_labeled_reference=FALSE): result must have LABEL column"
)
expect_true(
    "L" %in% as.character(result_turnover$LABEL),
    info = ".runTukey(is_labeled_reference=FALSE): result must contain L rows"
)
expect_true(
    "H" %in% as.character(result_turnover$LABEL),
    info = ".runTukey(is_labeled_reference=FALSE): result must contain H rows"
)
expect_true(
    "LogIntensities" %in% colnames(result_turnover),
    info = ".runTukey(is_labeled_reference=FALSE): result must have LogIntensities column"
)

# SRM path (is_labeled_reference=TRUE) must NOT return LABEL column
result_srm <- MSstats:::.runTukey(
    tukey_sf, is_labeled_reference = TRUE,
    censored_symbol = NULL, remove50missing = FALSE
)
expect_false(
    "LABEL" %in% colnames(result_srm),
    info = ".runTukey(is_labeled_reference=TRUE): SRM path must not return LABEL column"
)

# --- .fitTukey(is_labeled_reference=FALSE): result includes LABEL column -----

# Multi-feature input — exercises .fitTukey
tukey_mf <- data.table::data.table(
    PROTEIN      = rep("P1", 16),
    FEATURE      = factor(rep(c("F1","F2"), each = 8)),
    LABEL        = rep(rep(c("L","L","L","L","H","H","H","H"), 1), 2),
    RUN          = factor(rep(c("R1","R2","R3","R4","R1","R2","R3","R4"), 2)),
    newABUNDANCE = c(10,11,12,13, 8,9,10,11, 10.2,11.2,12.2,13.2, 8.2,9.2,10.2,11.2)
)

result_tukey_turnover <- MSstats:::.fitTukey(tukey_mf, is_labeled_reference = FALSE)

expect_true(
    "LABEL" %in% colnames(result_tukey_turnover),
    info = ".fitTukey(is_labeled_reference=FALSE): result must include LABEL column"
)
expect_true(
    "LogIntensities" %in% colnames(result_tukey_turnover),
    info = ".fitTukey(is_labeled_reference=FALSE): result must have LogIntensities column"
)

result_tukey_srm <- MSstats:::.fitTukey(tukey_mf, is_labeled_reference = TRUE)
expect_false(
    "LABEL" %in% colnames(result_tukey_srm),
    info = ".fitTukey(is_labeled_reference=TRUE): SRM path must not return LABEL column"
)

# --- .getNonMissingFilterStats: no LABEL == "L" guard -----------------------

nonmiss_input <- data.table::data.table(
    LABEL        = c("L","L","H","H"),
    newABUNDANCE = c(10, NA, 8, NA),
    INTENSITY    = c(100L, NA_integer_, 80L, NA_integer_),
    censored     = c(FALSE, TRUE, FALSE, TRUE),
    n_obs_run    = c(2L, 2L, 2L, 2L),
    n_obs        = c(4L, 4L, 4L, 4L)
)

filter_result <- MSstats:::.getNonMissingFilterStats(nonmiss_input, censored_symbol = "NA")

# Both L and H non-missing rows should be flagged TRUE
expect_true(
    filter_result[nonmiss_input$LABEL == "L" & !is.na(nonmiss_input$newABUNDANCE)],
    info = ".getNonMissingFilterStats: L non-missing row must be TRUE"
)
expect_true(
    filter_result[nonmiss_input$LABEL == "H" & !is.na(nonmiss_input$newABUNDANCE)],
    info = ".getNonMissingFilterStats: H non-missing row must also be TRUE (no LABEL guard)"
)

# --- Regression: SRMRawData summarization unchanged --------------------------

quant_srm <- dataProcess(SRMRawData, use_log_file = FALSE)

expect_true(
    nrow(quant_srm$ProteinLevelData) > 0,
    info = "regression: SRMRawData must still produce ProteinLevelData after per-label changes"
)
expect_true(
    "LogIntensities" %in% colnames(quant_srm$ProteinLevelData),
    info = "regression: ProteinLevelData must still have LogIntensities"
)
