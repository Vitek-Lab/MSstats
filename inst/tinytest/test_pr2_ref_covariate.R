# Tests for PR 2: Rename ref regression covariate → ref_covariate
#
# Pure symbol rename across utils_imputation.R, utils_summarization.R,
# dataProcess.R, and utils_summarization_prepare.R.  No logic changes.
#
# Tests verify:
#   1. .fitLinearModel uses the column named 'ref_covariate' in its formula
#      (coefficient names contain 'ref_covariate', not 'ref').
#   2. Supplying a column named 'ref' (the old name) instead of 'ref_covariate'
#      triggers an error because the formula references 'ref_covariate'.
#   3. Regression: dataProcess(SRMRawData) still produces valid ProteinLevelData.

# Helper -----------------------------------------------------------------------

# Build a single-protein data.table matching the structure that
# MSstatsSummarizeSingleLinear passes to .fitLinearModel in the labeled path.
# H rows carry ref_covariate = "0" (reference level); L rows carry the run id.
make_labeled_prot <- function(col_name = "ref_covariate") {
    dt <- data.table::data.table(
        PROTEIN       = rep("P1", 12),
        FEATURE       = factor(rep(c("F1", "F2", "F3"), each = 4)),
        RUN           = factor(rep(c("R1","R2","R3","R4"), 3)),
        LABEL         = rep(c("H","H","L","L"), 3),
        ABUNDANCE     = c(10, 11, 14, 15,
                          10.2, 11.2, 14.2, 15.2,
                          9.8, 10.8, 13.8, 14.8),
        newABUNDANCE  = c(10, 11, 14, 15,
                          10.2, 11.2, 14.2, 15.2,
                          9.8, 10.8, 13.8, 14.8),
        weights       = rep(NA_real_, 12)
    )
    # Add the covariate column under the requested name
    dt[[col_name]] <- factor(
        rep(c("0","0","R3","R4"), 3),
        levels = c("0","R3","R4")
    )
    dt
}

# --- Unit: .fitLinearModel uses 'ref_covariate' column name ------------------

lm_fit <- MSstats:::.fitLinearModel(
    make_labeled_prot("ref_covariate"),
    is_single_feature = FALSE,
    is_labeled        = TRUE,
    equal_variances   = TRUE
)

expect_true(
    any(grepl("ref_covariate", names(coef(lm_fit)))),
    info = "ref_covariate rename: model coefficients must include 'ref_covariate' terms"
)
expect_false(
    any(grepl("^ref[^_]", names(coef(lm_fit)))),
    info = "ref_covariate rename: coefficient names must NOT start with bare 'ref'"
)

# --- Providing 'ref' (old name) instead of 'ref_covariate' must error --------

expect_error(
    MSstats:::.fitLinearModel(
        make_labeled_prot("ref"),     # wrong column name
        is_single_feature = FALSE,
        is_labeled        = TRUE,
        equal_variances   = TRUE
    ),
    info = "ref_covariate rename: formula references 'ref_covariate'; 'ref' column must cause an error"
)

# --- Regression: SRMRawData still summarizes after the rename ----------------

quant_srm <- dataProcess(SRMRawData, use_log_file = FALSE)

expect_true(
    nrow(quant_srm$ProteinLevelData) > 0,
    info = "regression: SRMRawData must still produce ProteinLevelData after ref→ref_covariate rename"
)
expect_true(
    "LogIntensities" %in% colnames(quant_srm$ProteinLevelData),
    info = "regression: ProteinLevelData must have LogIntensities column"
)
expect_true(
    "GROUP" %in% colnames(quant_srm$ProteinLevelData),
    info = "regression: ProteinLevelData must have GROUP column"
)
