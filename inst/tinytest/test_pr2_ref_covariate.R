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


expect_true(
    "LogIntensities" %in% colnames(quant_srm$ProteinLevelData),
    info = "regression: ProteinLevelData must have LogIntensities column"
)
expect_true(
    "GROUP" %in% colnames(quant_srm$ProteinLevelData),
    info = "regression: ProteinLevelData must have GROUP column"
)
