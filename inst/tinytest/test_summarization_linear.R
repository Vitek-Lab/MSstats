# Coverage for the linear summarization path in utils_summarization.R:
# .fitLinearModel (both single- and multi-feature) and .updateUnequalVariances
# (equalFeatureVar = FALSE). The default dataProcess uses TMP, so none of the
# linear-model code was previously exercised by the test suite.

# Test 1: linear summarization with equal feature variance ------------------
q_equal = dataProcess(DDARawData, summaryMethod = "linear",
                      equalFeatureVar = TRUE, use_log_file = FALSE)
expect_true(is.list(q_equal))
expect_true(nrow(q_equal$ProteinLevelData) > 0)
expect_true("LogIntensities" %in% colnames(q_equal$ProteinLevelData))

# Test 2: linear summarization with unequal feature variance ----------------
#         (this routes through .updateUnequalVariances / loess reweighting)
q_unequal = dataProcess(DDARawData, summaryMethod = "linear",
                        equalFeatureVar = FALSE, use_log_file = FALSE)
expect_true(nrow(q_unequal$ProteinLevelData) > 0)
expect_true(all(is.finite(q_unequal$ProteinLevelData$LogIntensities)))

# Test 3: both variance options summarize the same set of proteins ----------
expect_equal(sort(unique(as.character(q_equal$ProteinLevelData$Protein))),
             sort(unique(as.character(q_unequal$ProteinLevelData$Protein))))

# Test 4: .fitLinearModel directly, single feature vs multi-feature ---------
# Build a minimal single-protein feature-level table with the columns the
# linear model fitter expects.
make_feature_data = function(n_features) {
    runs = paste0("R", 1:4)
    features = paste0("F", seq_len(n_features))
    grid = expand.grid(RUN = runs, FEATURE = features,
                       stringsAsFactors = FALSE)
    data.table::data.table(
        FEATURE = factor(grid$FEATURE),
        RUN = factor(grid$RUN),
        newABUNDANCE = 20 + rnorm(nrow(grid)),
        weights = NA_real_
    )
}

set.seed(1)
single = make_feature_data(1)
fit_single = MSstats:::.fitLinearModel(single, is_single_feature = TRUE,
                                       is_labeled = FALSE, equal_variances = TRUE)
expect_true(inherits(fit_single, "lm"))

multi = make_feature_data(3)
fit_multi = MSstats:::.fitLinearModel(multi, is_single_feature = FALSE,
                                      is_labeled = FALSE, equal_variances = TRUE)
expect_true(inherits(fit_multi, "lm"))

# Test 5: .fitLinearModel with equal_variances = FALSE reweights via
#         .updateUnequalVariances and still returns a fitted lm
fit_unequal = MSstats:::.fitLinearModel(multi, is_single_feature = FALSE,
                                        is_labeled = FALSE, equal_variances = FALSE)
expect_true(inherits(fit_unequal, "lm"))
