# Coverage for .finalizeLinear in utils_output.R. This helper computes
# feature-level summary statistics for linear-model summarization; it is
# reachable only by direct call, so it is exercised here explicitly.
MSstatsConvert::MSstatsLogsSettings(FALSE)

make_input = function() {
    data.table::data.table(
        PROTEIN = factor("P1"),
        RUN = factor(rep(c("r1", "r2"), each = 2)),
        FEATURE = factor(rep(c("f1", "f2"), 2)),
        LABEL = "L",
        newABUNDANCE = c(10, 11, 10.5, NA),
        INTENSITY = c(1024, 2048, 1500, NA),
        censored = c(FALSE, FALSE, FALSE, TRUE),
        n_obs = c(2, 2, 2, 2),
        n_obs_run = c(2, 2, 2, 1),
        total_features = 2
    )
}

# Test 1: with a censored symbol, summary columns are added and NumImputedFeature
#         is always 0 for the linear method
res = MSstats:::.finalizeLinear(make_input(), "NA")
expect_true(all(c("NumMeasuredFeature", "MissingPercentage", "more50missing",
                  "nonmissing_orig", "NumImputedFeature") %in% colnames(res)))
expect_true(all(res$NumImputedFeature == 0))

# Test 2: the run with a censored/missing feature is flagged as >50% missing
run2 = res[RUN == "r2", ]
expect_true(all(run2$more50missing))
run1 = res[RUN == "r1", ]
expect_false(any(run1$more50missing))

# Test 3: censored_symbol = NULL takes the INTENSITY-based branch and still
#         produces the measured-feature statistics
res_null = MSstats:::.finalizeLinear(make_input(), NULL)
expect_true("NumMeasuredFeature" %in% colnames(res_null))
expect_equal(res_null[RUN == "r1", ]$NumMeasuredFeature, c(2, 2))
