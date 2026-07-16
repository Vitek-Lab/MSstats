# Additional coverage for utils_groupcomparison_checks.R:
# .checkContrastMatrix and .checkGroupComparisonInput.

MSstatsConvert::MSstatsLogsSettings(FALSE)

# .checkGroupComparisonInput -------------------------------------------------

# Test 1: a data.table with all required columns passes through unchanged
good_input = data.table::data.table(
    RUN = 1, Protein = "P1", LogIntensities = 10, originalRUN = "r1",
    GROUP = "A", SUBJECT = 1, more50missing = FALSE, NumMeasuredFeature = 2
)
res = MSstats:::.checkGroupComparisonInput(good_input)
expect_true(data.table::is.data.table(res))
expect_equal(nrow(res), 1)

# Test 2: missing a required column (not processed by dataProcess) -> error
bad_input = data.table::data.table(
    RUN = 1, Protein = "P1", LogIntensities = 10
)
expect_error(MSstats:::.checkGroupComparisonInput(bad_input))

# .checkContrastMatrix -------------------------------------------------------

summarized = data.table::data.table(
    GROUP = factor(c("A", "B", "C"), levels = c("A", "B", "C"))
)

# Test 3: contrast matrix with a column per group passes
cm_ok = matrix(c(1, -1, 0), nrow = 1)
rownames(cm_ok) = "A-B"
colnames(cm_ok) = c("A", "B", "C")
res_cm = MSstats:::.checkContrastMatrix(cm_ok, summarized)
expect_equal(ncol(res_cm), 3)

# Test 4: wrong number of columns (fewer than number of groups) -> error
cm_wrong = matrix(c(1, -1), nrow = 1)
rownames(cm_wrong) = "A-B"
colnames(cm_wrong) = c("A", "B")
expect_error(MSstats:::.checkContrastMatrix(cm_wrong, summarized))
