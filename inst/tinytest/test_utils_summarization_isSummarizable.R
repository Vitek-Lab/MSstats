# Coverage for .isSummarizable in utils_summarization.R: the guard that decides
# whether a protein can be summarized with TMP. Each early-return branch
# corresponds to a different degenerate protein, exercised here directly.
MSstatsConvert::MSstatsLogsSettings(FALSE)

base = data.table::data.table(
    PROTEIN = "P1", RUN = c("r1", "r2"),
    newABUNDANCE = c(10, 11),
    n_obs = c(2, 2), n_obs_run = c(2, 2),
    prop_features = c(1, 1)
)

# Test 1: a well-behaved protein is returned unchanged
res_ok = MSstats:::.isSummarizable(data.table::copy(base), TRUE)
expect_true(!is.null(res_ok))
expect_equal(nrow(res_ok), 2)

# Test 2: all measurements missing/zero -> NULL
b_na = data.table::copy(base); b_na$newABUNDANCE = NA_real_
expect_true(is.null(MSstats:::.isSummarizable(b_na, TRUE)))

# Test 3: all n_obs are zero -> NULL
b_nobs0 = data.table::copy(base); b_nobs0$n_obs = 0
expect_true(is.null(MSstats:::.isSummarizable(b_nobs0, TRUE)))

# Test 4: every feature has only a single measurement across runs -> NULL
b_nobs1 = data.table::copy(base); b_nobs1$n_obs = 1
expect_true(is.null(MSstats:::.isSummarizable(b_nobs1, TRUE)))

# Test 5: remove50missing = TRUE and all runs >50% missing -> NULL
b_50 = data.table::copy(base); b_50$prop_features = c(0.2, 0.3)
expect_true(is.null(MSstats:::.isSummarizable(b_50, TRUE)))

# Test 6: same data with remove50missing = FALSE is summarizable (branch skipped)
b_50b = data.table::copy(base); b_50b$prop_features = c(0.2, 0.3)
expect_true(!is.null(MSstats:::.isSummarizable(b_50b, FALSE)))

# Test 7: runs with no observations are filtered out, others retained
b_missing_run = data.table::data.table(
    PROTEIN = "P1", RUN = c("r1", "r2", "r3"),
    newABUNDANCE = c(10, 11, 12),
    n_obs = c(2, 2, 2), n_obs_run = c(2, 2, 0),
    prop_features = c(1, 1, 1)
)
res_filtered = MSstats:::.isSummarizable(b_missing_run, TRUE)
expect_equal(nrow(res_filtered), 2)
expect_false("r3" %in% res_filtered$RUN)
