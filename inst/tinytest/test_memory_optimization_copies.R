# Tests for in-place data.table operations in the dataProcess pipeline:
# column drops, keyed-lookup joins, targeted [i, j := v] writes, and
# the predicted_survival extraction in MSstatsSummarizationOutput.


# --- .normalizeMedian: column drops are in place ------------------------------
#
# .normalizeMedian() adds two temp columns (ABUNDANCE_RUN, ABUNDANCE_FRACTION)
# via :=, uses them, then removes them with data.table::set(j=, value=NULL)
# so the data.table is modified in-place. We verify the address is preserved.

MSstatsConvert::MSstatsLogsSettings(FALSE)

# Test 2a: .normalizeMedian preserves data.table identity ---
#
# data.table::address() returns the hex pointer of the data.table object.
# An in-place modification preserves the address; a copy changes it.

norm_input = data.table::data.table(
    PROTEIN = factor(rep(c("P1", "P2"), each = 20)),
    PEPTIDE = factor(rep(c("pep1", "pep2"), each = 10, times = 2)),
    FEATURE = factor(rep(c("feat1", "feat2"), each = 10, times = 2)),
    LABEL = rep("L", 40),
    RUN = factor(rep(1:10, 4)),
    FRACTION = rep(1L, 40),
    ABUNDANCE = rnorm(40, mean = 20, sd = 2),
    INTENSITY = 2^rnorm(40, mean = 20, sd = 2),
    GROUP_ORIGINAL = rep(c("A", "B"), each = 5, times = 4),
    SUBJECT_ORIGINAL = rep(1:10, 4),
    TRANSITION = factor(rep("y3_1", 40))
)

addr_before = data.table::address(norm_input)
norm_result = MSstats:::.normalizeMedian(norm_input)
addr_after = data.table::address(norm_result)

# The address should be the same — in-place modification, no copy
expect_equal(addr_before, addr_after,
             info = paste(".normalizeMedian should modify in-place (same address).",
                          "Before:", addr_before, "After:", addr_after))

# Temp columns should not exist in the result
expect_false("ABUNDANCE_RUN" %in% colnames(norm_result),
             info = "Temp column ABUNDANCE_RUN should be removed")
expect_false("ABUNDANCE_FRACTION" %in% colnames(norm_result),
             info = "Temp column ABUNDANCE_FRACTION should be removed")

# ABUNDANCE should still exist and be modified (normalized)
expect_true("ABUNDANCE" %in% colnames(norm_result),
            info = "ABUNDANCE column should still exist after normalization")


# Test 2b: .normalizeMedian does not allocate a full-table copy ---
#
# With in-place set(), the only allocations should be the two temp column
# vectors (~2x nrow doubles), not the whole table.

n_rows = 100000
big_input = data.table::data.table(
    PROTEIN = factor(rep(paste0("P", 1:100), each = n_rows / 100)),
    PEPTIDE = factor(rep(paste0("pep", 1:200), each = n_rows / 200)),
    FEATURE = factor(rep(paste0("feat", 1:200), each = n_rows / 200)),
    LABEL = rep("L", n_rows),
    RUN = factor(rep(1:50, n_rows / 50)),
    FRACTION = rep(1L, n_rows),
    ABUNDANCE = rnorm(n_rows, mean = 20, sd = 2),
    INTENSITY = 2^rnorm(n_rows, mean = 20, sd = 2),
    GROUP_ORIGINAL = rep(c("A", "B"), each = 25, times = n_rows / 50),
    SUBJECT_ORIGINAL = rep(1:50, n_rows / 50),
    TRANSITION = factor(rep("y3_1", n_rows))
)

# Measure the size of the input table and one column for reference.
# We use these to reason about whether the allocation during
# .normalizeMedian is "2 temp columns worth" vs "full table copy worth".
table_size = as.numeric(object.size(big_input))
one_col_size = as.numeric(object.size(big_input$ABUNDANCE))
ncols = ncol(big_input)

# --- Memory measurement ---
#
# R manages data in "Vcells" (vector cells) — this is where column data
# lives. gc() returns a matrix with memory stats; row 2 = Vcells,
# column 6 = "max used (Mb)" = peak Vcells since last reset.
#
# gc(reset = TRUE) resets the "max used" counter to the current level,
# giving us a clean window to measure just the .normalizeMedian call.
gc(reset = TRUE)
big_result = MSstats:::.normalizeMedian(big_input)
gc_after = gc()

# Extract peak vector memory (Vcells, max used, in MB).
# We store this for informational purposes but do NOT assert on it
# because gc() accounting is non-deterministic:
#   - R may GC the temp columns before we call gc(), or may not
#   - Other small allocations (formula parsing, string interning)
#     contaminate the measurement
#   - Two identical runs can give different numbers
# The address() check below is the reliable, deterministic test.
peak_vcells_mb = gc_after[2, 6]

# Sanity check on the test data itself (not on the function):
# Verify that 2 temp columns are much smaller than the full table.
# If they were similar sizes, we couldn't distinguish "allocated 2 temp
# columns" from "copied the whole table" and the test would be meaningless.
# For our 11-column, 100K-row table: 2 cols ≈ 1.6 MB, full table ≈ 6 MB.
expect_true(one_col_size * 2 < table_size * 0.5,
            info = paste("Sanity check: 2 columns should be much less than",
                         "50% of the table. 2 cols:", one_col_size * 2,
                         "50% table:", table_size * 0.5))

# The real test: same address = same R object = no copy was made.
# This is deterministic and reliable, unlike gc() measurements.
expect_equal(data.table::address(big_input), data.table::address(big_result),
             info = "Large table: .normalizeMedian should modify in-place")
expect_false("ABUNDANCE_RUN" %in% colnames(big_result),
             info = "Large table: temp column ABUNDANCE_RUN should be removed")


# --- .finalizeTMP: in-place update from predicted_survival --------------------
#
# .finalizeTMP replaces the pre-imputation newABUNDANCE column in the full
# feature-level table with imputed values from predicted_survival, using:
#   set(input, j = "newABUNDANCE", value = NULL)  — remove old column
#   input[, c(...) := .(NA_real_, NA_real_)]       — add NA placeholders
#   input[predicted_survival, ... := ..., on = .]  — update matched rows
#
# These tests verify:
# (a) The data.table is modified in-place (address preserved)
# (b) Matched rows get imputed values from predicted_survival
# (c) Unmatched rows are left as NA (equivalent to all.x = TRUE)
# (d) Allocation stays small relative to a full-table copy

# Test 3a: .finalizeTMP modifies input in-place with correct values ---

# Build a feature-level input table mimicking post-summarization state.
# Two proteins, 2 features each, 3 runs = 12 rows per label.
n = 12
finalize_input = data.table::data.table(
    PROTEIN = factor(rep(c("P1", "P2"), each = 6)),
    PEPTIDE = factor(rep(c("pep1", "pep2"), each = 3, times = 2)),
    FEATURE = factor(rep(c("feat1", "feat2"), each = 3, times = 2)),
    LABEL = rep("L", n),
    GROUP = factor(rep(c("A", "B", "C"), times = 4)),
    RUN = factor(rep(1:3, times = 4)),
    SUBJECT = factor(rep(1:3, times = 4)),
    FRACTION = rep(1L, n),
    INTENSITY = 2^rnorm(n, 20, 2),
    ABUNDANCE = rnorm(n, 20, 2),
    # Pre-imputation newABUNDANCE: some are NA (censored)
    newABUNDANCE = c(10.1, NA, 10.5, 11.2, 10.8, NA,
                     9.5, NA, 10.2, 11.0, NA, 10.6),
    cen = c(1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1),
    censored = c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE,
                 FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
    originalRUN = factor(rep(paste0("Run", 1:3), times = 4)),
    n_obs = rep(3, n),
    n_obs_run = rep(2, n),
    total_features = rep(2, n),
    nonmissing = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE,
                   TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
)

# Build predicted_survival: simulates what rbindlist of per-protein survival
# tables looks like. Only contains rows that passed the n_obs/n_obs_run filter
# in MSstatsSummarizeSingleTMP — so fewer rows than input.
# Rows 2, 6, 8, 11 are censored — only some of them appear in predicted_survival
# (e.g., rows filtered out by n_obs <= 1 would be absent).
# Here we include most rows but deliberately EXCLUDE rows for RUN=2, FEATURE=feat2
# of P2 to test the unmatched (NA) case.
predicted_survival = data.table::data.table(
    newABUNDANCE = c(10.1, 8.5, 10.5, 11.2, 10.8, 7.9,
                     9.5, 8.2, 10.2, 11.0, 10.6),
    cen = c(1, 0, 1, 1, 1, 0,
            1, 0, 1, 1, 1),
    RUN = factor(c(1, 2, 3, 1, 2, 3,
                   1, 2, 3, 1, 3)),
    FEATURE = factor(c(rep("feat1", 3), rep("feat2", 3),
                       rep("feat1", 3), rep("feat2", 2))),
    predicted = c(NA, 8.5, NA, NA, NA, 7.9,
                  NA, 8.2, NA, NA, NA)
)

addr_before = data.table::address(finalize_input)

# Call .finalizeTMP with impute = TRUE.
# This should replace newABUNDANCE and add predicted via in-place join.
result = MSstats:::.finalizeTMP(finalize_input, censored_symbol = "NA",
                                impute = TRUE,
                                predicted_survival = predicted_survival)

addr_after = data.table::address(result)

# (a) Address preserved — the data.table was modified in-place.
expect_equal(addr_before, addr_after,
             info = paste(".finalizeTMP should modify input in-place.",
                          "Before:", addr_before, "After:", addr_after))

# (b) Matched rows should have newABUNDANCE from predicted_survival.
# Any row whose (cen, RUN, FEATURE) key exists in predicted_survival
# should get a non-NA newABUNDANCE value.
matched_count = sum(!is.na(result$newABUNDANCE))
expect_true(matched_count > 0,
            info = paste("Matched rows should have non-NA newABUNDANCE.",
                         "Found", matched_count, "non-NA values"))

# (c) Unmatched rows should have NA.
# predicted_survival has 11 rows but input has 12. The missing key
# (cen=0, RUN=2, FEATURE=feat2) has no match, so that row should be NA.
# We check by looking for the key that's absent in predicted_survival.
unmatched = result[cen == 0 & RUN == 2 & FEATURE == "feat2", newABUNDANCE]
if (length(unmatched) > 0) {
    expect_true(is.na(unmatched[1]),
                info = "Unmatched row should get NA for newABUNDANCE")
}

# (d) The predicted column should exist.
expect_true("predicted" %in% colnames(result),
            info = ".finalizeTMP should add the 'predicted' column")


# Test 3b: .finalizeTMP at scale — memory stays low ---
#
# At ~50K rows, .finalizeTMP should modify the input in place; peak Vcells
# should stay close to the input size, not balloon to a multi-table merge.

n_big = 50000
n_runs = 50
n_features = 20
n_proteins = n_big / (n_runs * n_features)  # = 50

big_finalize_input = data.table::data.table(
    PROTEIN = factor(rep(paste0("P", 1:n_proteins), each = n_runs * n_features)),
    FEATURE = factor(rep(paste0("feat", 1:n_features), each = n_runs, times = n_proteins)),
    LABEL = rep("L", n_big),
    RUN = factor(rep(1:n_runs, times = n_proteins * n_features)),
    INTENSITY = 2^rnorm(n_big, 20, 2),
    ABUNDANCE = rnorm(n_big, 20, 2),
    newABUNDANCE = rnorm(n_big, 20, 2),
    cen = sample(c(0, 1), n_big, replace = TRUE, prob = c(0.2, 0.8)),
    censored = FALSE,
    originalRUN = factor(rep(paste0("Run", 1:n_runs), times = n_proteins * n_features)),
    GROUP = factor(rep("A", n_big)),
    SUBJECT = factor(rep(1:n_runs, times = n_proteins * n_features)),
    n_obs = rep(n_runs, n_big),
    n_obs_run = rep(n_features, n_big),
    total_features = rep(n_features, n_big),
    nonmissing = rep(TRUE, n_big)
)
big_finalize_input[cen == 0, censored := TRUE]
big_finalize_input[cen == 0, newABUNDANCE := NA]

# predicted_survival: exclude specific (RUN, FEATURE) combinations to
# guarantee some rows in input have no match. Random row sampling doesn't
# work because multiple rows can share the same join key, so excluded rows
# may still match via a duplicate key in the kept rows.
excluded_features = c("feat1", "feat2")
excluded_run = "1"
big_predicted = big_finalize_input[
    !(FEATURE %in% excluded_features & RUN == excluded_run),
    .(newABUNDANCE = ifelse(cen == 0, 7.0, newABUNDANCE),
      cen, RUN, FEATURE,
      predicted = ifelse(cen == 0, 7.0, NA_real_))]

input_size = as.numeric(object.size(big_finalize_input))
addr_big_before = data.table::address(big_finalize_input)

# Reset GC, run .finalizeTMP, measure peak.
gc(reset = TRUE)
big_finalize_result = MSstats:::.finalizeTMP(
    big_finalize_input, censored_symbol = "NA", impute = TRUE,
    predicted_survival = big_predicted)
gc_after_finalize = gc()

# Peak Vcells (MB) — informational.
peak_mb_finalize = gc_after_finalize[2, 6]

# Address preserved = in-place modification, no full-table copy.
expect_equal(addr_big_before, data.table::address(big_finalize_result),
             info = paste("Large .finalizeTMP: should modify in-place.",
                          "Input size:", round(input_size / 1e6, 1), "MB.",
                          "Peak Vcells:", peak_mb_finalize, "MB"))

# newABUNDANCE and predicted columns should both exist.
expect_true("newABUNDANCE" %in% colnames(big_finalize_result),
            info = "Large .finalizeTMP: newABUNDANCE should exist after join")
expect_true("predicted" %in% colnames(big_finalize_result),
            info = "Large .finalizeTMP: predicted should exist after join")

# Rows with excluded (RUN, FEATURE) keys should have NA for newABUNDANCE,
# because they have no match in predicted_survival.
unmatched_mask = big_finalize_result$FEATURE %in% excluded_features &
                 big_finalize_result$RUN == excluded_run
na_count = sum(is.na(big_finalize_result$newABUNDANCE[unmatched_mask]))
expect_true(na_count > 0,
            info = paste("Rows with excluded keys should have NA newABUNDANCE.",
                         "Found", na_count, "NAs out of",
                         sum(unmatched_mask), "unmatched rows"))


# --- MSstatsSummarizationOutput / .finalizeTMP: predicted_survival contract ---
#
# MSstatsSummarizationOutput splits the nested 'summarized' list upfront:
#   predicted_survival = rbindlist(lapply(summarized, function(x) x[[2]]))
#   protein_summaries  = lapply(summarized, function(x) x[[1]])
#   rm(summarized)
# .finalizeTMP receives predicted_survival directly (a combined data.table),
# not the full nested list.
#
# These tests verify:
# (a) .finalizeTMP accepts predicted_survival as a data.table
# (b) MSstatsSummarizationOutput correctly decomposes the nested list
# (c) The nested list does not persist through .finalizeTMP

# Test 4a: .finalizeTMP takes a combined data.table, not a nested list ---
#
# .finalizeTMP's 4th parameter is 'predicted_survival', a data.table.

# Rebuild a finalize_input fixture since .finalizeTMP modified the previous
# one in-place above.
finalize_input_4 = data.table::data.table(
    PROTEIN = factor(rep(c("P1", "P2"), each = 6)),
    FEATURE = factor(rep(c("feat1", "feat2"), each = 3, times = 2)),
    LABEL = rep("L", 12),
    RUN = factor(rep(1:3, times = 4)),
    INTENSITY = 2^rnorm(12, 20, 2),
    ABUNDANCE = rnorm(12, 20, 2),
    newABUNDANCE = c(10.1, NA, 10.5, 11.2, 10.8, NA,
                     9.5, NA, 10.2, 11.0, NA, 10.6),
    cen = c(1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1),
    censored = c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE,
                 FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
    originalRUN = factor(rep(paste0("Run", 1:3), times = 4)),
    GROUP = factor(rep(c("A", "B", "C"), times = 4)),
    SUBJECT = factor(rep(1:3, times = 4)),
    n_obs = rep(3, 12),
    n_obs_run = rep(2, 12),
    total_features = rep(2, 12),
    nonmissing = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE,
                   TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
)

# predicted_survival as a plain data.table (the new interface).
pred_surv_4 = data.table::data.table(
    newABUNDANCE = c(10.1, 8.5, 10.5, 11.2, 10.8, 7.9,
                     9.5, 8.2, 10.2, 11.0, 10.6, 8.0),
    cen = c(1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0),
    RUN = factor(rep(1:3, times = 4)),
    FEATURE = factor(rep(c("feat1", "feat2"), each = 3, times = 2)),
    predicted = c(NA, 8.5, NA, NA, NA, 7.9, NA, 8.2, NA, NA, NA, 8.0)
)

# This should work — .finalizeTMP now expects a data.table, not a list.
result_4 = MSstats:::.finalizeTMP(finalize_input_4, censored_symbol = "NA",
                                   impute = TRUE,
                                   predicted_survival = pred_surv_4)
expect_true(data.table::is.data.table(result_4) || is.data.frame(result_4),
            info = ".finalizeTMP should return a data.frame/data.table when given predicted_survival directly")
expect_true("newABUNDANCE" %in% colnames(result_4),
            info = ".finalizeTMP should produce newABUNDANCE from a data.table input")
expect_true("predicted" %in% colnames(result_4),
            info = ".finalizeTMP should produce predicted column from a data.table input")


# Test 4b: MSstatsSummarizationOutput correctly splits the nested list ---
#
# Build a minimal 'summarized' in the nested format that the summarization
# loop produces: list(list(protein_result, survival), ...).
# Verify that MSstatsSummarizationOutput decomposes it correctly and
# produces valid output.

# Use the real pipeline on a small built-in dataset to get a realistic
# summarized list, then verify the output structure.
MSstatsConvert::MSstatsLogsSettings(FALSE)
small_input = MSstatsPrepareForDataProcess(SRMRawData, 2, NULL)
small_input = MSstatsNormalize(small_input, "EQUALIZEMEDIANS")
small_input = MSstatsMergeFractions(small_input)
small_input = MSstatsHandleMissing(small_input, "TMP", TRUE, "NA", 0.999)
small_input = MSstatsSelectFeatures(small_input, "all")
processed = getProcessed(small_input)
small_input = MSstatsPrepareForSummarization(small_input, "TMP", TRUE, "NA", FALSE)
summarized_list = MSstatsSummarizeWithSingleCore(small_input, "TMP", TRUE, "NA",
                                                  FALSE, TRUE, 90)

# Verify summarized_list is a nested list (the format we're testing).
expect_true(is.list(summarized_list),
            info = "Summarized should be a list")
expect_true(is.list(summarized_list[[1]]),
            info = "Each element should be a list(protein_result, survival)")
expect_equal(length(summarized_list[[1]]), 2,
             info = "Each element should have exactly 2 components")

# MSstatsSummarizationOutput should decompose the nested list and produce
# valid FeatureLevelData / ProteinLevelData.
output = MSstatsSummarizationOutput(small_input, summarized_list, processed,
                                     "TMP", TRUE, "NA")

expect_true(!is.null(output$FeatureLevelData),
            info = "Output should have FeatureLevelData")
expect_true(!is.null(output$ProteinLevelData),
            info = "Output should have ProteinLevelData")
expect_true(nrow(output$FeatureLevelData) > 0,
            info = "FeatureLevelData should have rows")
expect_true(nrow(output$ProteinLevelData) > 0,
            info = "ProteinLevelData should have rows")
expect_true("newABUNDANCE" %in% colnames(output$FeatureLevelData),
            info = "FeatureLevelData should contain newABUNDANCE (from imputation)")
expect_true("predicted" %in% colnames(output$FeatureLevelData),
            info = "FeatureLevelData should contain predicted column")


# Test 4c: predicted_survival alone is smaller than the full nested list ---
#
# .finalizeTMP receives only predicted_survival (the survival half of the
# nested summarized list), not the full list. We verify the combined
# predicted_survival table is smaller than the full nested list.

# Measure sizes of the two halves from the real summarized list.
survival_tables = lapply(summarized_list, function(x) x[[2]])
protein_tables = lapply(summarized_list, function(x) x[[1]])

size_full_list = as.numeric(object.size(summarized_list))
size_survival_only = as.numeric(object.size(survival_tables))
size_protein_only = as.numeric(object.size(protein_tables))
size_combined_survival = as.numeric(
    object.size(data.table::rbindlist(survival_tables, fill = TRUE)))

# The combined predicted_survival table (what .finalizeTMP receives now)
# should be smaller than the full nested list (what it received before).
expect_true(size_combined_survival < size_full_list,
            info = paste("Combined predicted_survival should be smaller than",
                         "the full nested list.",
                         "predicted_survival:", size_combined_survival, "bytes.",
                         "Full nested list:", size_full_list, "bytes."))
