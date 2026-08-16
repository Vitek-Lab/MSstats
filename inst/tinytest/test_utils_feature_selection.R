# Setup --------------------------------------------------------------------
MSstatsConvert::MSstatsLogsSettings(FALSE)
raw_input = MSstatsPrepareForDataProcess(DDARawData, 2, NULL)
raw_input = MSstatsNormalize(raw_input, "EQUALIZEMEDIANS")
raw_input = MSstatsMergeFractions(raw_input)
raw_input = MSstatsHandleMissing(raw_input, "TMP", TRUE, "NA", 0.999)

# Test 1: method = "all" keeps all features (no "remove" column added) -------
input_all = MSstatsSelectFeatures(data.table::copy(raw_input), "all")
expect_true(data.table::is.data.table(input_all))
expect_equal(nrow(input_all), nrow(raw_input))

# Test 2: method = "topN" flags features outside the top N per protein ------
input_topn = MSstatsSelectFeatures(data.table::copy(raw_input), "topN", top_n = 2)
expect_true("remove" %in% colnames(input_topn))
kept_per_protein = input_topn[!(remove),
                              list(n_features = length(unique(FEATURE))),
                              by = "PROTEIN"]
expect_true(all(kept_per_protein$n_features <= 2))

# Test 3: top3 is equivalent to topN with default top_n --------------------
input_top3 = MSstatsSelectFeatures(data.table::copy(raw_input), "top3")
expect_true("remove" %in% colnames(input_top3))

# Test 4: method = "highQuality" flags feature quality and outliers ---------
input_hq = MSstatsSelectFeatures(data.table::copy(raw_input), "highQuality",
                                  min_feature_count = 2)
expect_true(all(c("feature_quality", "is_outlier") %in% colnames(input_hq)))
expect_true(all(input_hq$feature_quality %in% c("Informative", "Noisy", "Uninformative")))
expect_true(is.logical(input_hq$is_outlier))
expect_equal(nrow(input_hq), nrow(raw_input))

# Test 5: invalid method errors ----------------------------------------------
expect_error(MSstatsSelectFeatures(data.table::copy(raw_input), "invalid_method"))

# Test 6: full pipeline still summarizes with highQuality selection ----------
selected = MSstatsSelectFeatures(data.table::copy(raw_input), "highQuality")
processed = getProcessed(selected)
prepared = MSstatsPrepareForSummarization(selected, "TMP", TRUE, "NA", FALSE)
summarized = MSstatsSummarizeWithSingleCore(prepared, "TMP", TRUE, "NA", FALSE, TRUE, 90)
output = MSstatsSummarizationOutput(prepared, summarized, processed, "TMP", TRUE, "NA")
expect_true(nrow(output$ProteinLevelData) > 0)

# Test 7: method = "highQuality-lsqr" flags feature quality and outliers ----
# (LSQR-based robust fit in place of "highQuality"'s MASS::rlm/QR-based one)
input_hq_lsqr = MSstatsSelectFeatures(data.table::copy(raw_input), "highQuality-lsqr",
                                      min_feature_count = 2)
expect_true(all(c("feature_quality", "is_outlier") %in% colnames(input_hq_lsqr)))
expect_true(all(input_hq_lsqr$feature_quality %in% c("Informative", "Noisy", "Uninformative")))
expect_true(is.logical(input_hq_lsqr$is_outlier))
expect_equal(nrow(input_hq_lsqr), nrow(raw_input))

# The two fitting methods should agree closely (not necessarily bit-identical,
# given different IRLS numerics) on the same real fixture.
expect_equal(mean(input_hq$feature_quality == input_hq_lsqr$feature_quality), 1)
expect_equal(mean(input_hq$is_outlier == input_hq_lsqr$is_outlier), 1)

# Test 8: full pipeline still summarizes with highQuality-lsqr selection ----
selected_lsqr = MSstatsSelectFeatures(data.table::copy(raw_input), "highQuality-lsqr")
processed_lsqr = getProcessed(selected_lsqr)
prepared_lsqr = MSstatsPrepareForSummarization(selected_lsqr, "TMP", TRUE, "NA", FALSE)
summarized_lsqr = MSstatsSummarizeWithSingleCore(prepared_lsqr, "TMP", TRUE, "NA", FALSE, TRUE, 90)
output_lsqr = MSstatsSummarizationOutput(prepared_lsqr, summarized_lsqr, processed_lsqr,
                                        "TMP", TRUE, "NA")
expect_true(nrow(output_lsqr$ProteinLevelData) > 0)
