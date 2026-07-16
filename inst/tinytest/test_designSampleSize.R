# Setup ------------------------------------------------------------------
QuantData = dataProcess(SRMRawData, use_log_file = FALSE)
comparison1 = matrix(c(-1,0,1,0,0,0,0,0,0,0), nrow=1)
comparison2 = matrix(c(-1,0,0,0,0,0,1,0,0,0), nrow=1)
comparison3 = matrix(c(-1,0,0,0,0,0,0,0,1,0), nrow=1)
comparison = rbind(comparison1, comparison2, comparison3)
row.names(comparison) = c("T3-T1", "T7-T1", "T9-T1")
groups = levels(QuantData$ProteinLevelData$GROUP)
colnames(comparison) = groups[order(as.numeric(groups))]

testResultMultiComparisons = groupComparison(contrast.matrix = comparison,
                                              data = QuantData,
                                              use_log_file = FALSE)

# Test 1: designSampleSize with numSample = TRUE (calculate sample size) -----
result_sample = designSampleSize(data = testResultMultiComparisons$FittedModel,
                                  numSample = TRUE,
                                  desiredFC = c(1.25, 1.75), FDR = 0.05, power = 0.8,
                                  use_log_file = FALSE)
expect_true(is.data.frame(result_sample))
expect_true(all(c("desiredFC", "numSample", "FDR", "power", "CV") %in%
                colnames(result_sample)))
expect_true(all(result_sample$numSample > 0))

# Test 2: designSampleSize with power = TRUE (calculate power) ---------------
result_power = designSampleSize(data = testResultMultiComparisons$FittedModel,
                                 numSample = 2,
                                 desiredFC = c(1.25, 1.75), FDR = 0.05, power = TRUE,
                                 use_log_file = FALSE)
expect_true(is.data.frame(result_power))
expect_true(all(c("desiredFC", "numSample", "FDR", "power", "CV") %in%
                colnames(result_power)))
expect_true(all(result_power$power >= 0 & result_power$power <= 1))
expect_true(all(result_power$numSample == 2))

# Test 3: increasing desired fold change reduces the required sample size ---
result_sample_small_fc = designSampleSize(data = testResultMultiComparisons$FittedModel,
                                           numSample = TRUE,
                                           desiredFC = c(1.05, 1.075), FDR = 0.05, power = 0.8,
                                           use_log_file = FALSE)
expect_true(max(result_sample_small_fc$numSample) >= min(result_sample$numSample))

# Test designSampleSizePlots (ggplot2 branch, isPlotly = FALSE) ---------------

tmp_dir = tempfile("msstats_designsamplesize_")
dir.create(tmp_dir)
pdf(file.path(tmp_dir, "plots.pdf"))

# Test 4: numSample plot (base graphics) does not error
expect_silent(designSampleSizePlots(data = result_sample))

# Test 5: power plot (base graphics) does not error
expect_silent(designSampleSizePlots(data = result_power))

dev.off()
unlink(tmp_dir, recursive = TRUE)

# Test designSampleSizePlots (plotly branch, isPlotly = TRUE) -----------------

# Test 6: numSample plot returns a plotly object
plotly_sample = designSampleSizePlots(data = result_sample, isPlotly = TRUE)
expect_true(inherits(plotly_sample, "plotly"))

# Test 7: power plot returns a plotly object
plotly_power = designSampleSizePlots(data = result_power, isPlotly = TRUE)
expect_true(inherits(plotly_power, "plotly"))

# Test .getVarComponent / .getMedianSigmaSubject / .calculatePower / .getNumSample

# Test 8: .getVarComponent extracts Error/Subject/GroupBySubject per protein
var_component = MSstats:::.getVarComponent(testResultMultiComparisons$FittedModel)
expect_true(all(c("Error", "Subject", "GroupBySubject", "Protein") %in%
                colnames(var_component)))
expect_equal(nrow(var_component), length(testResultMultiComparisons$FittedModel))

# Test 9: .getMedianSigmaSubject returns 0 when all components are NA --------
all_na_component = data.frame(Subject = NA_real_, GroupBySubject = NA_real_)
expect_equal(MSstats:::.getMedianSigmaSubject(all_na_component), 0)
