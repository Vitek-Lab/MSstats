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
comparison_result = testResultMultiComparisons$ComparisonResult

# Test 1: invalid type errors -------------------------------------------------
expect_error(groupComparisonPlots(data = comparison_result, type = "INVALID"))

# Test 2: invalid which.Comparison label errors -------------------------------
expect_error(
    groupComparisonPlots(data = comparison_result, type = "VolcanoPlot",
                          which.Comparison = "NotARealComparison")
)

tmp_dir = tempfile("msstats_groupcomparisonplots_")
dir.create(tmp_dir)
address_prefix = paste0(tmp_dir, "/")

# Test 3: Heatmap (ggplot2) creates a pdf file and always warns --------------
expect_warning(
    groupComparisonPlots(data = comparison_result, type = "Heatmap",
                          logBase.pvalue = 2, address = address_prefix)
)
expect_true(file.exists(paste0(address_prefix, "Heatmap.pdf")))

# Test 4: Heatmap with FCcutoff also runs without error ----------------------
expect_warning(
    groupComparisonPlots(data = comparison_result, type = "Heatmap",
                          FCcutoff = 70, logBase.pvalue = 2, address = address_prefix)
)

# Test 5: VolcanoPlot (ggplot2) creates a pdf file, no warning required ------
invisible(capture.output(
    groupComparisonPlots(data = comparison_result, type = "VolcanoPlot",
                          logBase.pvalue = 2, address = address_prefix)
))
expect_true(file.exists(paste0(address_prefix, "VolcanoPlot.pdf")))

# Test 6: VolcanoPlot with FCcutoff and no protein names ---------------------
invisible(capture.output(
    groupComparisonPlots(data = comparison_result, type = "VolcanoPlot",
                          FCcutoff = 70, logBase.pvalue = 2, ylimUp = 100,
                          ProteinName = FALSE, address = address_prefix)
))

# Test 7: ComparisonPlot (ggplot2) creates a pdf file and always warns -------
expect_warning(
    groupComparisonPlots(data = comparison_result, type = "ComparisonPlot",
                          address = address_prefix)
)
expect_true(file.exists(paste0(address_prefix, "ComparisonPlot.pdf")))

unlink(tmp_dir, recursive = TRUE)

# plotly (isPlotly = TRUE) branches --------------------------------------------

# Test 8: VolcanoPlot (plotly), address = FALSE, returns a list of plotly objects
# (address = FALSE requires a single comparison to be selected)
plotly_volcano = invisible(capture.output(
    result_volcano <- groupComparisonPlots(data = comparison_result, type = "VolcanoPlot",
                                            which.Comparison = "T3-T1",
                                            logBase.pvalue = 2, address = FALSE,
                                            isPlotly = TRUE)
))
expect_true(is.list(result_volcano))
expect_true(length(result_volcano) > 0)
expect_true(inherits(result_volcano[[1]], "plotly"))

# Test 9: Heatmap (plotly), saved as a zipped HTML file ----------------------
tmp_dir2 = tempfile("msstats_groupcomparisonplots_heatmap_")
dir.create(tmp_dir2)
address_prefix2 = paste0(tmp_dir2, "/")
invisible(capture.output(suppressWarnings(
    groupComparisonPlots(data = comparison_result, type = "Heatmap",
                          logBase.pvalue = 2, address = address_prefix2,
                          isPlotly = TRUE)
)))
expect_true(any(grepl("Heatmap.*\\.zip$", list.files(tmp_dir2))))
unlink(tmp_dir2, recursive = TRUE)

# Test 10: ComparisonPlot (plotly, all proteins), saved as a zipped HTML file
tmp_dir3 = tempfile("msstats_groupcomparisonplots_comparison_")
dir.create(tmp_dir3)
address_prefix3 = paste0(tmp_dir3, "/")
invisible(capture.output(suppressWarnings(
    groupComparisonPlots(data = comparison_result, type = "ComparisonPlot",
                          which.Protein = "all", address = address_prefix3,
                          isPlotly = TRUE)
)))
expect_true(any(grepl("ComparisonPlot.*\\.zip$", list.files(tmp_dir3))))
unlink(tmp_dir3, recursive = TRUE)
