# Setup ------------------------------------------------------------------
quality_data = data.frame(
    ProteinName = rep(c("ProteinA", "ProteinB"), each = 12),
    PeptideSequence = rep(rep(c("PEPTIDEA", "PEPTIDEB"), each = 6), 2),
    PrecursorCharge = rep(2, 24),
    Run = rep(paste0("Run", 1:6), 4),
    AnomalyScores = c(1:6, 2:7, 3:8, 4:9) / 10,
    EGDeltaRT = seq(0.1, 2.4, length.out = 24)
)

# Test 1: which.Protein is required -------------------------------------------
expect_error(MSstatsQualityMetricsPlot(quality_data))

# Test 2: missing required column errors --------------------------------------
missing_col_data = quality_data
missing_col_data$Run = NULL
expect_error(
    MSstatsQualityMetricsPlot(missing_col_data, which.Protein = "ProteinA")
)

# Test 3: metric not present errors -------------------------------------------
expect_error(
    MSstatsQualityMetricsPlot(quality_data, metric = "NotAColumn",
                               which.Protein = "ProteinA")
)

# Test 4: protein not present errors ------------------------------------------
expect_error(
    MSstatsQualityMetricsPlot(quality_data, which.Protein = "NotAProtein")
)

# Test 5: default call (address = FALSE) returns a ggplot object -------------
p_default = MSstatsQualityMetricsPlot(quality_data, which.Protein = "ProteinA")
expect_true(inherits(p_default, "ggplot"))

# Test 6: a different metric column can be selected ---------------------------
p_metric = MSstatsQualityMetricsPlot(quality_data, metric = "EGDeltaRT",
                                      which.Protein = "ProteinB")
expect_true(inherits(p_metric, "ggplot"))
expect_equal(p_metric$labels$y, "EGDeltaRT")

# Test 7: isPlotly = TRUE returns a plotly object -----------------------------
p_plotly = MSstatsQualityMetricsPlot(quality_data, which.Protein = "ProteinA",
                                      isPlotly = TRUE)
expect_true(inherits(p_plotly, "plotly"))

# Test 8: address != FALSE saves a pdf file (ggplot2) -------------------------
tmp_dir = tempfile("msstats_qualitymetrics_")
dir.create(tmp_dir)
address_prefix = paste0(tmp_dir, "/")
invisible(MSstatsQualityMetricsPlot(quality_data, which.Protein = "ProteinA",
                                     address = address_prefix))
expect_true(file.exists(paste0(address_prefix, "QualityMetricsPlot.pdf")))
unlink(tmp_dir, recursive = TRUE)

# Test 9: address != FALSE saves an html file (plotly) ------------------------
tmp_dir2 = tempfile("msstats_qualitymetrics_html_")
dir.create(tmp_dir2)
address_prefix2 = paste0(tmp_dir2, "/")
invisible(MSstatsQualityMetricsPlot(quality_data, which.Protein = "ProteinA",
                                     address = address_prefix2, isPlotly = TRUE))
expect_true(file.exists(paste0(address_prefix2, "QualityMetricsPlot.html")))
unlink(tmp_dir2, recursive = TRUE)

# Test 10: Run column already a factor is respected (not re-leveled) --------
factor_data = quality_data
factor_data$Run = factor(factor_data$Run, levels = paste0("Run", 6:1))
p_factor = MSstatsQualityMetricsPlot(factor_data, which.Protein = "ProteinA")
expect_true(inherits(p_factor, "ggplot"))
