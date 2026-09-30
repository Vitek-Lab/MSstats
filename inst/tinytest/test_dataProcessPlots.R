# Setup ------------------------------------------------------------------
QuantData = dataProcess(SRMRawData, use_log_file = FALSE)
protein_name = as.character(unique(QuantData$ProteinLevelData$Protein))[1]

# Test 1: invalid type errors -------------------------------------------------
expect_error(dataProcessPlots(QuantData, type = "INVALID"))

# Test 2: address = FALSE with which.Protein = "all" errors -------------------
expect_error(
    dataProcessPlots(QuantData, type = "ProfilePlot", which.Protein = "all",
                      address = FALSE)
)

# Test 3: address = FALSE with multiple proteins errors -----------------------
expect_error(
    dataProcessPlots(QuantData, type = "ProfilePlot", which.Protein = c(1, 2),
                      address = FALSE)
)

# ggplot2 (isPlotly = FALSE) branches, saved to a tempdir ---------------------

tmp_dir = tempfile("msstats_dataprocessplots_")
dir.create(tmp_dir)
address_prefix = paste0(tmp_dir, "/")

# Test 4: ProfilePlot (ggplot2) creates a pdf file, and always warns ----------
expect_warning(
    dataProcessPlots(QuantData, type = "ProfilePlot", which.Protein = protein_name,
                      address = address_prefix,
                      remove_uninformative_feature_outlier = TRUE)
)
expect_true(file.exists(paste0(address_prefix, "ProfilePlot.pdf")))

# Test 5: QCPlot (ggplot2) creates a pdf file ---------------------------------
expect_warning(
    dataProcessPlots(QuantData, type = "QCPlot", which.Protein = protein_name,
                      address = address_prefix)
)
expect_true(file.exists(paste0(address_prefix, "QCPlot.pdf")))

# Test 6: ConditionPlot (ggplot2) creates a pdf file --------------------------
expect_warning(
    dataProcessPlots(QuantData, type = "ConditionPlot", which.Protein = protein_name,
                      address = address_prefix)
)
expect_true(file.exists(paste0(address_prefix, "ConditionPlot.pdf")))

unlink(tmp_dir, recursive = TRUE)

# plotly (isPlotly = TRUE) branches, address = FALSE (no file saving) ---------

# Test 7: ProfilePlot (plotly) returns a list of plotly objects --------------
plotly_profile = suppressWarnings(
    dataProcessPlots(QuantData, type = "ProfilePlot", which.Protein = protein_name,
                      address = FALSE, isPlotly = TRUE)
)
expect_true(is.list(plotly_profile))
expect_true(length(plotly_profile) > 0)
expect_true(inherits(plotly_profile[[1]], "plotly"))

# Test 8: QCPlot (plotly) returns a list of plotly objects -------------------
plotly_qc = suppressWarnings(
    dataProcessPlots(QuantData, type = "QCPlot", which.Protein = protein_name,
                      address = FALSE, isPlotly = TRUE)
)
expect_true(is.list(plotly_qc))
expect_true(length(plotly_qc) > 0)
expect_true(inherits(plotly_qc[[1]], "plotly"))

# Test 9: ConditionPlot (plotly), saved as a zipped HTML file ----------------
# Note: which.Protein must be "all" here (rather than a single protein as
# above) - selecting a subset together with isPlotly = TRUE hits a pre-existing
# bug in dataProcessPlots() where NULL placeholders for unselected proteins
# are still passed to the ggplot-to-plotly conversion step.
tmp_dir2 = tempfile("msstats_dataprocessplots_condition_")
dir.create(tmp_dir2)
address_prefix2 = paste0(tmp_dir2, "/")
invisible(capture.output(suppressWarnings(
    dataProcessPlots(QuantData, type = "ConditionPlot", which.Protein = "all",
                      address = address_prefix2, isPlotly = TRUE)
)))
expect_true(any(grepl("ConditionPlot.*\\.zip$", list.files(tmp_dir2))))
unlink(tmp_dir2, recursive = TRUE)

# Test 10: the Plotly legend is mounted on the right ------------------------

legend_spec = function() {
    plot = suppressWarnings(
        dataProcessPlots(QuantData, type = "ProfilePlot",
                          which.Protein = protein_name, summaryPlot = FALSE,
                          address = FALSE, isPlotly = TRUE)
    )[[1]]
    plotly::plotly_build(plot)$x$layout
}

spec_right = legend_spec()
expect_true(spec_right$showlegend)
expect_equal(spec_right$legend$orientation, "v")
expect_true(spec_right$legend$x > 1)

# Test 11: width sizes the Plotly canvas ------------------------------------

plotly_width = function(...) {
    out = NULL
    invisible(capture.output(suppressWarnings(
        out <- dataProcessPlots(QuantData, type = "QCPlot",
                                which.Protein = protein_name, address = FALSE,
                                isPlotly = TRUE, ...)
    )))
    out[[1]]$x$layout$width
}
expect_equal(plotly_width(), 800)
expect_equal(plotly_width(width = 1400), 1400)

# Test 12: the saved HTML container is sized to the plot it holds ------------

html = as.character(MSstats:::.getPlotlyPlotHTML(list("plot"), 1400, 600))
expect_true(grepl("width:1400px; height:600px", html))
