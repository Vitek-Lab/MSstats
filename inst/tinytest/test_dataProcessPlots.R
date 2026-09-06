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

# Test 10: every documented legend.position is honoured ----------------------
# Regression test. The placement used to be hard-coded to a right-side vertical
# legend, so "left", "top" and "bottom" were accepted and then ignored -- and
# because the ggplot theme *was* set to the requested side, ggplotly reserved a
# band there that nothing occupied, squeezing the panel. plotly::layout() defers
# into layoutAttrs, so these have to be checked after plotly_build().

legend_spec = function(position) {
    plot = suppressWarnings(
        dataProcessPlots(QuantData, type = "ProfilePlot",
                          which.Protein = protein_name, summaryPlot = FALSE,
                          address = FALSE, isPlotly = TRUE,
                          legend.position = position)
    )[[1]]
    plotly::plotly_build(plot)$x$layout
}

spec_right = legend_spec("right")
expect_true(spec_right$showlegend)
expect_equal(spec_right$legend$orientation, "v")
expect_true(spec_right$legend$x > 1)

spec_left = legend_spec("left")
expect_true(spec_left$showlegend)
expect_equal(spec_left$legend$orientation, "v")
expect_true(spec_left$legend$x < 0)

spec_top = legend_spec("top")
expect_true(spec_top$showlegend)
expect_equal(spec_top$legend$orientation, "h")
expect_true(spec_top$legend$y > 1)

spec_bottom = legend_spec("bottom")
expect_true(spec_bottom$showlegend)
expect_equal(spec_bottom$legend$orientation, "h")
expect_true(spec_bottom$legend$y < 0)

# The four placements have to differ from one another, which is the thing the
# original bug got wrong while still looking correct for the default.
expect_false(isTRUE(all.equal(spec_right$legend, spec_left$legend)))
expect_false(isTRUE(all.equal(spec_top$legend, spec_bottom$legend)))
expect_false(isTRUE(all.equal(spec_right$legend, spec_top$legend)))

# Test 11: legend.position = "none" hides the legend outright ----------------
# Both the layout flag and every trace, because the post-processing helpers turn
# individual traces back on after conversion. .fixCensoredPointsLegendProfile-
# PlotsPlotly() in particular re-enables the detected and censored entries, so a
# layout flag on its own leaves traces marked visible underneath it.
plot_none = suppressWarnings(
    dataProcessPlots(QuantData, type = "ProfilePlot",
                      which.Protein = protein_name, summaryPlot = FALSE,
                      address = FALSE, isPlotly = TRUE,
                      legend.position = "none")
)[[1]]
built_none = plotly::plotly_build(plot_none)
expect_false(built_none$x$layout$showlegend)
expect_true(all(!vapply(built_none$x$data,
                        function(trace) isTRUE(trace$showlegend), logical(1))))

# Test 12: an undocumented legend.position is rejected rather than ignored ----
expect_error(
    dataProcessPlots(QuantData, type = "ProfilePlot",
                      which.Protein = protein_name, address = FALSE,
                      isPlotly = TRUE, legend.position = "middle")
)
