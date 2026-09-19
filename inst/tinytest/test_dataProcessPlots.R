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
# Regression test. The placement used to disagree with the ggplot theme, so
# ggplotly reserved a band that nothing occupied and squeezed the panel into the
# corner. plotly::layout() defers into layoutAttrs, so this has to be checked
# after plotly_build().

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

# Test 11: text.angle no longer suppresses the Plotly label layout -----------
# ggplotly() does not carry geom_text() rotation through, so a rotated Plotly
# plot is drawn horizontally and has exactly the crowding problem the layout
# exists to solve. The layout therefore runs whatever text.angle says, and the
# untruncated name stays on hover.

QuantDataLong = QuantData
long_group = function(x) factor(paste0("Study_Tissue_Timepoint_", x))
QuantDataLong$FeatureLevelData$GROUP = long_group(QuantDataLong$FeatureLevelData$GROUP)
QuantDataLong$ProteinLevelData$GROUP = long_group(QuantDataLong$ProteinLevelData$GROUP)

condition_label_trace = function(text.angle) {
    plot = suppressWarnings(
        dataProcessPlots(QuantDataLong, type = "ProfilePlot",
                          which.Protein = protein_name, summaryPlot = FALSE,
                          address = FALSE, isPlotly = TRUE,
                          text.angle = text.angle)
    )[[1]]
    traces = plotly::plotly_build(plot)$x$data
    Filter(function(trace) identical(trace$mode, "text"), traces)[[1]]
}

for (angle in c(0, 90)) {
    trace = condition_label_trace(angle)
    expect_true(all(nchar(trace$text) < nchar(trace$hovertext)))
    expect_true(all(grepl("^Study_Tissue_Timepoint_", trace$hovertext)))
}

# Test 12: the ggplot2/PDF path keeps the full names and honours text.angle --
# Nothing is laid out there, so condition.layout is NULL. This exercises that
# branch: the labels fall back to the full condition names and no headroom is
# added.

tmp_dir3 = tempfile("msstats_dataprocessplots_pdf_")
dir.create(tmp_dir3)
expect_silent(suppressWarnings(
    dataProcessPlots(QuantDataLong, type = "ProfilePlot",
                      which.Protein = protein_name, summaryPlot = FALSE,
                      address = paste0(tmp_dir3, "/"), text.angle = 90)
))
expect_true(any(grepl("ProfilePlot.*\\.pdf$", list.files(tmp_dir3))))
unlink(tmp_dir3, recursive = TRUE)

# Test 13: there is no separate Plotly width argument ------------------------
# The Plotly canvas width is an internal constant, not something the caller
# sizes; width is the PDF page.

expect_false("width.plotly" %in% names(formals(dataProcessPlots)))

# Test 14: the saved HTML container is sized to the plot it holds ------------
# Regression test. The container was pinned at 800px while the widget inside it
# was 1400px wide, so the right-hand side of every saved plot -- which is where
# the feature legend is mounted -- fell outside the box.

tmp_dir4 = tempfile("msstats_dataprocessplots_html_")
dir.create(tmp_dir4)
invisible(capture.output(suppressWarnings(
    dataProcessPlots(QuantData, type = "ProfilePlot", which.Protein = protein_name,
                      summaryPlot = FALSE, address = paste0(tmp_dir4, "/"),
                      isPlotly = TRUE)
)))
zip_path = list.files(tmp_dir4, pattern = "\\.zip$", full.names = TRUE)[1]
unzip(zip_path, exdir = file.path(tmp_dir4, "unzipped"))
html_path = list.files(file.path(tmp_dir4, "unzipped"), pattern = "\\.html$",
                       full.names = TRUE, recursive = TRUE)[1]
html = paste(readLines(html_path, warn = FALSE), collapse = "\n")

# The container div spaces its declarations, the widget div does not, so the
# two are told apart by the space after the semicolon.
px = function(pattern) {
    as.integer(sub("^width:([0-9]+)px.*", "\\1",
                   regmatches(html, regexpr(pattern, html))))
}
container_width = px("width:[0-9]+px; height:[0-9]+px; margin")
widget_width = px("width:[0-9]+px;height:[0-9]+px")

expect_equal(container_width, widget_width)
expect_true(container_width >= 1400L)
unlink(tmp_dir4, recursive = TRUE)
