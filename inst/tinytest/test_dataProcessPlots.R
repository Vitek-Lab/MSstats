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

# plotly (isPlotly = TRUE) does not render plots onto a graphics device --------
old_wd = getwd()
tmp_dir3 = tempfile("msstats_dataprocessplots_nodevice_")
dir.create(tmp_dir3)
setwd(tmp_dir3)

# Test 10: no device is opened and no Rplots.pdf is written -------------------
graphics.off() # start with no open device, so any stray print() would open one
devices_before = dev.list()
invisible(capture.output(suppressWarnings({
    dataProcessPlots(QuantData, type = "ProfilePlot", which.Protein = protein_name,
                      address = FALSE, isPlotly = TRUE)
    dataProcessPlots(QuantData, type = "QCPlot", which.Protein = protein_name,
                      address = FALSE, isPlotly = TRUE)
    dataProcessPlots(QuantData, type = "ConditionPlot", which.Protein = "all",
                      address = paste0(tmp_dir3, "/"), isPlotly = TRUE)
})))
expect_identical(dev.list(), devices_before)
expect_false(file.exists(file.path(tmp_dir3, "Rplots.pdf")))

# Test 11: QCPlot (plotly, saved) leaves the user's open device alone ---------
png(file.path(tmp_dir3, "user_device.png"))
user_device = dev.cur()
invisible(capture.output(suppressWarnings(
    dataProcessPlots(QuantData, type = "QCPlot", which.Protein = protein_name,
                      address = paste0(tmp_dir3, "/"), isPlotly = TRUE)
)))
# Checks the device is still open rather than current: plotly::ggplotly() itself
# can switch the current device when several are open.
expect_true(user_device %in% dev.list())
if (user_device %in% dev.list()) dev.off(user_device)

# tinytest opens a null pdf device per test file and closes it when the file
# ends; Test 10 closed it, so open one again for that cleanup to succeed.
graphics.off()
grDevices::pdf(file = nullfile())
setwd(old_wd)
unlink(tmp_dir3, recursive = TRUE)

# QC plot for PDF output: flat layers built from precomputed boxplot statistics
# Reference: ggplot2::stat_boxplot via layer_data() of a geom_boxplot plot.
# Built like .makeQCPlot's plotly version: one box per RUN within each LABEL.
boxplot_reference = function(input, y_limits) {
    p = ggplot2::ggplot(input, ggplot2::aes(RUN, ABUNDANCE)) +
        ggplot2::facet_grid(~LABEL) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = LABEL)) +
        ggplot2::scale_y_continuous(limits = y_limits)
    ref = suppressWarnings(ggplot2::layer_data(p))
    ref = ref[order(ref$PANEL, ref$x), ]
    list(stats = as.matrix(ref[, c("ymin", "lower", "middle", "upper", "ymax")]),
         x = as.integer(ref$x),
         outliers = sort(unlist(ref$outliers)))
}
expect_stats_match = function(input, y_limits) {
    ref = boxplot_reference(input, y_limits)
    new = MSstats:::.qcBoxStats(data.table::as.data.table(input),
                                y_limits[1], y_limits[2])
    boxes = new$boxes[order(new$boxes$LABEL, new$boxes$x)]
    expect_equal(boxes$x, ref$x)
    expect_equal(unname(as.matrix(boxes[, c("ymin", "lower", "middle", "upper", "ymax")])),
                 unname(ref$stats))
    expect_equal(sort(new$outliers$ABUNDANCE), ref$outliers)
}
qc_input = data.table::as.data.table(QuantData$FeatureLevelData)
qc_input[, RUN := factor(RUN)]
qc_input[, LABEL := factor(LABEL)]
default_limits = c(-1, ceiling(max(qc_input$ABUNDANCE, na.rm = TRUE) + 3))

# Test 12: statistics match ggplot2's boxplot on real data --------------------
expect_stats_match(qc_input, default_limits)

# Test 13: user y-limits drop out-of-range values before the statistics -------
expect_stats_match(qc_input, c(10, 20))

# Test 14: a run with no values and a run with a single value ------------------
edge_input = data.table::copy(qc_input)
runs = levels(edge_input$RUN)
edge_input[RUN == runs[1], ABUNDANCE := NA]
edge_input[RUN == runs[2], ABUNDANCE := c(15, rep(NA, .N - 1))]
expect_stats_match(edge_input, default_limits)

# Test 15: the plotly path keeps geom_boxplot; the PDF path does not ----------
qc_args = list(qc_input, TRUE, default_limits[1], default_limits[2], 10, 10, 4,
               0, 7, c("darkseagreen1", "lightblue"), c(3, 6),
               data.frame(RUN = 2, ABUNDANCE = 1,
               Name = "g"), 3, "Log2-intensities")
layer_geoms = function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
plotly_plot = do.call(MSstats:::.makeQCPlot, c(qc_args, isPlotly = TRUE))
pdf_plot = do.call(MSstats:::.makeQCPlot, c(qc_args, isPlotly = FALSE))
expect_true("GeomBoxplot" %in% layer_geoms(plotly_plot))
expect_false("GeomBoxplot" %in% layer_geoms(pdf_plot))

# Test 16: both versions use the same x-axis range ----------------------------
x_range = function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range
expect_equal(suppressWarnings(x_range(pdf_plot)),
             suppressWarnings(x_range(plotly_plot)))
