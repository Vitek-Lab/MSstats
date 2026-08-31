# Condition label layout ------------------------------------------------------
#
# .layoutConditionLabels() decides whether condition names still fit inside the
# panel. Crowding is name width relative to the slot a condition gets, which is
# panel_width / n_conditions -- runs per condition do not affect it.

mk_group_name = function(nms) {
    data.frame(RUN = seq_along(nms), ABUNDANCE = rep(20, length(nms)),
               Name = nms, stringsAsFactors = FALSE)
}
layout_of = function(nms, n_facets = 2, width = 800, text.size = 4,
                     text.angle = 0, adjust = TRUE) {
    MSstats:::.layoutConditionLabels(mk_group_name(nms), n_facets, width,
                                     text.size, text.angle, adjust)
}
long_names = paste0("Timepoint_", sprintf("%02d", 1:10), "_Treated_High_Dose")

# Test 1: short names keep the in-panel layout (the no-regression path)
expect_false(layout_of(as.character(1:10))$use_facets)

# Test 2: eight short names over two facets stay in-panel. This is the shape of
# the dataset the original report pointed at, which renders correctly today.
expect_false(layout_of(c("C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"))$use_facets)

# Test 3: long names switch to facet strips
expect_true(layout_of(long_names)$use_facets)

# Test 4: an explicit text.angle is the caller's choice, so leave the layout be
expect_false(layout_of(long_names, text.angle = 45)$use_facets)

# Test 5: the opt-out is honoured even when crowded
expect_false(layout_of(long_names, adjust = FALSE)$use_facets)

# Test 6: a single condition cannot collide with anything
expect_false(layout_of("A_Single_Very_Long_Condition_Name")$use_facets)

# Test 7: a degenerate canvas width must not error
expect_false(layout_of(long_names, width = 0)$use_facets)

# Test 8: two facets halve the panel, so wrapping is at least as tight as one
expect_true(layout_of(long_names, n_facets = 2)$wrap_chars <=
                layout_of(long_names, n_facets = 1)$wrap_chars)

# Test 9: more conditions means less room, so a smaller strip font
many = paste0("Timepoint_", sprintf("%02d", 1:20), "_Treated_High_Dose")
expect_true(layout_of(many)$strip_size < layout_of(long_names)$strip_size)

# Test 10: a strip never shouts louder than the in-panel label it replaces
expect_true(layout_of(long_names, text.size = 2)$strip_size <= 2 * ggplot2::.pt)

# Test 11: the strip font is floored rather than shrunk into illegibility
crowded = paste0("Condition_", sprintf("%02d", 1:60), "_Treated_High_Dose")
expect_true(layout_of(crowded)$strip_size >= 4)

# Condition label wrapping ----------------------------------------------------

# Test 12: wrapping is lossless at every width, not just the one in use.
# strwrap() breaks on whitespace, so separators are turned into break points and
# the injected spaces removed again; the lines must rejoin to the original name.
lossy_widths = Filter(function(w) {
    any(gsub("\n", "", MSstats:::.wrapConditionLabels(long_names, w),
             fixed = TRUE) != long_names)
}, 4:40)
expect_equal(length(lossy_widths), 0L)

# Test 13: wrapping does insert breaks when the name cannot fit on one line
expect_true(any(grepl("\n", MSstats:::.wrapConditionLabels(long_names, 11),
                      fixed = TRUE)))

# Test 14: a name that already fits is returned untouched
expect_equal(MSstats:::.wrapConditionLabels("Ctrl", 20), "Ctrl")

# Profile plot construction ---------------------------------------------------

prep_profile = function(raw) {
    quant = dataProcess(raw, use_log_file = FALSE)
    processed = data.table::as.data.table(quant$FeatureLevelData)
    summarized = data.table::as.data.table(quant$ProteinLevelData)
    processed$PROTEIN = factor(processed$PROTEIN)
    summarized$Protein = factor(summarized$Protein)
    list(processed = processed, summarized = summarized)
}
# .plotProfile() returns its ggplot objects only when isPlotly = TRUE; the
# conversion to plotly happens later, in dataProcessPlots(), so these are ggplots.
build_profile = function(d, protein, featureName = "Transition",
                         legend.position = "top", legend.ncol = NULL,
                         max.legend.entries = 30, condition.label.adjust = TRUE) {
    plots = MSstats:::.plotProfile(
        d$processed, d$summarized, featureName, FALSE, FALSE, 10, 10, 4, 0, 7, 2,
        800, 600, protein, TRUE, FALSE, FALSE, FALSE, TRUE,
        legend.position, legend.ncol, max.legend.entries, condition.label.adjust)
    plots[["original_plot"]][["plot 1"]]
}
colour_guide_of = function(p) {
    guides = p$guides$guides
    if (!is.null(guides$colour)) guides$colour else guides$color
}
has_geom = function(p, cls) {
    any(vapply(p$layers, function(l) inherits(l$geom, cls), logical(1)))
}

grDevices::pdf(NULL)  # .plotProfile() print()s each plot; swallow the output
dia = prep_profile(DIARawData)
big_protein = "RNA  helicase exp9"  # 149 transitions, the worst bundled case

# Test 15: past max.legend.entries the feature legend is dropped
expect_true(is.null(colour_guide_of(build_profile(dia, big_protein))) ||
                identical(colour_guide_of(build_profile(dia, big_protein)), "none"))

# Test 16: and the count is moved into the title, so the loss is not silent
expect_true(grepl("149 features; legend omitted",
                  build_profile(dia, big_protein)$labels$title, fixed = TRUE))

# Test 17: the comparison is strict -- a protein exactly at the cap keeps its legend
expect_false(grepl("legend omitted",
                   build_profile(dia, big_protein, max.legend.entries = 149)$labels$title))

# Test 18: one below the count, and it is dropped
expect_true(grepl("legend omitted",
                  build_profile(dia, big_protein, max.legend.entries = 148)$labels$title))

# Test 19: an infinite cap never omits
expect_false(grepl("legend omitted",
                   build_profile(dia, big_protein, max.legend.entries = Inf)$labels$title))

# Test 20: featureName = "NA" is the pre-existing full suppression, and must not
# gain an omission note for a legend the user asked not to have
na_plot = build_profile(dia, big_protein, featureName = "NA")
expect_false(grepl("legend omitted", na_plot$labels$title))

# Test 21: the title is the protein name, not a factor level index. PROTEIN is a
# factor and ggplotly() renders a factor title as "1" unless it is coerced.
expect_true(grepl("RNA", build_profile(dia, big_protein)$labels$title))

# Test 22: legend.position reaches the theme
expect_equal(build_profile(dia, big_protein, legend.position = "right")$theme$legend.position,
             "right")

# Test 23: legend.ncol overrides the default of three columns
expect_equal(colour_guide_of(build_profile(dia, "FabG", legend.ncol = 5,
                                           max.legend.entries = Inf))$params$ncol, 5)

# Test 24: and the default is still three
expect_equal(colour_guide_of(build_profile(dia, "FabG",
                                           max.legend.entries = Inf))$params$ncol, 3)

# Facet layout for crowded condition names ------------------------------------

srm_short = prep_profile(SRMRawData)
srm_raw_long = SRMRawData
srm_raw_long$Condition = long_names[as.integer(as.character(srm_raw_long$Condition))]
srm_long = prep_profile(srm_raw_long)

short_plot = build_profile(srm_short, "IDHC")
long_plot = build_profile(srm_long, "IDHC")

# Test 25: uncrowded labels keep the in-panel geom_text layer
expect_true(has_geom(short_plot, "GeomText"))

# Test 26: and still render both isotope label panels
expect_equal(length(unique(ggplot2::ggplot_build(short_plot)$data[[1]]$PANEL)), 2L)

# Test 27: crowded labels become a facet on condition instead
expect_true("GROUP" %in% names(ggplot2::ggplot_build(long_plot)$layout$layout))

# Test 28: the facet replaces the in-panel text rather than doubling up with it
expect_false(has_geom(long_plot, "GeomText"))

# Test 29: the dashed condition separators go with it
expect_false(has_geom(long_plot, "GeomVline"))

# Test 30: opting out keeps the in-panel layout even when crowded
expect_true(has_geom(build_profile(srm_long, "IDHC",
                                   condition.label.adjust = FALSE), "GeomText"))

grDevices::dev.off()
