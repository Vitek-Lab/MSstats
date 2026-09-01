#' Get name for y-axis
#' @param temp data.table
#' @keywords internal
.getYaxis = function(temp) {
    INTENSITY = ABUNDANCE = NULL
    
    temp = temp[!is.na(INTENSITY) & !is.na(ABUNDANCE),]
    temp_abund = temp[1, "ABUNDANCE"]
    temp_inten = temp[1, "INTENSITY"]
    log2_diff = abs(log(temp_inten, 2) - temp_abund)
    log10_diff = abs(log(temp_inten, 10) - temp_abund)
    if (log2_diff < log10_diff) {
        "Log2-intensities"
    } else {
        "Log10-intensities"
    }
}

#' Get data for a single protein to plot
#' @param dataProcess output -> FeatureLevelData
#' @param all_proteins character, set of protein names
#' @param i integer, index of protein to use
#' @keywords internal
.getSingleProteinForProfile = function(processed, all_proteins, i) {
    FEATURE = SUBJECT = GROUP = PEPTIDE = NULL
    
    single_protein = processed[processed$PROTEIN == all_proteins[i], ]
    single_protein[, FEATURE := factor(FEATURE)]
    single_protein[, SUBJECT := factor(SUBJECT)]
    single_protein[, GROUP := factor(GROUP)]
    single_protein[, PEPTIDE := factor(PEPTIDE)]
    single_protein
}


#' Decide how condition labels should be laid out
#'
#' Condition names are drawn inside the panel by `geom_text()`, one per
#' condition, and conditions tile the panel evenly. Each name gets
#' `panel_width / n_conditions` of horizontal room no matter how many runs it
#' contains -- runs per condition do not affect label spacing at all, which is
#' why crowding is measured as name width relative to that slot and never as a
#' sample count.
#'
#' Once a name is wider than its slot, the in-panel layout cannot be rescued by
#' moving text about. Rotation does not survive `ggplotly()` and MSstatsShiny
#' renders through it; `ggrepel` drops labels outright for the same reason;
#' shrinking the font enough to fit 30-character names reaches ~2pt, which is
#' unreadable. Wrapping and staggering inside the panel do make the names
#' legible, but only by printing them on top of the data.
#'
#' So past the threshold the condition becomes a facet instead, and ggplot2
#' reserves a strip band for it -- collision with the data and clipping at the
#' panel edge stop being possible rather than being mitigated. Verified to
#' survive `ggplotly()`: the strips, `space = "free_x"` sizing, and newline
#' wrapping (which is translated to `<br />`).
#'
#' @param groupName data.frame of RUN, ABUNDANCE and Name, as built by
#'   `.plotProfile()` and `.plotQC()`
#' @param n_facets number of facet panels actually drawn. Pass
#'   `length(unique(input$LABEL))`, not `nlevels()`: LABEL is a factor over the
#'   whole table, so `nlevels()` reports 2 for a protein carrying only one label
#'   while `facet_grid()` draws a single panel, and the labels would then be
#'   given half the room they really have.
#' @param width width of the canvas in pixels
#' @param text.size size of the condition labels
#' @param text.angle angle of the condition labels
#' @param condition.label.adjust if FALSE, keep the in-panel layout unchanged
#' @param isPlotly TRUE when the plot is bound for `ggplotly()`. The canvas is
#'   the same number, but not the same size: on the pdf device `width` is points
#'   at 72dpi, while plotly treats it as CSS pixels at 96dpi, which is a quarter
#'   less room. Measuring both at 72dpi under-triggers the facet layout in the
#'   browser -- an 8-condition SILAC design scores 0.80 on the pdf geometry and
#'   1.11 on the real one, so it stayed in-panel and the labels collided.
#'
#' @return list with `use_facets`, and when that is TRUE the `wrap_chars`
#'   width to wrap strip labels at and the `strip_size` to draw them at
#' @keywords internal
.layoutConditionLabels = function(
    groupName, n_facets, width, text.size, text.angle,
    condition.label.adjust = TRUE, isPlotly = FALSE
) {
    inline = list(use_facets = FALSE, wrap_chars = NA_integer_,
                  strip_size = NA_real_)
    if (!isTRUE(condition.label.adjust)) {
        return(inline)
    }
    # A non-zero text.angle is a deliberate choice by the caller. Honour it and
    # leave the layout alone rather than overriding it.
    if (!isTRUE(all.equal(as.numeric(text.angle), 0))) {
        return(inline)
    }
    n_conditions = nrow(groupName)
    if (n_conditions < 2L || !is.numeric(width) || width <= 0) {
        return(inline)
    }
    # Horizontal room per condition. On the pdf device the canvas is width/72
    # inches (savePlot() converts the same way); plotly reads the same number as
    # CSS pixels at 96dpi. ~1.1in of it goes to the y-axis title, tick labels and
    # margins, and what is left is split across the facets.
    dpi = if (isTRUE(isPlotly)) 96 else 72
    panel_in = (width / dpi - 1.1) / max(n_facets, 1L)
    slot_in = panel_in / n_conditions
    if (slot_in <= 0) {
        return(inline)
    }
    # Label width is estimated from nchar, not measured. grid::stringWidth() is
    # exact but needs an open graphics device, which is not available while the
    # plot is being built -- measuring would make the result device-dependent
    # and this function untestable. 0.53 em per character is calibrated against
    # graphics::strwidth() on the reproduction data and lands within ~7%.
    char_in = text.size * ggplot2::.pt * 0.53 / 72
    if (max(nchar(as.character(groupName$Name))) * char_in < slot_in) {
        # The name is narrower than its slot. Every dataset that renders
        # correctly today takes this path and is left exactly as it was.
        return(inline)
    }
    # Size the strip text to the room a strip actually has. strwrap() cannot
    # break inside a word, so the binding constraint is the longest unbreakable
    # token, not the longest name: pick the largest font at which that token
    # still fits its slot. Bounded above by the caller's text.size (a strip
    # should never shout louder than the in-panel labels it replaces) and below
    # by 4pt, past which shrinking buys illegibility rather than fit.
    tokens = unlist(strsplit(gsub("([_.])", "\\1 ", as.character(groupName$Name)),
                             "[[:space:]]+"), FALSE, FALSE)
    longest_token = max(nchar(tokens), 1L)
    fitted_pt = slot_in * 72 / (longest_token * 0.53)
    strip_size = max(4, min(text.size * ggplot2::.pt, fitted_pt))
    # Wrap at the size the strip will actually be drawn at, not the in-panel one.
    strip_char_in = strip_size * 0.53 / 72
    list(use_facets = TRUE,
         wrap_chars = max(1L, floor(slot_in / strip_char_in)),
         strip_size = strip_size)
}


#' Wrap condition names for use as facet strip labels
#'
#' `strwrap()` breaks only at whitespace and condition names are usually
#' underscore-delimited, so separators are turned into break opportunities and
#' the injected spaces are removed again afterwards. The lines therefore rejoin
#' to the original name at any wrap width.
#'
#' @param x character vector of condition names
#' @param wrap_chars width to wrap at, in characters
#'
#' @return character vector with newlines inserted
#' @keywords internal
.wrapConditionLabels = function(x, wrap_chars) {
    vapply(as.character(x), function(nm) {
        parts = strwrap(gsub("([_.])", "\\1 ", nm), width = wrap_chars + 1L)
        paste(gsub("([_.]) ", "\\1", parts), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}


#' Facet, separator and label pieces for one condition layout
#'
#' Returns the parts that differ between the in-panel layout and the facet
#' layout so the plot builders can add them unconditionally: ggplot2 treats the
#' addition of NULL as a no-op, which keeps one pipeline instead of two.
#'
#' @param layout result of `.layoutConditionLabels()`
#' @param groupName data.frame of label positions, for the in-panel layout
#' @param lineNameAxis positions of the condition separator lines
#' @param text.size size of the in-panel condition labels
#' @param text.angle angle of the in-panel condition labels
#' @param x.axis.size size of the x axis text
#'
#' @return list of ggplot2 objects or NULL, named facet, vline, text and strip
#' @keywords internal
.conditionLayoutSpecs = function(layout, groupName, lineNameAxis, text.size,
                                 text.angle, x.axis.size) {
    GROUP = LABEL = NULL
    
    if (!isTRUE(layout$use_facets)) {
        return(list(
            facet = facet_grid(~LABEL),
            vline = geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey",
                               linetype = "longdash"),
            text = geom_text(data = groupName,
                             aes(x = .data$RUN, y = .data$ABUNDANCE,
                                 label = .data$Name),
                             size = text.size, angle = text.angle,
                             color = "black"),
            strip = NULL))
    }
    # The dashed separators and the in-panel names are exactly what the facet
    # replaces, so both are dropped rather than drawn twice.
    list(
        # scales = "free_x" gives each condition its own run axis. space =
        # "free_x" would also size panels by run count, but it computes a
        # non-finite panel width when a protein is absent from a condition and
        # the x scale is discrete (the QC plot), so grid fails with
        # "non-finite location and/or size for viewport". Equal-width panels
        # cost nothing here and cannot produce that.
        facet = facet_grid(LABEL ~ GROUP, scales = "free_x",
                           labeller = labeller(GROUP = function(x)
                               .wrapConditionLabels(x, layout$wrap_chars))),
        vline = NULL,
        text = NULL,
        strip = theme(strip.text.x = element_text(size = layout$strip_size),
                      panel.spacing.x = unit(1.5, "pt"),
                      axis.text.x = element_text(size = max(x.axis.size - 5, 4))))
}


#' Create profile plot
#' @inheritParams dataProcessPlots
#' @param input data.table
#' @param is_censored TRUE if censored values were imputed
#' @keywords internal
.makeProfilePlot = function(
    input, is_censored, featureName, y.limdown, y.limup, x.axis.size, 
    y.axis.size, text.size, text.angle, legend.size, dot.size.profile, 
    ss, s, cumGroupAxis, yaxis.name, lineNameAxis, groupNametemp, dot_colors,
    legend.position = "top", legend.ncol = NULL, max.legend.entries = 30,
    width = 800, condition.label.adjust = TRUE, isPlotly = FALSE
) {
    RUN = ABUNDANCE = Name = NULL
    
    if (is_censored) {
        input$is_censored = factor(input$is_censored, 
                                   levels = c("FALSE", "TRUE"))
    }
    featureName = toupper(featureName)
    if (featureName == "TRANSITION") {
        type_color = "FEATURE"
    } else {
        type_color = "PEPTIDE"
    }
    
    # ggplot2 allocates legend space before panel space, and the number of
    # features per protein is unbounded, so a large legend can shrink the panel
    # to nothing: DIARawData's 149-feature protein renders as a full page of
    # legend with no plot on it at all. Past `max.legend.entries` the colour key
    # is dropped and the count moves into the title, so the loss is stated
    # rather than silent -- the plotly path already loses entries silently, and
    # reproducing that would be no fix.
    #
    # Why 30: over 105 proteins of a real Spectronaut export (median 11 features
    # per protein, p90 37, max 107) a cap of 30 leaves ~86% of proteins
    # untouched; 40 leaves ~91%, 50 leaves ~95%. The comparison is strict, so 30
    # entries keep the legend and 31 lose it.
    #
    # The counter-argument, reviewed and rejected 2026-08-30: at the boundary the
    # legend is not yet doing harm. DDARawData's "rabbit" (31 features) renders a
    # readable panel *with* its legend; omitting it buys ~44% panel height and
    # costs all 31 feature identities, so 40 or 50 would also have been
    # defensible. 30 was kept deliberately. If you are here because a reviewer
    # asked to raise it, the numbers above are the argument -- it is a one-line
    # change and nothing else depends on the value.
    #
    # Note the cap rarely binds on the plotly path: ggplotly() flattens the
    # legend to one column and truncates it to ~10 entries before the cap is
    # reached. There the title note, not the cap, is what fixes the bug.
    n_legend_entries = data.table::uniqueN(input[[type_color]])
    omit_legend = featureName != "NA" && is.numeric(max.legend.entries) &&
        n_legend_entries > max.legend.entries
    # as.character matters: PROTEIN is a factor, and ggplotly() renders a factor
    # title as its level index ("1") rather than the protein name.
    plot_title = as.character(unique(input$PROTEIN))
    if (omit_legend) {
        plot_title = paste0(plot_title, " (", n_legend_entries,
                            " features; legend omitted)")
    }
    legend_ncol = if (is.null(legend.ncol)) 3 else legend.ncol
    
    cond_specs = .conditionLayoutSpecs(
        .layoutConditionLabels(groupNametemp, length(unique(input$LABEL)), width,
                               text.size, text.angle, condition.label.adjust,
                               isPlotly),
        groupNametemp, lineNameAxis, text.size, text.angle, x.axis.size)
    
    profile_plot = ggplot(data = input, aes(x = .data$RUN, y = .data$newABUNDANCE,
                                            color = .data[[type_color]], linetype = .data$FEATURE)) +
        cond_specs$facet +
        geom_line(linewidth = 0.5)
    
    if (is_censored) {
        profile_plot = profile_plot +
        geom_point(aes(x = .data$RUN, y = .data$newABUNDANCE, color = .data[[type_color]], shape = .data$censored),
                   data = input,
                   size = dot.size.profile) +
        scale_shape_manual(values = c(16, 1),
                           labels = c("Detected data", "Censored missing data"))
    } else {
        profile_plot = profile_plot +
            geom_point(size = dot.size.profile) +
            scale_shape_manual(values = c(16))
    }
    
    
    if (featureName == "TRANSITION") {
        profile_plot = profile_plot +
            scale_colour_manual(values = dot_colors[s])
    } else if (featureName == "PEPTIDE") {
        profile_plot = profile_plot +
            scale_colour_manual(values = dot_colors[seq_along(unique(s))])
    } else if (featureName == "NA") {
        if (is_censored) {
            profile_plot = profile_plot +
                scale_colour_manual(values = dot_colors[seq_along(unique(s))])
        } else {
            profile_plot = profile_plot +
                scale_colour_manual(values = dot_colors[s])
        }
    }
    
    profile_plot = profile_plot + scale_linetype_manual(values = ss, guide = "none") 
    profile_plot = profile_plot +
        scale_x_continuous('MS runs', breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        cond_specs$vline +
        labs(title = plot_title) +
        cond_specs$text +
        theme_msstats("PROFILEPLOT", x.axis.size, y.axis.size, legend.size,
                      legend_position = legend.position) +
        cond_specs$strip
    
    if (featureName == "TRANSITION") {
        color_guide = guide_legend(order=1,
                                   override.aes = list(size=1.2,
                                                       linetype = ss),
                                   title = paste("# peptide:", nlevels(input$PEPTIDE)), 
                                   title.theme = element_text(size = 13, angle = 0),
                                   keywidth = 0.25,
                                   keyheight = 0.1,
                                   default.unit = 'inch',
                                   ncol = legend_ncol)
    } else if (featureName == "PEPTIDE") {
        color_guide = guide_legend(order=1,
                                   title = paste("# peptide:", nlevels(input$PEPTIDE)), 
                                   title.theme = element_text(size = 13, angle = 0),
                                   keywidth = 0.25,
                                   keyheight = 0.1,
                                   default.unit = 'inch',
                                   ncol = legend_ncol)
    }
    shape_guide = guide_legend(order=2,
                               title = NULL,
                               label.theme = element_text(size = 11, angle = 0),
                               keywidth = 0.1,
                               keyheight = 0.1,
                               default.unit = 'inch')
    if (omit_legend) {
        # Keep the (small) censoring legend; only the per-feature colour key goes.
        if (is_censored) {
            profile_plot = profile_plot + guides(color = "none",
                                                 shape = shape_guide)
        } else {
            profile_plot = profile_plot + guides(color = "none")
        }
    } else if (is_censored) {
        if (featureName == "NA") {
            profile_plot = profile_plot + guides(color = "none",
                                                 shape = shape_guide)
        } else {
            profile_plot = profile_plot + guides(color = color_guide,
                                                 shape = shape_guide)
        } 
    } else {
        profile_plot = profile_plot + guides(color = color_guide)
    }
    profile_plot    
}


#' Make summary profile plot
#' @inheritParams dataProcessPlots
#' @inheritParams .makeProfilePlot
#' @keywords internal
.makeSummaryProfilePlot = function(
    input, is_censored, y.limdown, y.limup, x.axis.size, y.axis.size, 
    text.size, text.angle, legend.size, dot.size.profile, cumGroupAxis, 
    yaxis.name, lineNameAxis, groupNametemp, legend.position = "top",
    width = 800, condition.label.adjust = TRUE, isPlotly = FALSE
) {
    RUN = ABUNDANCE = Name = NULL
    
    cond_specs = .conditionLayoutSpecs(
        .layoutConditionLabels(groupNametemp, length(unique(input$LABEL)), width,
                               text.size, text.angle, condition.label.adjust,
                               isPlotly),
        groupNametemp, lineNameAxis, text.size, text.angle, x.axis.size)
    
    num_features = data.table::uniqueN(input$FEATURE)
    profile_plot = ggplot(data = input, 
                          aes(x = .data$RUN, y = .data$newABUNDANCE, 
                                     color = .data$analysis, linetype = .data$FEATURE, 
                                     size = .data$analysis)) +
        cond_specs$facet +
        geom_line(linewidth = 0.5)
    
    if (is_censored) { # splitting into two layers to keep red above grey
        profile_plot = profile_plot +
            geom_point(data = input[input$PEPTIDE != "Run summary"], 
                       aes(x = .data$RUN, y = .data$newABUNDANCE, 
                                  color = .data$analysis, size = .data$analysis, 
                                  shape = .data$censored)) +
            geom_point(data = input[input$PEPTIDE == "Run summary"], 
                       aes(x = .data$RUN, y = .data$newABUNDANCE, 
                                  color = .data$analysis, size = .data$analysis, 
                                  shape = .data$censored)) +
            geom_errorbar(data = input[input$PEPTIDE == "Run summary"],
                          aes(x = .data$RUN, 
                              ymin = .data$LOWERBOUND, 
                              ymax = .data$UPPERBOUND,
                              color = .data$analysis),
                          width = 0.3,
                          linewidth = 0.5,
                          linetype = "solid") + 
            scale_shape_manual(values = c(16, 1), 
                               labels = c("Detected data",
                                          "Censored missing data"))
    } else {
        profile_plot = profile_plot +         
            geom_point(size = dot.size.profile) +
            scale_shape_manual(values = c(16))
    }
    
    profile_plot  =  profile_plot +
        scale_colour_manual(values = c("lightgray", "darkred")) +
        scale_size_manual(values = c(1.7, 2), guide = "none") +
        scale_linetype_manual(values = c(rep(1, times = num_features - 1), 2), 
                              guide = "none") +
        scale_x_continuous("MS runs", breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        cond_specs$vline +
        labs(title = unique(input$PROTEIN)) +
        cond_specs$text +
        theme_msstats("PROFILEPLOT", x.axis.size, y.axis.size, 
                      legend.size, legend_position = legend.position,
                      legend.title = element_blank()) +
        cond_specs$strip
    color_guide  =  guide_legend(order = 1,
                                 title = NULL,
                                 label.theme = element_text(size = 10, angle = 0))
    shape_guide  =  guide_legend(order = 2, 
                                 title = NULL,
                                 label.theme = element_text(size = 10, angle = 0))
    if (is_censored) {
        profile_plot = profile_plot +
            guides(color = color_guide, shape = shape_guide)
    } else {
        profile_plot = profile_plot +
            guides(color = color_guide) +
            geom_point(aes(x = .data$RUN, y = .data$newABUNDANCE, size = .data$analysis,
                                  color = .data$analysis), data = input)
    }
    profile_plot
}


#' Make QC plot
#' @inherit dataProcessPlots
#' @param input data.table
#' @param all_proteins character vector of protein names
#' @keywords internal
.makeQCPlot = function(
    input, all_proteins, y.limdown, y.limup, x.axis.size, y.axis.size, 
    text.size, text.angle, legend.size, label.color, cumGroupAxis, groupName,
    lineNameAxis, yaxis.name, width = 800, condition.label.adjust = TRUE,
    isPlotly = FALSE
) { 
    RUN = ABUNDANCE = Name = NULL
    
    cond_specs = .conditionLayoutSpecs(
        .layoutConditionLabels(groupName, length(unique(input$LABEL)), width,
                               text.size, text.angle, condition.label.adjust,
                               isPlotly),
        groupName, lineNameAxis, text.size, text.angle, x.axis.size)
    
    if (all_proteins) {
        plot_title = "All"
    } else {
        plot_title = unique(input$PROTEIN)
    }
    
    ggplot(input, aes(x = .data$RUN, y = .data$ABUNDANCE)) +
        cond_specs$facet +
        geom_boxplot(aes(fill = .data$LABEL), outlier.shape = 1,
                     outlier.size = 1.5) +
        scale_fill_manual(values = label.color, guide = "none") +
        scale_x_discrete("MS runs", breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        cond_specs$vline +
        labs(title  =  plot_title) +
        cond_specs$text +
        theme_msstats("QCPLOT", x.axis.size, y.axis.size,
                      legend_size = NULL) +
        cond_specs$strip
    
}


#' Make condition plot
#' @inheritParams dataProcessPlots
#' @param input data.table
#' @param single_protein data.table
#' @keywords internal
.makeConditionPlot = function(
    input, scale, single_protein, y.limdown, y.limup, x.axis.size, y.axis.size, 
    text.size, text.angle, legend.size, dot.size.condition, yaxis.name
) {
    Mean = ciw = NULL
    
    colnames(input)[colnames(input) == "GROUP"] = "Label"
    if (scale) {
        input$Label = as.numeric(gsub("\\D", "", unique(input$Label)))
    }
    
    plot = ggplot(aes(x = .data$Label, y = .data$Mean), data = input) +
        geom_errorbar(aes(ymax = .data$Mean + .data$ciw, ymin = .data$Mean - .data$ciw),
                      data = input, width = 0.1, colour = "red") +
        geom_point(size = dot.size.condition, colour = "darkred")
    
    if (!scale) {
        plot = plot + scale_x_discrete("Condition")
    } else {
        plot = plot + scale_x_continuous("Condition", breaks = input$Label, 
                                         labels = input$Label)
    }
    
    plot = plot +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_hline(yintercept = 0, linetype = "twodash", 
                   colour = "darkgrey", linewidth = 0.6) +
        labs(title = unique(single_protein$PROTEIN)) +
        theme_msstats("CONDITIONPLOT", x.axis.size, y.axis.size, 
                      text_angle = text.angle)
    plot
}
