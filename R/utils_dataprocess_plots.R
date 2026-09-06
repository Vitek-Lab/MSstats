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


#' Drop the prefix that every condition name shares
#'
#' Condition names in real designs are usually built from a common stem plus a
#' distinguishing tail -- "Cyno_Colon_Timepoint_0hr", "Cyno_Colon_Timepoint_12hrs".
#' Only the tail identifies the block, but the shared stem is what consumes the
#' horizontal room, so it is dropped from the in-panel label and reported once in
#' the x-axis title instead. No information is lost from the static image.
#'
#' @param names character, condition names in plotting order
#' @return list with `labels` (shortened) and `prefix` (what was removed, "" when
#'   nothing is shared)
#' @keywords internal
.stripCommonAffix = function(names) {
    names = as.character(names)
    unchanged = list(labels = names, prefix = "")
    if (length(unique(names)) < 2L) {
        return(unchanged)
    }
    # Split after each separator so the separator stays with the token it follows
    # and the pieces can simply be pasted back together.
    tokens = strsplit(names, "(?<=[_.[:space:]-])", perl = TRUE)
    n_shared = 0L
    repeat {
        # Never consume a name entirely; a condition with no label left would be
        # indistinguishable from its neighbours.
        nth = vapply(tokens, function(x) {
            if (length(x) > n_shared + 1L) x[n_shared + 1L] else NA_character_
        }, character(1))
        if (anyNA(nth) || length(unique(nth)) != 1L) {
            break
        }
        n_shared = n_shared + 1L
    }
    if (n_shared == 0L) {
        return(unchanged)
    }
    list(labels = vapply(tokens, function(x) {
             paste(x[-seq_len(n_shared)], collapse = "")
         }, character(1)),
         prefix = paste(tokens[[1]][seq_len(n_shared)], collapse = ""))
}


#' Number of characters that fit in one condition's slot
#'
#' Conditions tile the panel evenly, so each name gets `panel_width /
#' n_conditions` of room no matter how many runs it covers -- which is why
#' crowding is a function of name length and condition count, and never of the
#' number of samples per condition.
#'
#' Width is estimated from `nchar` rather than measured. `grid::stringWidth()` is
#' exact but needs an open graphics device, which is not available while the plot
#' is being built; measuring would make the layout device-dependent and this
#' function untestable. 0.53 em per character is calibrated against
#' `graphics::strwidth()` and lands within ~7%.
#'
#' @param n_conditions number of conditions
#' @param n_facets number of facet panels actually drawn. Pass
#'   `length(unique(input$LABEL))`, not `nlevels()`: LABEL is a factor over the
#'   whole table, so `nlevels()` reports 2 for a protein carrying only one label
#'   while `facet_grid()` draws a single panel.
#' @param width width of the canvas in pixels, read as CSS pixels at 96dpi
#' @param text.size size of the condition labels
#' @return integer, at least 1
#' @keywords internal
.conditionSlotChars = function(n_conditions, n_facets, width, text.size) {
    if (!is.numeric(width) || length(width) != 1L || is.na(width) ||
        width <= 0 || n_conditions < 1L) {
        return(.Machine$integer.max)
    }
        return(.Machine$integer.max)
    }
    # ~1.1in of the canvas goes to the y-axis title, tick labels and margins;
    # what is left is split across the facets and then across the conditions.
    panel_in = (width / 96 - 1.1) / max(n_facets, 1L)
    # Only fill part of the slot: a label that fills it exactly touches its
    # neighbours, and the first and last labels overhang the panel edge because
    # they are centred on their block.
    slot_in = 0.85 * panel_in / n_conditions
    char_in = text.size * ggplot2::.pt * 0.53 / 72
    if (slot_in <= 0 || char_in <= 0) {
        return(1L)
    }
    max(1L, as.integer(floor(slot_in / char_in)))
}


#' Wrap condition names onto several lines so they fit their slot
#'
#' Used only for what `.stripCommonAffix()` and shrinking the font cannot fix.
#' `strwrap()` breaks only at whitespace and condition names are usually
#' underscore-delimited, so separators are turned into break opportunities here.
#' A single token wider than the slot cannot be broken and is truncated; the
#' untruncated name stays available in the Plotly hover.
#'
#' @param names character, condition names
#' @param chars maximum characters per line
#' @return character, `names` unchanged when they all already fit
#' @keywords internal
.wrapConditionLabels = function(names, chars) {
    names = as.character(names)
    if (all(nchar(names) <= chars)) {
        return(names)
    }
    vapply(names, function(name) {
        tokens = regmatches(name, gregexpr("[^_.[:space:]-]+[_.[:space:]-]*",
                                           name))[[1]]
        if (length(tokens) == 0L) {
            tokens = name
        }
        tokens = vapply(tokens, function(token) {
            if (nchar(token) > chars) {
                paste0(substr(token, 1L, max(1L, chars - 3L)), "...")
            } else {
                token
            }
        }, character(1), USE.NAMES = FALSE)
        lines = character(0)
        current = ""
        for (token in tokens) {
            candidate = paste0(current, token)
            if (nchar(trimws(candidate)) > chars && nzchar(current)) {
                lines = c(lines, current)
                current = token
            } else {
                current = candidate
            }
        }
        paste(c(lines, current), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}


#' Lay out condition labels so they do not overlap
#'
#' Applies the three mitigations in order of how much they cost the reader:
#' drop the shared stem, then shrink the font, then wrap. Each is a no-op when
#' the labels already fit, so a plot that renders correctly today is unchanged.
#'
#' @inheritParams .conditionSlotChars
#' @param names character, condition names in plotting order
#' @param text.angle angle of the labels. A non-zero value is a deliberate choice
#'   by the caller, so the layout is left alone. Note that rotation is not carried
#'   through by `ggplotly()`, so it does not help the MSstatsShiny output.
#' @return list with `labels`, the `size` to draw them at, the `n_lines` they
#'   occupy, and the `xaxis` title to use
#' @keywords internal
.layoutConditionLabels = function(names, n_facets, width, text.size,
                                  text.angle = 0) {
    labels = as.character(names)
    unchanged = list(labels = labels, size = text.size, n_lines = 1L,
                     xaxis = "MS runs")
    if (!isTRUE(all.equal(as.numeric(text.angle), 0))) {
        return(unchanged)
    }
    n_conditions = length(labels)
    if (n_conditions < 2L) {
        return(unchanged)
    }
    if (max(nchar(labels)) <=
        .conditionSlotChars(n_conditions, n_facets, width, text.size)) {
        return(unchanged)
    }
    xaxis = "MS runs"
    stripped = .stripCommonAffix(labels)
    if (nzchar(stripped$prefix)) {
        labels = stripped$labels
        xaxis = paste0("MS runs   (conditions: ", stripped$prefix, "*)")
    }
    # Shrink before wrapping: one legible line beats two cramped ones. The floor
    # is where shrinking stops buying fit and starts buying illegibility.
    size = text.size
    repeat {
        chars = .conditionSlotChars(n_conditions, n_facets, width, size)
        if (max(nchar(labels)) <= chars || size <= 2.5) {
            break
        }
        size = size - 0.25
    }
    labels = .wrapConditionLabels(labels, chars)
    list(labels = labels, size = size,
         n_lines = max(lengths(strsplit(labels, "\n", fixed = TRUE))),
         xaxis = xaxis)
}


#' Accessors for the condition label layout
#'
#' The builders are also called with `condition.layout = NULL` (nothing computed
#' a layout), in which case they fall back to the historical behaviour.
#' @param layout result of `.layoutConditionLabels()`, or NULL
#' @keywords internal
.conditionXlab = function(layout) {
    if (is.null(layout$xaxis)) "MS runs" else layout$xaxis
}

#' @rdname dot-conditionXlab
#' @param text.size size to fall back to
#' @keywords internal
.conditionTextSize = function(layout, text.size) {
    if (is.null(layout$size)) text.size else layout$size
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
    condition.layout = NULL
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
    
    profile_plot = ggplot(data = input, aes(x = .data$RUN, y = .data$newABUNDANCE,
                                            color = .data[[type_color]], linetype = .data$FEATURE)) +
        facet_grid(~LABEL) +
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
        scale_x_continuous(.conditionXlab(condition.layout), breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey", linetype = "longdash") +
        labs(title = unique(input$PROTEIN)) +
        geom_text(data = groupNametemp, aes(x = .data$RUN, y = .data$ABUNDANCE, label = .data$Label), 
                  size = .conditionTextSize(condition.layout, text.size), 
                  angle = text.angle, 
                  vjust = 1,
                  color = "black") +
        theme_msstats("PROFILEPLOT", x.axis.size, y.axis.size, legend.size)
    
    if (featureName == "TRANSITION") {
        color_guide = guide_legend(order=1,
                                   override.aes = list(size=1.2,
                                                       linetype = ss),
                                   title = paste("# peptide:", nlevels(input$PEPTIDE)), 
                                   title.theme = element_text(size = 13, angle = 0),
                                   keywidth = 0.25,
                                   keyheight = 0.1,
                                   default.unit = 'inch',
                                   ncol = 3)
    } else if (featureName == "PEPTIDE") {
        color_guide = guide_legend(order=1,
                                   title = paste("# peptide:", nlevels(input$PEPTIDE)), 
                                   title.theme = element_text(size = 13, angle = 0),
                                   keywidth = 0.25,
                                   keyheight = 0.1,
                                   default.unit = 'inch',
                                   ncol = 3)
    }
    shape_guide = guide_legend(order=2,
                               title = NULL,
                               label.theme = element_text(size = 11, angle = 0),
                               keywidth = 0.1,
                               keyheight = 0.1,
                               default.unit = 'inch')
    if (is_censored) {
        if (featureName == "NA") {
            profile_plot = profile_plot + guides(color = FALSE,
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
    yaxis.name, lineNameAxis, groupNametemp, condition.layout = NULL
) {
    RUN = ABUNDANCE = Name = NULL
    
    num_features = data.table::uniqueN(input$FEATURE)
    profile_plot = ggplot(data = input, 
                          aes(x = .data$RUN, y = .data$newABUNDANCE, 
                                     color = .data$analysis, linetype = .data$FEATURE, 
                                     size = .data$analysis)) +
        facet_grid(~LABEL) +
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
        scale_x_continuous(.conditionXlab(condition.layout), breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(xintercept = lineNameAxis + 0.5, 
                   colour = "grey", linetype = "longdash") +
        labs(title = unique(input$PROTEIN)) +
        geom_text(data = groupNametemp, aes(x = .data$RUN, y = .data$ABUNDANCE, label = .data$Label), 
                  size = .conditionTextSize(condition.layout, text.size), 
                  angle = text.angle, 
                  vjust = 1,
                  color = "black") +
        theme_msstats("PROFILEPLOT", x.axis.size, y.axis.size, 
                      legend.size, legend.title = element_blank())
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
    lineNameAxis, yaxis.name, condition.layout = NULL
) { 
    RUN = ABUNDANCE = Name = NULL
    
    if (all_proteins) {
        plot_title = "All"
    } else {
        plot_title = unique(input$PROTEIN)
    }
    
    ggplot(input, aes(x = .data$RUN, y = .data$ABUNDANCE)) +
        facet_grid(~LABEL) +
        geom_boxplot(aes(fill = .data$LABEL), outlier.shape = 1,
                     outlier.size = 1.5) +
        scale_fill_manual(values = label.color, guide = "none") +
        scale_x_discrete(.conditionXlab(condition.layout), breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey",
                   linetype = "longdash") +
        labs(title  =  plot_title) +
        geom_text(data = groupName, aes(x = .data$RUN, y = .data$ABUNDANCE, label = .data$Label),
                  size = .conditionTextSize(condition.layout, text.size),
                  angle = text.angle, vjust = 1, color = "black") +
        theme_msstats("QCPLOT", x.axis.size, y.axis.size,
                      legend_size = NULL)
    
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
