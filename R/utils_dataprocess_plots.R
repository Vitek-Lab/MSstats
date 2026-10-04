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


#' Create profile plot
#' @inheritParams dataProcessPlots
#' @param input data.table
#' @param is_censored TRUE if censored values were imputed
#' @keywords internal
.makeProfilePlot = function(
    input, is_censored, featureName, y.limdown, y.limup, x.axis.size, 
    y.axis.size, text.size, text.angle, legend.size, dot.size.profile, 
    ss, s, cumGroupAxis, yaxis.name, lineNameAxis, groupNametemp, dot_colors
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
        scale_x_continuous('MS runs', breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey", linetype = "longdash") +
        labs(title = unique(input$PROTEIN)) +
        geom_text(data = groupNametemp, aes(x = .data$RUN, y = .data$ABUNDANCE, label = .data$Name), 
                  size = text.size, 
                  angle = text.angle, 
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
    yaxis.name, lineNameAxis, groupNametemp
) {
    RUN = ABUNDANCE = Name = NULL
    
    num_features = data.table::uniqueN(input$FEATURE)
    profile_plot = ggplot(data = input, 
                          aes(x = .data$RUN, y = .data$newABUNDANCE, 
                                     color = .data$analysis, linetype = .data$FEATURE)) +
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
        scale_x_continuous("MS runs", breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(xintercept = lineNameAxis + 0.5, 
                   colour = "grey", linetype = "longdash") +
        labs(title = unique(input$PROTEIN)) +
        geom_text(data = groupNametemp, aes(x = .data$RUN, y = .data$ABUNDANCE, label = .data$Name), 
                  size = text.size, 
                  angle = text.angle, 
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
#' @param isPlotly TRUE if the plot will be converted with ggplotly
#' @keywords internal
.makeQCPlot = function(
    input, all_proteins, y.limdown, y.limup, x.axis.size, y.axis.size, 
    text.size, text.angle, legend.size, label.color, cumGroupAxis, groupName,
    lineNameAxis, yaxis.name, isPlotly = FALSE
) { 
    RUN = ABUNDANCE = Name = NULL
    LABEL = x = lower = upper = ymin = ymax = NULL
    
    if (all_proteins) {
        plot_title = "All"
    } else {
        plot_title = unique(input$PROTEIN)
    }
    
    if (isPlotly) {
        # ggplotly converts geom_boxplot into a single native plotly box trace,
        # so the per-box drawing cost below does not apply here
        qc_plot = ggplot(input, aes(x = .data$RUN, y = .data$ABUNDANCE)) +
        facet_grid(~LABEL) +
        geom_boxplot(aes(fill = .data$LABEL), outlier.shape = 1,
                     outlier.size = 1.5) +
            scale_x_discrete("MS runs", breaks = cumGroupAxis)
    } else {
        # geom_boxplot draws every box (one per run) as a separate grob tree,
        # which grows linearly with the number of runs. Draw the same
        # statistics as a few vectorized layers instead.
        box_stats = .qcBoxStats(input, y.limdown, y.limup)
        boxes = box_stats$boxes
        whiskers = rbind(boxes[, list(LABEL, x, y = upper, yend = ymax)],
                         boxes[, list(LABEL, x, y = lower, yend = ymin)])
        break_positions = match(as.character(cumGroupAxis), box_stats$x_levels)
        has_break = !is.na(break_positions)
        qc_plot = ggplot(boxes) +
            facet_grid(~LABEL) +
            geom_segment(data = whiskers,
                         aes(x = .data$x, xend = .data$x,
                             y = .data$y, yend = .data$yend),
                         colour = "#333333", linewidth = 0.5) +
            geom_rect(aes(xmin = .data$x - 0.375, xmax = .data$x + 0.375,
                          ymin = .data$lower, ymax = .data$upper,
                          fill = .data$LABEL),
                      colour = "#333333", linewidth = 0.5) +
            geom_segment(aes(x = .data$x - 0.375, xend = .data$x + 0.375,
                             y = .data$middle, yend = .data$middle),
                         colour = "#333333", linewidth = 1) +
            geom_point(data = box_stats$outliers,
                       aes(x = .data$x, y = .data$ABUNDANCE),
                       shape = 1, size = 1.5, colour = "#333333") +
            # expansion matches the discrete axis: 0.4 to n + 0.6
            scale_x_continuous("MS runs",
                               breaks = break_positions[has_break],
                               labels = as.character(cumGroupAxis)[has_break],
                               expand = expansion(add = 0.225))
    }

    qc_plot +
        scale_fill_manual(values = label.color, guide = "none") +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey",
                   linetype = "longdash") +
        labs(title  =  plot_title) +
        geom_text(data = groupName, aes(x = .data$RUN, y = .data$ABUNDANCE, label = .data$Name),
                  size = text.size, angle = text.angle, color = "black") +
        theme_msstats("QCPLOT", x.axis.size, y.axis.size,
                      legend_size = NULL)
    
}


#' Boxplot statistics for the QC plot, matching ggplot2::stat_boxplot
#' @param input data.table with RUN (factor), LABEL and ABUNDANCE columns
#' @param y.limdown,y.limup y-axis limits; values outside them are dropped
#' before the statistics are computed, as scale_y_continuous(limits) does
#' @return list with `boxes` (one row per LABEL and run), `outliers` and
#' `x_levels` (the RUN levels present, in x-axis order)
#' @keywords internal
.qcBoxStats = function(input, y.limdown, y.limup) {
    RUN = LABEL = ABUNDANCE = x = lower = upper = NULL

    # discrete axis positions: index among RUN levels present in the data
    x_levels = levels(droplevels(input$RUN))
    values = data.table::data.table(
        LABEL = input$LABEL,
        x = match(as.character(input$RUN), x_levels),
        ABUNDANCE = input$ABUNDANCE)
    values = values[!is.na(ABUNDANCE) & ABUNDANCE >= y.limdown &
                        ABUNDANCE <= y.limup]
    boxes = values[, {
        q = as.numeric(stats::quantile(ABUNDANCE, c(0, 0.25, 0.5, 0.75, 1)))
        iqr = q[4] - q[2]
        is_outlier = ABUNDANCE < q[2] - 1.5 * iqr | ABUNDANCE > q[4] + 1.5 * iqr
        if (any(is_outlier)) {
            q[c(1, 5)] = range(c(q[2:4], ABUNDANCE[!is_outlier]))
        }
        list(ymin = q[1], lower = q[2], middle = q[3], upper = q[4],
             ymax = q[5])
    }, by = c("LABEL", "x")]
    outliers = values[boxes, on = c("LABEL", "x")][
        ABUNDANCE < lower - 1.5 * (upper - lower) |
            ABUNDANCE > upper + 1.5 * (upper - lower),
        list(LABEL, x, ABUNDANCE)]
    list(boxes = boxes, outliers = outliers, x_levels = x_levels)
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
