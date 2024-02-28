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
    
    profile_plot = ggplot(input, aes_string(x = "RUN", y = "newABUNDANCE",
                                            color = type_color, linetype = "FEATURE")) +
        facet_grid(~LABEL) +
        geom_line(size = 0.5)
    
    if (is_censored) {
        profile_plot = profile_plot +
        geom_point(aes_string(x = "RUN", y = "newABUNDANCE", color = type_color, shape = "censored"),
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
        geom_text(data = groupNametemp, aes(x = RUN, y = ABUNDANCE, label = Name), 
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
                          aes_string(x = "RUN", y = "newABUNDANCE", 
                                     color = "analysis", linetype = "FEATURE", 
                                     size = "analysis")) +
        facet_grid(~LABEL) +
        geom_line(size = 0.5)
    
    if (is_censored) { # splitting into two layers to keep red above grey
        profile_plot = profile_plot +
            geom_point(data = input[input$PEPTIDE != "Run summary"], 
                       aes_string(x = "RUN", y = "newABUNDANCE", 
                                  color = "analysis", size = "analysis", 
                                  shape = "censored")) +
            geom_point(data = input[input$PEPTIDE == "Run summary"], 
                       aes_string(x = "RUN", y = "newABUNDANCE", 
                                  color = "analysis", size = "analysis", 
                                  shape = "censored")) +
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
        geom_text(data = groupNametemp, aes(x = RUN, y = ABUNDANCE, label = Name), 
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
            geom_point(aes_string(x = "RUN", y = "newABUNDANCE", size = "analysis",
                                  color = "analysis"), data = input)
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
    lineNameAxis, yaxis.name
) { 
    RUN = ABUNDANCE = Name = NULL
    
    if (all_proteins) {
        plot_title = "All proteins"
    } else {
        plot_title = unique(input$PROTEIN)
    }
    
    ggplot(input, aes_string(x = "RUN", y = "ABUNDANCE")) +
        facet_grid(~LABEL) +
        geom_boxplot(aes_string(fill = "LABEL"), outlier.shape = 1,
                     outlier.size = 1.5) +
        scale_fill_manual(values = label.color, guide = "none") +
        scale_x_discrete("MS runs", breaks = cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey",
                   linetype = "longdash") +
        labs(title  =  plot_title) +
        geom_text(data = groupName, aes(x = RUN, y = ABUNDANCE, label = Name),
                  size = text.size, angle = text.angle, color = "black") +
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
    
    plot = ggplot(aes_string(x = "Label", y = "Mean"), data = input) +
        geom_errorbar(aes(ymax = Mean + ciw, ymin = Mean - ciw),
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
                   colour = "darkgrey", size = 0.6) +
        labs(title = unique(single_protein$PROTEIN)) +
        theme_msstats("CONDITIONPLOT", x.axis.size, y.axis.size, 
                      text_angle = text.angle)
    plot
}
