#' Check groupComparisonPlots parameters
#' @param type type of a plot: HEATMAP/VOLCANOPLOT/COMPARISONPLOT
#' @param log_base 2 or 10
#' @param selected_labels character vector of contrast labels
#' @param all_labels character vector of all contrast labels
#' @keywords internal
.checkGCPlotsInput = function(type, log_base, selected_labels, all_labels) {
    checkmate::assertChoice(type, c("HEATMAP", "VOLCANOPLOT", "COMPARISONPLOT"))
    checkmate::assertChoice(log_base, c(2, 10))
    if (selected_labels != "all") {
        if (is.character(selected_labels)) {
            chosen_labels = selected_labels
            print("labels")
            print(chosen_labels)
            wrong_labels = setdiff(chosen_labels, all_labels)
            if (length(wrong_labels) > 0) {
                msg_1 = paste("Please check labels of comparisons.",
                              "Result does not have the following comparisons:")
                msg_2 = paste(wrong_labels, sep = ", ", collapse = ", ")
                msg = paste(msg_1, msg_2)
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            }
        }
        if (is.numeric(selected_labels)) {
            n_labels = length(all_labels)
            if (n_labels < max(selected_labels)) {
                msg = paste("Please check your selection of comparisons. There are",
                            n_labels, "comparisons in this result.")
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            } else {
                chosen_labels = all_labels[selected_labels]
            }
        }  
    } else {
        chosen_labels = all_labels
    }
    chosen_labels
}


#' @importFrom stats quantile dist
#' @keywords internal
.getOrderedMatrix = function(input, type) {
    input_tmp = input
    input_tmp[is.na(input)] = 50
    if (toupper(type) == "PROTEIN") {
        input = input[hclust(dist(input_tmp), method = "ward.D")$order, ]
    } else if (toupper(type) == "COMPARISON") {
        input = input[, hclust(dist(t(input_tmp)), method = "ward.D")$order]
    } else if (toupper(type) == "BOTH") {
        input = input[hclust(dist(input_tmp), method = "ward.D")$order,
                      hclust(dist(t(input)), method = "ward.D")$order]
    }
    input
}


colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)

#' Create colorkey for ggplot2 heatmap
#' @param my.colors blocks
#' @keywords internal
.getColorKeyGGPlot2 = function(my.colors, blocks) {
    x.at = seq(-0.05, 1.05, length.out = 14)
    par(mar = c(3, 3, 3, 3), mfrow = c(3, 1), oma = c(3, 0, 3, 0))
    plot.new()
    image(z = matrix(seq(seq_len(length(my.colors) -1)), ncol = 1),
          col = my.colors,
          xaxt = "n",
          yaxt = "n")
    mtext("Color Key", side = 3,line = 1, cex = 3)
    mtext("(sign) Adjusted p-value", side = 1, line = 3, at = 0.5, cex = 1.7)
    mtext(blocks, side = 1, line = 1, at = x.at, cex = 1)
}

#' Create colorkey for plotly heatmap
#' @param my.colors blocks
#' @keywords internal
.getColorKeyPlotly = function(my.colors, blocks) {
    color.key.plot <- plotly::layout(
        plot_ly(type = "image", z = list(my.colors)),
        xaxis = list(
            dtick = 0,
            ticktext = as.character(blocks),
            tickmode = "array",
            tickvals = -0.5:length(blocks),
            tickangle = 0,
            title = "(sign) Adjusted p-value"
        ),
        yaxis = list(
            ticks = "",
            showticklabels = FALSE
        )
    )
    
    color.key.plot <- plotly::style(color.key.plot, hoverinfo = "none")
    color.key.plot
}

#' Create heatmap
#' @param input data.table
#' @inheritParams groupComparisonPlots
#' @keywords internal
.makeHeatmapPlotly = function(input, my.colors, my.breaks, x.axis.size, y.axis.size, height, numProtein) {
    input <- input[1:pmin(numProtein, nrow(input)), ,drop=F]
    par(oma = c(3, 0, 0, 4))
    label_formatter <- list(
        title = "",
        # titlefont = f1,
        showticklabels = TRUE,
        tickangle = 45,
        # tickfont = f2,
        exponentformat = "E")
    
    # adjust my.breaks
    x = my.breaks
    dltx <- diff(x)[1]
    x <- sort(c(x,-dltx/16,dltx/16))
    x <- x[x!=0]
    x.resc <- (x-min(x))/(max(x)-min(x))

    # get color scale
    cols = my.colors
    colorScale <- data.frame(
        z = c(0,rep(x.resc[2:(length(x.resc)-1)],each=2),1),
        col=rep(cols,each=2)
    )
    
    # Creating the custom hover text matrix
    row_names <- rownames(input)
    col_names <- colnames(input)
    hover_text_matrix <- matrix("", nrow = nrow(input), ncol = ncol(input))
    for (i in 1:nrow(input)) {
        for (j in 1:ncol(input)) {
            hover_text_matrix[i, j] <- sprintf("Comparison: %s<br>Protein: %s<br>Value: %.2f", 
                                               col_names[j], 
                                               row_names[i], 
                                               input[i, j])
        }
    }
    
    heatmap_plot = plot_ly(z = as.matrix(input),
                           zmin = x[1],
                           zmax = x[length(x)],
                           x = colnames(input),
                           xgap = 0,
                           y = rownames(input),
                           ygap = 0,
                           type = "heatmap",
                           hoverinfo = "text",
                           text=hover_text_matrix,
                           showlegend = FALSE, 
                           showscale = FALSE,
                           colorscale = colorScale,
                           colorbar = list(ypad = 520, tick0 = x[1], dtick = dltx, len = 1, orientation = "h"), 
                           width = 800)
    
    heatmap_plot <- plotly::layout(heatmap_plot, 
                                   xaxis = label_formatter, 
                                   plot_bgcolor = "grey",
                                   height = height)
    heatmap_plot
}

.makeHeatmapGgplot2 = function(input, my.colors, my.breaks, x.axis.size, y.axis.size,height) {
    par(oma = c(3, 0, 0, 4))
    heatmap.2(as.matrix(input),
              col = my.colors,
              Rowv = FALSE, Colv = FALSE,
              dendrogram = "none", breaks = my.breaks,
              trace = "none", na.color = "grey",
              cexCol = (x.axis.size / 10),
              cexRow = (y.axis.size / 10),
              key = FALSE,
              lhei = c(0.1, 0.9), lwid = c(0.1, 0.9))
}

#' Create a volcano plot
#' @inheritParams groupComparisonPlots
#' @param input data.table
#' @param label_name contrast label
#' @param log_base_FC 2 or 10
#' @param log_base_pval 2 or 10
#' @keywords internal
.makeVolcano = function(
    input, label_name, log_base_FC, log_base_pval, x.lim, ProteinName, dot.size,
    y.limdown, y.limup, text.size, FCcutoff, sig, x.axis.size, y.axis.size,
    legend.size, log_adjp
) {
    Protein = NULL
    plot = ggplot(aes_string(x = "logFC", 
                             y = log_adjp,
                             color = "colgroup",
                             label = "Protein"),
                  data = input) +
        geom_point(size = dot.size) +
        scale_colour_manual(values = c("gray65", "blue", "red"), 
                            limits = c("black", "blue", "red"), 
                            breaks = c("black", "blue", "red"), 
                            labels = c("No regulation", "Down-regulated", "Up-regulated")) +
        scale_y_continuous(paste0("-Log", log_base_pval, " (adjusted p-value)"),
                           limits = c(y.limdown, y.limup)) +
        labs(title = unique(label_name))
    plot = plot +
        scale_x_continuous(paste0("Log", log_base_pval, " fold change"),
                           limits = c(-x.lim, x.lim))
    if (ProteinName) {
        if (!(length(unique(input$colgroup)) == 1 & any(unique(input$colgroup) == "black"))) {
            plot = plot +
                geom_text_repel(data = input[input$colgroup != "black", ],
                                aes(label = Protein),
                                size = text.size,
                                col = "black")
        }
    } 
    if (!FCcutoff | is.numeric(FCcutoff)) {
        l = ifelse(!FCcutoff, 20, 10)
        sigcut = data.table::setnames(
            data.table::data.table("sigline", 
                                   seq(-x.lim, x.lim, length.out = l),
                                   (-log(sig, base = log_base_pval)),
                                   "twodash"), 
            c("Protein", "logFC", log_adjp, "line"))
    }
    if (!FCcutoff) {
        plot = plot +
            geom_line(data = sigcut,
                      aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                      colour = "darkgrey",
                      size = 0.6,
                      show.legend = TRUE) +
            scale_linetype_manual(values = c("twodash" = 6),
                                  labels = c(paste0("Adj p-value cutoff (", sig, ")"))) +
            guides(colour = guide_legend(override.aes = list(linetype = 0)),
                   linetype = guide_legend())
    }
    if (is.numeric(FCcutoff)) {
        FCcutpos = data.table::setnames(data.table("sigline", 
                                                   log(FCcutoff, log_base_pval), 
                                                   seq(y.limdown, y.limup, length.out = 10), 
                                                   "dotted"),
                                        c("Protein", "logFC", log_adjp, "line"))
        FCcutneg = data.table::setnames(data.table("sigline", 
                                                   (-log(FCcutoff, log_base_pval)), 
                                                   seq(y.limdown, y.limup, length.out = 10), 
                                                   "dotted"),
                                        c("Protein", "logFC", log_adjp, "line"))
        plot = plot +
            geom_line(data = sigcut, 
                      aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                      colour = "darkgrey",
                      size = 0.6,
                      show.legend = TRUE) +
            geom_line(data = FCcutpos,
                      aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                      colour = "darkgrey",
                      size = 0.6,
                      show.legend = TRUE) +
            geom_line(data = FCcutneg,
                      aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                      colour = "darkgrey",
                      size = 0.6) +
            scale_linetype_manual(values = c("dotted" = 3, "twodash" = 6),
                                  labels = c(paste0("Fold change cutoff (", FCcutoff, ")"),
                                             paste0("Adj p-value cutoff (", sig, ")"))) +
            guides(colour = guide_legend(override.aes = list(linetype = 0)),
                   linetype = guide_legend())
    }  
    plot = plot +
        theme_msstats("VOLCANOPLOT", x.axis.size, y.axis.size,
                      legend.size, strip_background = element_rect(),
                      strip_text_x = element_text(),
                      legend_position = "bottom", legend.title = element_blank())
    plot
}


#' Create comparison plot
#' @param input data.table
#' @param log_base 2 or 10
#' @inheritParams groupComparisonPlots
#' @keywords internal
.makeComparison = function(
    input, log_base, dot.size, x.axis.size, y.axis.size, 
    text.angle, hjust, vjust, y.limdown, y.limup
) {
    logFC = ciw = NULL
    
    protein = unique(input$Protein)
    plot = ggplot(input, aes_string(x = "Label", y = "logFC")) +
        geom_errorbar(aes(ymax = logFC + ciw, ymin = logFC - ciw),
                      data = input,
                      width = 0.1,
                      colour = "red") +
        geom_point(size = dot.size, 
                   colour = "darkred") +
        scale_x_discrete('Comparison') +
        geom_hline(yintercept = 0, 
                   linetype = "twodash", 
                   colour = "darkgrey", 
                   size = 0.6) +
        labs(title = protein) +
        theme_msstats("COMPARISONPLOT", x.axis.size, y.axis.size, 
                      text_angle = text.angle, text_hjust = hjust, 
                      text_vjust = vjust)
    plot = plot +
        scale_y_continuous(paste0("Log", log_base, "-Fold Change"),
                           limits = c(y.limdown, y.limup))
    plot
}
