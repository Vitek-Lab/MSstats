.checkGCPlotsInput = function(type, log_base, selected_labels, all_labels) {
    checkmate::assertChoice(type, c("HEATMAP", "VOLCANOPLOT", "COMPARISONPLOT"))
    checkmate::assertChoice(log_base, c(2, 10))
    if (selected_labels != "all") {
        if (is.character(selected_labels)) {
            chosen_labels = selected_labels
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


.makeHeatmap = function(input, my.colors, my.breaks, 
                        x.axis.size, y.axis.size) {
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

.makeVolcano = function(
    input, label_name, log_base_FC, log_base_pval, x.lim, ProteinName, dot.size,
    y.limdown, y.limup, text.size, FCcutoff, sig, x.axis.size, y.axis.size,
    legend.size, log_adjp
) {
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



.makeComparison = function(
    input, log_base, dot.size, x.axis.size, y.axis.size, 
    text.angle, hjust, vjust, y.limdown, y.limup
) {
    protein = unique(input$Protein)
    plot = ggplot(input, aes_string(x = 'Label', y = 'logFC')) +
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
