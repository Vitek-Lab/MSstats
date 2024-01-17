#' Visualization for model-based analysis and summarizing differentially abundant proteins
#' 
#' @description To summarize the results of log-fold changes and adjusted p-values for differentially abundant proteins, 
#' groupComparisonPlots takes testing results from function (\code{\link{groupComparison}}) as input and 
#' automatically generate three types of figures in pdf files as output : 
#' (1) volcano plot (specify "VolcanoPlot" in option type) for each comparison separately; 
#' (2) heatmap (specify "Heatmap" in option type) for multiple comparisons ; 
#' (3) comparison plot (specify "ComparisonPlot" in option type) for multiple comparisons per protein.
#' 
#' @param data 'ComparisonResult' in testing output from function groupComparison.
#' @param type choice of visualization. "VolcanoPlot" represents volcano plot of log fold changes and adjusted p-values for each comparison separately. "Heatmap" represents heatmap of adjusted p-values for multiple comparisons. "ComparisonPlot" represents comparison plot of log fold changes for multiple comparisons per protein.
#' @param sig FDR cutoff for the adjusted p-values in heatmap and volcano plot. level of significance for comparison plot. 100(1-sig)\% confidence interval will be drawn.  sig=0.05 is default.
#' @param FCcutoff for volcano plot or heatmap, whether involve fold change cutoff or not. FALSE (default) means no fold change cutoff is applied for significance analysis. FCcutoff = specific value means specific fold change cutoff is applied.
#' @param logBase.pvalue for volcano plot or heatmap, (-) logarithm transformation of adjusted p-value with base 2 or 10(default).
#' @param ylimUp for all three plots, upper limit for y-axis. FALSE (default) for volcano plot/heatmap use maximum of -log2 (adjusted p-value) or -log10 (adjusted p-value). FALSE (default) for comparison plot uses maximum of log-fold change + CI.
#' @param ylimDown for all three plots, lower limit for y-axis. FALSE (default) for volcano plot/heatmap use minimum of -log2 (adjusted p-value) or -log10 (adjusted p-value). FALSE (default) for comparison plot uses minimum of log-fold change - CI.
#' @param xlimUp for Volcano plot, the limit for x-axis. FALSE (default) for  use maximum for absolute value of log-fold change or 3 as default if maximum for absolute value of log-fold change is less than 3.
#' @param x.axis.size size of axes labels, e.g. name of the comparisons in heatmap, and in comparison plot. Default is 10.
#' @param y.axis.size size of axes labels, e.g. name of targeted proteins in heatmap. Default is 10.
#' @param dot.size size of dots in volcano plot and comparison plot. Default is 3.
#' @param text.size size of ProteinName label in the graph for Volcano Plot. Default is 4.
#' @param text.angle angle of x-axis labels represented each comparison at the bottom of graph in comparison plot. Default is 0.
#' @param legend.size size of legend for color at the bottom of volcano plot.  Default is 7.
#' @param ProteinName for volcano plot only, whether display protein names or not. TRUE (default) means protein names, which are significant, are displayed next to the points. FALSE means no protein names are displayed.
#' @param colorkey TRUE(default) shows colorkey.
#' @param numProtein For ggplot2: The number of proteins which will be presented in each heatmap. Default is 100. Maximum possible number of protein for one heatmap is 180.
#' For Plotly: use this parameter to adjust the number of proteins to be displayed on the heatmap 
#' @param clustering Determines how to order proteins and comparisons. Hierarchical cluster analysis with Ward method(minimum variance) is performed. 'protein' means that protein dendrogram is computed and reordered based on protein means (the order of row is changed). 'comparison' means comparison dendrogram is computed and reordered based on comparison means (the order of comparison is changed). 'both' means to reorder both protein and comparison. Default is 'protein'.
#' @param width width of the saved file. Default is 10.
#' @param height height of the saved file. Default is 10.
#' @param which.Comparison list of comparisons to draw plots. List can be labels of comparisons or order numbers of comparisons from levels(data$Label), such as levels(testResultMultiComparisons$ComparisonResult$Label). Default is "all", which generates all plots for each protein.
#' @param which.Protein Protein list to draw comparison plots. List can be names of Proteins or order numbers of Proteins from levels(testResultMultiComparisons$ComparisonResult$Protein). Default is "all", which generates all comparison plots for each protein.
#' @param address the name of folder that will store the results. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory. An output pdf file is automatically created with the default name of "VolcanoPlot.pdf" or "Heatmap.pdf" or "ComparisonPlot.pdf". The command address can help to specify where to store the file as well as how to modify the beginning of the file name. If address=FALSE, plot will be not saved as pdf file but showed in window.
#' @param isPlotly This parameter is for MSstatsShiny application for plotly 
#' render, this cannot be used for saving PDF files as plotly do not have 
#' suppprt for PDFs currently. address and isPlotly cannot be set as TRUE at the
#' same time.
#' 
#' 
#' @details 
#' \itemize{
#' \item{Volcano plot : illustrate actual log-fold changes and adjusted p-values for each comparison separately with all proteins. The x-axis is the log fold change. The base of logarithm transformation is the same as specified in "logTrans" from \code{\link{dataProcess}}. The y-axis is the negative log2 or log10 adjusted p-values. The horizontal dashed line represents the FDR cutoff. The points below the FDR cutoff line are non-significantly abundant proteins (colored in black). The points above the FDR cutoff line are significantly abundant proteins (colored in red/blue for up-/down-regulated). If fold change cutoff is specified (FCcutoff = specific value), the points above the FDR cutoff line but within the FC cutoff line are non-significantly abundant proteins (colored in black)/}
#' \item{Heatmap : illustrate up-/down-regulated proteins for multiple comparisons with all proteins. Each column represents each comparison of interest. Each row represents each protein. Color red/blue represents proteins in that specific comparison are significantly up-regulated/down-regulated proteins with FDR cutoff and/or FC cutoff. The color scheme shows the evidences of significance. The darker color it is, the stronger evidence of significance it has. Color gold represents proteins are not significantly different in abundance.}
#' \item{Comparison plot : illustrate log-fold change and its variation of multiple comparisons for single protein. X-axis is comparison of interest. Y-axis is the log fold change. The red points are the estimated log fold change from the model. The blue error bars are the confidence interval with 0.95 significant level for log fold change. This interval is only based on the standard error, which is estimated from the model. }
#' }
#' 
#' @importFrom gplots heatmap.2
#' @importFrom stats hclust
#' @importFrom ggrepel geom_text_repel
#' @importFrom marray maPalette
#' @importFrom plotly ggplotly style add_trace plot_ly subplot
#' 
#' @export
#' 
#' @examples
#' QuantData<-dataProcess(SRMRawData, use_log_file = FALSE)
#' head(QuantData$FeatureLevelData)
#' ## based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
#' comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
#' comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
#' comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
#' comparison<-rbind(comparison1,comparison2, comparison3)
#' row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")
#' groups = levels(QuantData$ProteinLevelData$GROUP)
#' colnames(comparison) <- groups[order(as.numeric(groups))]
#' testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,
#' data=QuantData, 
#' use_log_file = FALSE)
#' testResultMultiComparisons$ComparisonResult
#' # Volcano plot with FDR cutoff = 0.05 and no FC cutoff
#' groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="VolcanoPlot",
#' logBase.pvalue=2, address="Ex1_")
#' # Volcano plot with FDR cutoff = 0.05, FC cutoff = 70, upper y-axis limit = 100, 
#' # and no protein name displayed
#' # FCcutoff=70 is for demonstration purpose
#' groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="VolcanoPlot",
#' FCcutoff=70, logBase.pvalue=2, ylimUp=100, ProteinName=FALSE,address="Ex2_")
#' # Heatmap with FDR cutoff = 0.05
#' groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="Heatmap", 
#' logBase.pvalue=2, address="Ex1_")
#' # Heatmap with FDR cutoff = 0.05 and FC cutoff = 70
#' # FCcutoff=70 is for demonstration purpose
#' groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="Heatmap",
#' FCcutoff=70, logBase.pvalue=2, address="Ex2_")
#' # Comparison Plot
#' groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="ComparisonPlot",
#' address="Ex1_")
#' # Comparison Plot
#' groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="ComparisonPlot",
#' ylimUp=8, ylimDown=-1, address="Ex2_")
#' 
groupComparisonPlots = function(
    data, type, sig = 0.05, FCcutoff = FALSE, logBase.pvalue = 10, ylimUp = FALSE,
    ylimDown = FALSE, xlimUp = FALSE, x.axis.size = 10, y.axis.size = 10, 
    dot.size = 3, text.size = 4, text.angle = 0, legend.size = 13, 
    ProteinName = TRUE, colorkey = TRUE, numProtein = 100, clustering = "both", 
    width = 800, height = 600, which.Comparison = "all", which.Protein = "all",
    address = "", isPlotly=FALSE
) {
    Label = Protein = NULL
    
    type = toupper(type)
    input = data.table::as.data.table(data)
    all_labels = as.character(unique(data$Label))
    log_base_FC = ifelse(is.element("log2FC", colnames(data)), 2, 10)
    
    getOption("MSstatsLog")("INFO", "MSstats - groupComparisonPlots function")
    chosen_labels = .checkGCPlotsInput(type, logBase.pvalue, which.Comparison,
                                       all_labels)
    input = input[Label %in% chosen_labels]
    input[, Protein := factor(Protein)]
    input[, Label := factor(Label)]
    warning("Avoid plotting all proteins as it can take a large amount of time 
            to download the files")
    if(isPlotly & address != FALSE) {
        print("Plots will be saved as .HTML file as plotly is selected, set isPlotly = FALSE, if 
            you want to generate PDF using ggplot2")
    }
    
    if (type == "HEATMAP") { 
        plotly_plot <- .plotHeatmap(input, logBase.pvalue, ylimUp, FCcutoff, sig, clustering, 
                     numProtein, colorkey, width, height, log_base_FC,
                     x.axis.size, y.axis.size, address, isPlotly)
        if(isPlotly) {
            if(address != FALSE) {
                .savePlotlyPlotHTML(list(plotly_plot),address,"Heatmap" ,width, height)
            }
            plotly_plot
        }
    }
    else if (type == "VOLCANOPLOT") {
        plots <- .plotVolcano(input, which.Comparison, address, width, height, logBase.pvalue,
                     ylimUp, ylimDown, FCcutoff, sig, xlimUp, ProteinName, dot.size,
                     text.size, legend.size, x.axis.size, y.axis.size, log_base_FC, isPlotly)
        plotly_plots <- vector("list", length(plots))
        if(isPlotly) {
            for(i in seq_along(plots)) {
                plot <- plots[[i]]
                plotly_plot <- .convertGgplot2Plotly(plot,tips=c("Protein","logFC","log10adjp","log2adjp"))
                plotly_plot <- .fixLegendPlotlyPlotsVolcano(plotly_plot)
                plotly_plots[[i]] = list(plotly_plot)
            }
            if(address != FALSE) {
                .savePlotlyPlotHTML(plotly_plots,address,"VolcanoPlot" ,width, height)
            }
            plotly_plots <- unlist(plotly_plots, recursive = FALSE)
            plotly_plots
        }
    }
    else if (type == "COMPARISONPLOT") {
        plots <- .plotComparison(input, which.Protein, address, width, height, sig, ylimUp, 
                        ylimDown, text.angle, dot.size, x.axis.size, y.axis.size,
                        log_base_FC, isPlotly)
        plotly_plots <- vector("list", length(plots))
        if(isPlotly) {
            for(i in seq_along(plots)) {
                plot <- plots[[i]]
                plotly_plot <- .convertGgplot2Plotly(plot,tips=c("logFC"))
                plotly_plots[[i]] = list(plotly_plot)
            }
            if(address != FALSE) {
                .savePlotlyPlotHTML(plotly_plots,address,"ComparisonPlot" ,width, height)
            }
            plotly_plots <- unlist(plotly_plots, recursive = FALSE)
            plotly_plots
        }
    }
}


#' Prepare data for heatmaps and plot them
#' @inheritParams groupComparisonPlots
#' @param input data.table
#' @param log_base_pval log base for p-values
#' @param log_base_FC log base for log-fold changes - 2 or 10
#' @keywords internal
.plotHeatmap = function(
    input, log_base_pval, ylimUp, FCcutoff, sig, clustering, numProtein, colorkey, 
    width, height, log_base_FC, x.axis.size, y.axis.size, address, isPlotly
) {
    adj.pvalue = heat_val = NULL
    
    if (length(unique(input$Protein)) <= 1) {
        stop("At least two proteins are needed for heatmaps.")
    }
    if (length(unique(input$Label)) <= 1) {
        stop("At least two comparisons are needed for heatmaps.")
    }
    
    if (is.numeric(ylimUp)) {
        y.limUp = ylimUp 
    } else {
        y.limUp = ifelse(log_base_pval == 2, 30, 10)
        input[adj.pvalue < log_base_pval ^ (-y.limUp), adj.pvalue := log_base_pval ^ (-y.limUp)]
    }  
    
    if (is.numeric(FCcutoff)) {
        input$adj.pvalue = ifelse(input[, 3] < log(FCcutoff, log_base_FC) & input[, 3] > -log(FCcutoff, log_base_FC), 
                                  1, input$adj.pvalue)
    }
    
    input[, heat_val := -log(adj.pvalue, log_base_pval) * sign(input[, 3])]
    wide = data.table::dcast(input, Protein ~ Label,
                             value.var = "heat_val")
    proteins = wide$Protein
    wide = as.matrix(wide[, -1])
    rownames(wide) = proteins
    wide = wide[rowSums(!is.na(wide)) != 0, colSums(!is.na(wide)) != 0]
    wide = .getOrderedMatrix(wide, clustering)
    up = 10 
    temp = 10 ^ (-sort(ceiling(seq(2, up, length = 10)[c(1, 2, 3, 5, 10)]), decreasing = TRUE))
    breaks = c(temp, sig)
    neg.breaks = log(breaks, log_base_pval)
    my.breaks = c(neg.breaks, 0, -neg.breaks[6:1], 101)
    blocks = c(-breaks, 1, breaks[6:1])
    namepro = rownames(wide)
    totalpro = length(namepro)
    numheatmap = totalpro %/% numProtein + 1
    
    my.colors = maPalette(low = "blue", high = "red", mid = "black", k = 12)
    my.colors = c(my.colors,"grey")
    blocks = c(blocks,"NA")
    my.colors.rgb = lapply(my.colors, function(x) as.vector(col2rgb(x)))   
    color.key.plotly = .getColorKeyPlotly(my.colors.rgb, blocks)

    if(!isPlotly) {
        if(colorkey) {
            .getColorKeyGGPlot2(my.colors, blocks)
        }
        savePlot(address, "Heatmap", width, height)
    }
    blue.red.18 = maPalette(low = "blue", high = "red", mid = "black", k = 14)
    if(isPlotly) {
        # If plotly, plot all proteins on a single heatmap
        heatmap  = .makeHeatmapPlotly(wide, blue.red.18, my.breaks, x.axis.size, y.axis.size, height, numProtein)
    } else {
        for (j in seq_len(numheatmap)) {
            if (j != numheatmap) {
                partial_wide = wide[((j - 1) * numProtein + 1):(j * numProtein), ]
            } else {
                partial_wide = wide[((j - 1) * numProtein + 1):nrow(wide), ]
            }
            heatmap  = .makeHeatmapGgplot2(partial_wide, my.colors, my.breaks, x.axis.size, y.axis.size, height)
        }
    }
    
    
    
    if (address != FALSE) {
        dev.off()
    }
    if(isPlotly) {
        if(colorkey) {
            heatmap_and_color_key <- subplot(heatmap, color.key.plotly, nrows = 2)
            
            heatmap_and_color_key <- plotly::layout(
                heatmap_and_color_key,
                annotations = list(
                    list(
                        x = 0.5,
                        y = 1.1,
                        text = "Heatmap",
                        showarrow = FALSE,
                        xref = "paper",
                        yref = "paper",
                        font = list(
                            size = 18
                        )
                    ),
                    list(
                        x = 0.5,
                        y = 0.35,
                        text = "Color Key",
                        showarrow = FALSE,
                        xref = "paper",
                        yref = "paper",
                        font = list(
                            size = 18
                        )
                    )
                ),
                margin = list(l = 50, r = 50, b = 50, t = 50)
            )
            heatmap_and_color_key
        }
        else {
            heatmap
        } 
    }
    
} 




#' Preprocess data for volcano plots and create them
#' @inheritParams groupComparisonPlots
#' @keywords internal
.plotVolcano = function(
    input, which.Comparison, address, width, height, log_base_pval,
    ylimUp, ylimDown, FCcutoff, sig, xlimUp, ProteinName, dot.size,
    text.size, legend.size, x.axis.size, y.axis.size, log_base_FC, isPlotly
) {
    adj.pvalue = colgroup = logFC = Protein = issue = Label = newlogFC = NULL
    
    log_adjp = paste0("log", log_base_pval, "adjp")
    all_labels = unique(input$Label)
    input = input[!is.na(adj.pvalue), ]
    colname_log_fc = intersect(colnames(input), c("log2FC", "log10FC"))
    data.table::setnames(input, colname_log_fc, c("logFC"))
    
    if (address == FALSE) {
        if (which.Comparison == "all") {
            if (length(unique(input$Label)) > 1) {
                stop('** Cannnot generate all volcano plots in a screen. Please set one comparison at a time.')
            }
        } else if (length(which.Comparison) > 1) {
            stop( '** Cannnot generate multiple volcano plots in a screen. Please set one comparison at a time.' )
        }
    }
    
    if (is.numeric(ylimUp)) {
        y.limUp = ylimUp 
    } else {
        y.limUp = ifelse(log_base_pval == 2, 30, 10)
    }
    input[, adj.pvalue := ifelse(adj.pvalue < log_base_pval ^ (-y.limUp),
                                 log_base_pval ^ (-y.limUp), adj.pvalue)]
    
    if (!FCcutoff) { 
        logFC_cutoff = 0
    } else {
        logFC_cutoff = log(FCcutoff, log_base_FC)
    }
    input[, colgroup := ifelse(adj.pvalue >= sig, "black",
                               ifelse(logFC > logFC_cutoff,
                                      "red", "blue"))]
    input[, colgroup := factor(colgroup, levels = c("black", "blue", "red"))]
    input[, Protein := as.character(Protein)]
    input[!is.na(issue) & issue == "oneConditionMissing", 
          Protein := paste0("*", Protein)]
    if(!isPlotly) {
        savePlot(address, "VolcanoPlot", width, height)
    }
    plots <- vector("list", length(all_labels))
    for (i in seq_along(all_labels)) {
        label_name = all_labels[i]
        single_label = input[Label == label_name, ]
        
        y.limup = ceiling(max(-log(single_label[!is.na(single_label$adj.pvalue), "adj.pvalue"], log_base_pval)))
        if (y.limup < (-log(sig, log_base_pval))) {
            y.limup = (-log(sig, log_base_pval) + 1) ## for too small y.lim
        }
        y.limdown = ifelse(is.numeric(ylimDown), ylimDown, 0)
        x_ceiling = ceiling(max(abs(single_label[!is.na(single_label$logFC) & is.finite(single_label$logFC), logFC])))
        x.lim = ifelse(is.numeric(xlimUp), xlimUp, ifelse((x_ceiling < 3), 3, x_ceiling))
        
        single_label[[log_adjp]] = -log(single_label$adj.pvalue, log_base_pval)
        single_label$newlogFC = single_label$logFC
        single_label[!is.na(issue) &
                         issue == "oneConditionMissing" & 
                         logFC == Inf, newlogFC := (x.lim - 0.2)]
        single_label[!is.na(issue) & 
                         issue == "oneConditionMissing" & 
                         logFC == (-Inf), newlogFC := (x.lim - 0.2) * (-1)]
        plot = .makeVolcano(single_label, label_name, log_base_FC, log_base_pval, x.lim, ProteinName, dot.size,
                            y.limdown, y.limup, text.size, FCcutoff, sig, x.axis.size, y.axis.size,
                            legend.size, log_adjp)
        print(plot)
        plots[[i]] = plot
    }
    if (address != FALSE) {
        dev.off()
    }
    if (isPlotly) {
        plots
    }
}

#' Preprocess data for comparison plots and create them
#' @inheritParams groupComparisonPlots
#' @param input data.table
#' @param log_base_FC log base for log-fold changes - 2 or 10
#' @keywords internal
.plotComparison = function(
    input, proteins, address, width, height, sig, ylimUp, ylimDown,
    text.angle, dot.size, x.axis.size, y.axis.size, log_base_FC, isPlotly
) {
    adj.pvalue = Protein = ciw = NULL
    
    input = input[!is.na(adj.pvalue), ]
    all_proteins = unique(input$Protein)
    
    if (address == FALSE) {
        if (proteins == "all" | length(proteins) > 1) {
            stop("** Cannnot generate all comparison plots in a screen. Please set one protein at a time.")
        }
    }
    if (proteins != "all") {
        selected_proteins = getSelectedProteins(proteins, all_proteins)
        input = input[Protein %in% selected_proteins, ]
    }
    
    all_proteins = unique(input$Protein)
    input$Protein = factor(input$Protein)
    if(!isPlotly) {
        savePlot(address, "ComparisonPlot", width, height)
    }
    plots <- vector("list", length(all_proteins))
    log_fc_column = intersect(colnames(input), c("log2FC", "log10FC"))
    for (i in seq_along(all_proteins)) {
        single_protein = input[Protein == all_proteins[i], ] 		
        single_protein[, ciw := qt(1 - sig / (2 * nrow(single_protein)), single_protein$DF) * single_protein$SE]
        data.table::setnames(single_protein, log_fc_column, "logFC")
        y.limup = ifelse(is.numeric(ylimUp), ylimUp, ceiling(max(single_protein$logFC + single_protein$ciw)))
        y.limdown = ifelse(is.numeric(ylimDown), ylimDown, floor(min(single_protein$logFC - single_protein$ciw)))
        hjust = ifelse(text.angle != 0, 1, 0.5)
        vjust = ifelse(text.angle != 0, 1, 0.5)
        
        plot = .makeComparison(single_protein, log_base_FC, dot.size, x.axis.size,
                               y.axis.size, text.angle, hjust, vjust, y.limdown, 
                               y.limup)
        print(plot)
        plots[[i]] = plot
    }
    if (address != FALSE) {
        dev.off()
    }
    if (isPlotly) {
        plots
    }
}
