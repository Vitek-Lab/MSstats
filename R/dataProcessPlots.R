#' Visualization for explanatory data analysis
#' 
#' @description To illustrate the quantitative data after data-preprocessing and 
#' quality control of MS runs, dataProcessPlots takes the quantitative data from 
#' function (\code{\link{dataProcess}}) as input and automatically generate 
#' three types of figures in pdf files as output : 
#' (1) profile plot (specify "ProfilePlot" in option type), 
#' to identify the potential sources of variation for each protein; 
#' (2) quality control plot (specify "QCPlot" in option type), 
#' to evaluate the systematic bias between MS runs; 
#' (3) mean plot for conditions (specify "ConditionPlot" in option type), 
#' to illustrate mean and variability of each condition per protein. 
#' 
#' @param data name of the (output of dataProcess function) data set.
#' @param type choice of visualization. "ProfilePlot" represents profile plot of 
#' log intensities across MS runs. "QCPlot" represents quality control plot of 
#' log intensities across MS runs. "ConditionPlot" represents mean plot of log 
#' ratios (Light/Heavy) across conditions.
#' @param featureName for "ProfilePlot" only, "Transition" (default) means 
#' printing feature legend in transition-level; "Peptide" means printing feature 
#' legend in peptide-level; "NA" means no feature legend printing.
#' @param ylimUp upper limit for y-axis in the log scale. FALSE(Default) for 
#' Profile Plot and QC Plot use the upper limit as rounded off maximum of 
#' log2(intensities) after normalization + 3. FALSE(Default) for Condition Plot
#' is maximum of log ratio + SD or CI.
#' @param ylimDown lower limit for y-axis in the log scale. FALSE(Default) for 
#' Profile Plot and QC Plot is 0. FALSE(Default) for Condition Plot is minumum
#' of log ratio - SD or CI.
#' @param scale for "ConditionPlot" only, FALSE(default) means each conditional 
#' level is not scaled at x-axis according to its actual value (equal space at 
#' x-axis). TRUE means each conditional level is scaled at x-axis according to
#' its actual value (unequal space at x-axis).
#' @param interval for "ConditionPlot" only, "CI"(default) uses confidence 
#' interval with 0.95 significant level for the width of error bar. 
#' "SD" uses standard deviation for the width of error bar.
#' @param x.axis.size size of x-axis labeling for "Run" in Profile Plot and 
#' QC Plot, and "Condition" in Condition Plot. Default is 10.
#' @param y.axis.size size of y-axis labels. Default is 10.
#' @param text.size size of labels represented each condition at the top of 
#' graph in Profile Plot and QC plot. Default is 4.
#' @param text.angle angle of labels represented each condition at the top
#' of graph in Profile Plot and QC plot or x-axis labeling in Condition plot. 
#' Default is 0.
#' @param legend.size size of feature legend (transition-level or peptide-level)
#' above graph in Profile Plot. Default is 7.
#' @param dot.size.profile size of dots in profile plot. Default is 2.
#' @param dot.size.condition size of dots in condition plot. Default is 3.
#' @param width width of the saved file. Default is 10.
#' @param height height of the saved file. Default is 10.
#' @param which.Protein Protein list to draw plots. List can be names of Proteins
#' or order numbers of Proteins from levels(data$FeatureLevelData$PROTEIN).
#' Default is "all", which generates all plots for each protein. 
#' For QC plot, "allonly" will generate one QC plot with all proteins.
#' @param originalPlot TRUE(default) draws original profile plots.
#' @param summaryPlot TRUE(default) draws profile plots with 
#' summarization for run levels.
#' @param save_condition_plot_result TRUE saves the table with values 
#' using condition plots. Default is FALSE.
#' @param remove_uninformative_feature_outlier It only works after users used 
#' featureSubset="highQuality" in dataProcess. TRUE allows to remove
#' 1) the features are flagged in the column, feature_quality="Uninformative" 
#' which are features with bad quality, 
#' 2) outliers that are flagged in the column, is_outlier=TRUE in Profile plots. 
#' FALSE (default) shows all features and intensities in profile plots.
#' @param address prefix for the filename that will store the results. 
#' @param isPlotly Parameter to use Plotly or ggplot2. If set to TRUE, MSstats 
#' will save Plotly plots as HTML files. If set to FALSE MSstats will save ggplot2 plots
#' as PDF files
#' Default folder is the current working directory. 
#' The other assigned folder has to be existed under the current working directory.
#'  An output pdf file is automatically created with the default name of 
#'  "ProfilePlot.pdf" or "QCplot.pdf" or "ConditionPlot.pdf" or "ConditionPlot_value.csv". 
#'  The command address can help to specify where to store the file as well as 
#'  how to modify the beginning of the file name. 
#'  If address=FALSE, plot will be not saved as pdf file but showed in window.
#' 
#' @details
#' \itemize{
#' \item{Profile Plot : identify the potential sources of variation of each protein. QuantData$FeatureLevelData is used for plots. X-axis is run. Y-axis is log-intensities of transitions. Reference/endogenous signals are in the left/right panel. Line colors indicate peptides and line types indicate transitions. In summarization plots, gray dots and lines are the same as original profile plots with QuantData$FeatureLevelData. Dark dots and lines are for summarized intensities from QuantData$ProteinLevelData.}
#' \item{QC Plot : illustrate the systematic bias between MS runs. After normalization, the reference signals for all proteins should be stable across MS runs. QuantData$FeatureLevelData is used for plots. X-axis is run. Y-axis is log-intensities of transition. Reference/endogenous signals are in the left/right panel. The pdf file contains (1) QC plot for all proteins and (2) QC plots for each protein separately.}
#' \item{Condition Plot : illustrate the systematic difference between conditions. Summarized intensnties from QuantData$ProteinLevelData are used for plots. X-axis is condition. Y-axis is summarized log transformed intensity. If scale is TRUE, the levels of conditions is scaled according to its actual values at x-axis. Red points indicate the mean for each condition. If interval is "CI", blue error bars indicate the confidence interval with 0.95 significant level for each condition. If interval is "SD", blue error bars indicate the standard deviation for each condition.The interval is not related with model-based analysis.}
#' }
#' The input of this function is the quantitative data from function \code{\link{dataProcess}}. 
#' 
#' 
#' @import ggplot2
#' @importFrom graphics axis image legend mtext par plot.new title plot
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom plotly ggplotly style add_trace plot_ly subplot
#' @importFrom htmltools tagList div save_html
#' 
#' @export
#' 
#' @examples 
#' # Consider quantitative data (i.e. QuantData) from a yeast study with ten time points of interests, 
#' # three biological replicates, and no technical replicates which is a time-course experiment. 
#' # The goal is to provide pre-analysis visualization by automatically generate two types of figures 
#' # in two separate pdf files. 
#' # Protein IDHC (gene name IDP2) is differentially expressed in time point 1 and time point 7, 
#' # whereas, Protein PMG2 (gene name GPM2) is not.
#' 
#' QuantData<-dataProcess(SRMRawData, use_log_file = FALSE)
#' head(QuantData$FeatureLevelData)
#' # Profile plot
#' dataProcessPlots(data=QuantData,type="ProfilePlot")
#' # Quality control plot 
#' dataProcessPlots(data=QuantData,type="QCPlot")	
#' # Quantification plot for conditions
#' dataProcessPlots(data=QuantData,type="ConditionPlot")
#' 
dataProcessPlots = function(
  data, type, featureName = "Transition", ylimUp = FALSE, ylimDown = FALSE,
  scale = FALSE, interval = "CI", x.axis.size = 10, y.axis.size = 10,
  text.size = 4, text.angle = 0, legend.size = 7, dot.size.profile = 2,
  dot.size.condition = 3, width = 800, height = 600, which.Protein = "all",
  originalPlot = TRUE, summaryPlot = TRUE, save_condition_plot_result = FALSE,
  remove_uninformative_feature_outlier = FALSE, address = "", isPlotly = FALSE
) {
  PROTEIN = Protein = NULL
  
  type = toupper(type)
  processed = data.table::as.data.table(data$FeatureLevelData)
  summarized = data.table::as.data.table(data$ProteinLevelData)
  processed[, PROTEIN := factor(PROTEIN)]
  summarized[, Protein := factor(Protein)]
  
  checkmate::assertChoice(type, c("PROFILEPLOT", "QCPLOT", "CONDITIONPLOT"),
                          .var.name = "type")
  if (as.character(address) == "FALSE") {
    if (which.Protein == "all") {
      stop("** Cannnot generate all plots in a screen. Please set one protein at a time.")
    } else if (length(which.Protein) > 1) {
      stop("** Cannnot generate multiple plots in a screen. Please set one protein at a time.")
    }
  }
  warning("Avoid plotting all proteins as it can take a large amount of time 
            to download the files")
  if(isPlotly & address != FALSE) {
      print("Plots will be saved as .HTML file as plotly is selected, set isPlotly = FALSE, if 
            you want to generate PDF using ggplot2")
  }
  
  if (type == "PROFILEPLOT") {
      plots <- .plotProfile(processed, summarized, featureName, ylimUp, ylimDown,
                           x.axis.size, y.axis.size, text.size, text.angle, legend.size, 
                           dot.size.profile, width, height, which.Protein, originalPlot, 
                           summaryPlot, remove_uninformative_feature_outlier, address, isPlotly)
      plotly_plots = list()
      if(isPlotly) {
          og_plotly_plot = NULL
          summ_plotly_plot = NULL
          if("original_plot" %in% names(plots)) {
              for(i in seq_along(plots[["original_plot"]])) {
                  plot_i <- plots[["original_plot"]][[paste("plot",i)]]
                  og_plotly_plot <- .convertGgplot2Plotly(plot_i,tips=c("FEATURE","RUN","newABUNDANCE"))
                  og_plotly_plot = .fixLegendPlotlyPlotsDataprocess(og_plotly_plot)
                  og_plotly_plot = .fixCensoredPointsLegendProfilePlotsPlotly(og_plotly_plot)

                  if(toupper(featureName) == "NA") {
                      og_plotly_plot = .retainCensoredDataPoints(og_plotly_plot)
                  }
                  plotly_plots = c(plotly_plots, list(og_plotly_plot))
              }
          }
          if("summary_plot" %in% names(plots)) {
              for(i in seq_along(plots[["summary_plot"]])) {
                  plot_i <- plots[["summary_plot"]][[paste("plot",i)]]
                  summ_plotly_plot <- .convertGgplot2Plotly(plot_i,tips=c("FEATURE","RUN","ABUNDANCE"))
                  summ_plotly_plot = .fixLegendPlotlyPlotsDataprocess(summ_plotly_plot)
                  summ_plotly_plot = .fixCensoredPointsLegendProfilePlotsPlotly(summ_plotly_plot)
                  if(toupper(featureName) == "NA") {
                      summ_plotly_plot = .retainCensoredDataPoints(summ_plotly_plot)
                  }
                  plotly_plots = c(plotly_plots, list(summ_plotly_plot))
              }
          }
          
          if(address != FALSE) {
              .savePlotlyPlotHTML(plotly_plots,address,"ProfilePlot" ,width, height)
          }
          plotly_plots
      }
  }
    
  else if (type == "QCPLOT") {
      plots <- .plotQC(processed, featureName, ylimUp, ylimDown, x.axis.size, y.axis.size, 
                      text.size, text.angle, legend.size, dot.size.profile, width, height,
                      which.Protein, address, isPlotly)
      plotly_plots <- vector("list", length(plots))
      if(isPlotly) {
          for(i in seq_along(plots)) {
              plot <- plots[[i]]
              plotly_plot <- .convertGgplot2Plotly(plot)
              plotly_plot = .fixLegendPlotlyPlotsDataprocess(plotly_plot)
              plotly_plots[[i]] = list(plotly_plot)
          }
            if(address != FALSE) {
                .savePlotlyPlotHTML(plotly_plots,address,"QCPlot" ,width, height)
            }
          plotly_plots <- unlist(plotly_plots, recursive = FALSE)
          plotly_plots
      }
  }
    
  else if (type == "CONDITIONPLOT") {
      plots <- .plotCondition(processed, summarized, ylimUp, ylimDown, scale, interval,
                             x.axis.size, y.axis.size, text.size, text.angle, legend.size, 
                             dot.size.profile, dot.size.condition, width, height,
                             which.Protein, save_condition_plot_result, address, isPlotly)
      plotly_plots <- vector("list", length(plots))
      if(isPlotly) {
          for(i in seq_along(plots)) {
              plot <- plots[[i]]
              plotly_plot <- .convertGgplot2Plotly(plot)
              plotly_plot = .fixLegendPlotlyPlotsDataprocess(plotly_plot)
              plotly_plots[[i]] = list(plotly_plot)
          }
          if(address != FALSE) {
              .savePlotlyPlotHTML(plotly_plots,address,"ConditionPlot" ,width, height)
          }
          plotly_plots <- unlist(plotly_plots, recursive = FALSE)
          plotly_plots
      }
  }
}


#' @importFrom utils setTxtProgressBar
#' @importFrom stats xtabs
#' @keywords internal
.plotProfile = function(
  processed, summarized, featureName, ylimUp, ylimDown, x.axis.size, y.axis.size, 
  text.size, text.angle, legend.size, dot.size.profile, width, height, proteins, 
  originalPlot, summaryPlot, remove_uninformative_feature_outlier, address, isPlotly
) {
  ABUNDANCE = PROTEIN = feature_quality = is_outlier = Protein = GROUP = NULL
  SUBJECT = LABEL = RUN = xtabs = PEPTIDE = FEATURE = NULL
  LogIntensities = TRANSITION = FRACTION = censored = analysis = NULL
  
  yaxis.name = .getYaxis(processed)
  is_censored = is.element("censored", colnames(processed))
  all_proteins = as.character(unique(processed$PROTEIN))
  if (remove_uninformative_feature_outlier) {
    if (is.element("feature_quality", colnames(processed))) {
      processed[, ABUNDANCE := ifelse(
        feature_quality == "Noninformative" | is_outlier, NA, ABUNDANCE)]
      msg = "** Filtered out uninformative feature and outliers in the profile plots."
    } else {
      msg = "** To remove uninformative features or outliers, please use \"featureSubset == \"highQuality\" option in \"dataProcess\" function."
    }
    getOption("MSstatsMsg")("INFO", msg)
  }
  
  processed = processed[order(GROUP, SUBJECT, LABEL)]
  processed[, RUN := factor(RUN, levels = unique(RUN), 
                            labels = seq(1, length(unique(RUN))))]
  processed[, RUN := as.numeric(RUN)]
  
  summarized = summarized[order(GROUP, SUBJECT)]
  summarized[, RUN := factor(RUN, levels = unique(RUN), 
                            labels = seq(1, length(unique(RUN))))]
  summarized[, RUN := as.numeric(RUN)]
  
  ## Meena :due to GROUP=0 for labeled.. extra care required.
  tempGroupName = unique(processed[, c("GROUP", "RUN")])
  if (length(unique(processed$LABEL)) == 2) {
    tempGroupName = tempGroupName[GROUP != '0']
  } 
  tempGroupName = tempGroupName[order(RUN), ] ## Meena : should we order by GROUP or RUN? I guess by RUn, because x-axis is by RUN
  level.group = as.character(unique(tempGroupName$GROUP))
  tempGroupName$GROUP = factor(tempGroupName$GROUP,
                                levels = level.group) ## Meena : factor GROUP again, due to 1, 10, 2, ... if you have better way, please change
  
  groupAxis = as.numeric(xtabs(~GROUP, tempGroupName))
  cumGroupAxis = cumsum(groupAxis)
  lineNameAxis = cumGroupAxis[-nlevels(tempGroupName$GROUP)]

  if (proteins != "all") {
      selected_proteins = getSelectedProteins(proteins, all_proteins)
      processed = processed[PROTEIN %in% selected_proteins]
      summarized = summarized[Protein %in% selected_proteins]
      processed[, PROTEIN := factor(PROTEIN)]
      summarized[, PROTEIN := factor(Protein)]
  }
  
  y.limup = ifelse(is.numeric(ylimUp), ylimUp, ceiling(max(processed$ABUNDANCE, na.rm = TRUE) + 3))
  y.limdown = ifelse(is.numeric(ylimDown), ylimDown, -1)
  
  groupName = data.frame(RUN = c(0, lineNameAxis) + groupAxis / 2 + 0.5,
                         ABUNDANCE = rep(y.limup - 1, length(groupAxis)),
                         Name = levels(tempGroupName$GROUP))

  
  if (length(unique(processed$LABEL)) == 2) {
    processed[, LABEL := factor(LABEL, labels = c("Reference", "Endogenous"))]
  } else {
    if (unique(processed$LABEL) == "L") {
      processed[, LABEL := factor(LABEL, labels = c("Endogenous"))]
    } else {
      processed[, LABEL := factor(LABEL, labels = c("Reference"))]
    }
  }
  
  if ("feature_quality" %in% colnames(processed)) {
    processed[, feature_quality := NULL]
  }
  if ("is_outlier" %in% colnames(processed)) {
    processed[, is_outlier := NULL]
  }
  output_plots <- list()
  output_plots[["original_plot"]] = list()
  output_plots[["summary_plot"]] = list()
  all_proteins = levels(processed$PROTEIN)
  if (originalPlot) {
    if(!isPlotly) {
        savePlot(address, "ProfilePlot", width, height)
    }
    pb = utils::txtProgressBar(min = 0, max = length(all_proteins), style = 3)
    for (i in seq_along(all_proteins)) {
      single_protein = .getSingleProteinForProfile(processed, all_proteins, i)
      if (all(is.na(single_protein$ABUNDANCE))) {
        next()
      }
      
      pept_feat = unique(single_protein[, list(PEPTIDE, FEATURE)])
      counts = pept_feat[, list(N = .N), by = "PEPTIDE"]$N
      s = rep(seq_along(counts), times = counts)
      ss = unlist(lapply(counts, function(x) seq(1, x)), FALSE, FALSE)
      groupNametemp = data.frame(groupName,
                                 "FEATURE" = unique(single_protein$FEATURE)[1],
                                 "PEPTIDE" = unique(single_protein$PEPTIDE)[1])
      
      dot_colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      check_length = length(unique(s)) %/% length(dot_colors)
      if (check_length > 0) {
        dot_colors = rep(dot_colors, times = check_length + 1)
      }
      
      profile_plot = .makeProfilePlot(single_protein, is_censored, featureName, 
                                      y.limdown, y.limup,
                                      x.axis.size, y.axis.size, 
                                      text.size, text.angle, 
                                      legend.size, dot.size.profile, 
                                      ss, s, cumGroupAxis, yaxis.name,
                                      lineNameAxis, groupNametemp, dot_colors)
      
      setTxtProgressBar(pb, i)
      print(profile_plot)
      output_plots[["original_plot"]][[paste("plot",i)]] <- profile_plot
    }
    
    close(pb)
    
    if (address != FALSE & !isPlotly) {
      dev.off()
    } 
  }
  
  if (summaryPlot) {
    protein_by_run = expand.grid(Protein = unique(summarized$Protein), 
                                 RUN = unique(summarized$RUN))
    summarized = merge(summarized, protein_by_run, by = c("Protein", "RUN"),
                       all.x = TRUE, all.y = TRUE)
    if(!isPlotly) {
        savePlot(address, "ProfilePlot_wSummarization", width, height)
    }
    pb = utils::txtProgressBar(min = 0, max = length(all_proteins), style = 3)
    for (i in seq_along(all_proteins)) {
      single_protein = .getSingleProteinForProfile(processed, all_proteins, i)
      if (all(is.na(single_protein$ABUNDANCE))) {
        next()
      }
      
      pept_feat = unique(single_protein[, list(PEPTIDE, FEATURE)])
      counts = pept_feat[, list(N = .N), by = "PEPTIDE"]$N
      s = rep(seq_along(counts), times = counts)
      ss = unlist(lapply(counts, function(x) seq(1, x)), FALSE, FALSE)
      groupNametemp = data.frame(groupName,
                                 FEATURE = unique(single_protein$FEATURE)[1],
                                 analysis = "Run summary")
      
      single_protein_summ = summarized[Protein == all_proteins[i], ]
      quant = single_protein_summ[
        Protein == all_proteins[i],
        list(PROTEIN = unique(Protein), PEPTIDE = "Run summary",
             TRANSITION = "Run summary", FEATURE = "Run summary",
             LABEL = "Endogenous", RUN = RUN,
             ABUNDANCE = LogIntensities, FRACTION = 1)
        ]
      if (is_censored) {
        quant$censored = FALSE
      }
      quant$analysis = "Run summary"
      quant$newABUNDANCE = quant$ABUNDANCE
      single_protein$analysis = "Processed feature-level data"
      combined = rbind(single_protein[
        , 
        list(PROTEIN, PEPTIDE, TRANSITION, FEATURE, LABEL,
             RUN, ABUNDANCE, newABUNDANCE,FRACTION, censored, analysis)], quant,fill=TRUE)
      combined$analysis = factor(combined$analysis)
      combined$FEATURE = factor(combined$FEATURE)
      combined$RUN = as.numeric(combined$RUN)
      profile_plot = .makeSummaryProfilePlot(
        combined, is_censored, y.limdown, y.limup, x.axis.size, y.axis.size, 
        text.size, text.angle, legend.size, dot.size.profile, cumGroupAxis, 
        yaxis.name, lineNameAxis, groupNametemp
      )
      print(profile_plot)
      setTxtProgressBar(pb, i)
      output_plots[["summary_plot"]][[paste("plot",i)]] <- profile_plot
      
    }
    close(pb)
    
    if (address != FALSE & !isPlotly) {
      dev.off()
    } 
  }
  if(isPlotly) {
      output_plots
  }
  
}


#' @importFrom stats xtabs
#' @importFrom utils setTxtProgressBar
.plotQC = function(
  processed, featureName, ylimUp, ylimDown, x.axis.size, y.axis.size, text.size, 
  text.angle, legend.size, dot.size.profile, width, height, protein, address, isPlotly
) {
  GROUP = SUBJECT = RUN = LABEL = PROTEIN = NULL
  
  yaxis.name = .getYaxis(processed)
  y.limup = ifelse(is.numeric(ylimUp), ylimUp, ceiling(max(processed$ABUNDANCE, na.rm = TRUE) + 3))
  y.limdown = ifelse(is.numeric(ylimDown), ylimDown, -1)
  
  processed = processed[order(GROUP, SUBJECT), ]
  processed[, RUN := factor(RUN, levels = unique(processed$RUN),
                            labels = seq(1, data.table::uniqueN(processed$RUN)))]
  
  if (length(unique(processed$LABEL)) == 2) {
    processed[, LABEL := factor(LABEL, labels = c("Reference", "Endogenous"))]
    label.color = c("darkseagreen1", "lightblue")
  } else {
    if (unique(processed$LABEL) == "L") {
      processed[, LABEL := factor(LABEL, labels = c("Endogenous"))]
      label.color = c("lightblue")
    } else {
      processed[, LABEL := factor(LABEL, labels = c("Reference"))]
      label.color = c("darkseagreen1")
    }
  }
  
  processed = processed[order(LABEL, GROUP, SUBJECT)]
  
  ## Meena :due to GROUP=0 for labeled.. extra care required.
  tempGroupName = unique(processed[, list(GROUP, RUN)])
  if (length(unique(processed$LABEL)) == 2) {
    tempGroupName = tempGroupName[GROUP != '0']
  } 
  tempGroupName = tempGroupName[order(RUN), ] ## Meena : should we order by GROUP or RUN? I guess by RUn, because x-axis is by RUN
  level.group = as.character(unique(tempGroupName$GROUP))
  tempGroupName$GROUP = factor(tempGroupName$GROUP,
                                levels = level.group) ## Meena : factor GROUP again, due to 1, 10, 2, ... if you have better way, please change
  
  groupAxis = as.numeric(xtabs(~GROUP, tempGroupName))
  cumGroupAxis = cumsum(groupAxis)
  lineNameAxis = cumGroupAxis[-nlevels(tempGroupName$GROUP)]
  groupName = data.frame(RUN = c(0, lineNameAxis) + groupAxis / 2 + 0.5,
                         ABUNDANCE = rep(y.limup - 1, length(groupAxis)),
                         Name = levels(tempGroupName$GROUP))
  if (!isPlotly) {
      savePlot(address, "QCPlot", width, height)
  }
  all_proteins = as.character(levels(processed$PROTEIN))
  plots <- vector("list", length(all_proteins) + 1) # +1 for all/allonly plot
  if (protein %in% c("all", "allonly")) {
    qc_plot = .makeQCPlot(processed, TRUE, y.limdown, y.limup, x.axis.size, 
                          y.axis.size, text.size, text.angle, legend.size, 
                          label.color, cumGroupAxis, groupName, lineNameAxis, 
                          yaxis.name)
    print(qc_plot)
    plots[[1]] = qc_plot
  } 
  
  if (protein != "allonly") {
    all_proteins = as.character(levels(processed$PROTEIN))
    if (protein != "all") {
      selected_proteins = getSelectedProteins(protein, all_proteins)
      processed = processed[PROTEIN %in% selected_proteins]
      processed[, PROTEIN := factor(PROTEIN)]
    }
    pb = utils::txtProgressBar(min = 0, max = length(all_proteins), style = 3)
    for (i in seq_along(all_proteins)) {	
    
      single_protein = processed[processed$PROTEIN == all_proteins[i], ]
      single_protein = single_protein[order(LABEL, RUN)]
      if (all(is.na(single_protein$ABUNDANCE))) {
        next()
      }
      qc_plot = .makeQCPlot(single_protein, FALSE, y.limdown, y.limup, 
                            x.axis.size, y.axis.size, text.size, text.angle, 
                            legend.size, label.color, cumGroupAxis, groupName,
                            lineNameAxis, yaxis.name)
      print(qc_plot)
      plots[[i+1]] = qc_plot # to accomodate all proteins
      setTxtProgressBar(pb, i)
    } 
    close(pb)
  } 
  if (address != FALSE) {
    dev.off()
  }
  if (isPlotly) {
      plots <- Filter(function(x) !is.null(x), plots) # remove if protein was not "all"
      plots
  }
} 


#' @importFrom stats qt sd na.omit
#' @importFrom utils setTxtProgressBar
#' @keywords internal
.plotCondition = function(
  processed, summarized, ylimUp, ylimDown, scale, interval, x.axis.size, 
  y.axis.size, text.size, text.angle, legend.size, dot.size.profile, 
  dot.size.condition, width, height, protein, save_plot, address, isPlotly
) {
  adj.pvalue = Protein = ciw = PROTEIN = GROUP = SUBJECT = ABUNDANCE = NULL
  
  data.table::setnames(summarized, c("Protein", "LogIntensities"),
                       c("PROTEIN", "ABUNDANCE"))
  all_proteins = levels(summarized$PROTEIN)
  if (protein != "all") {
    proteins = getSelectedProteins(protein, all_proteins)
    summarized = summarized[PROTEIN %in% proteins, ]
    summarized[, PROTEIN := factor(PROTEIN)]
  }
  yaxis.name = .getYaxis(processed)
  
  results = vector("list", length(all_proteins))
  if(!isPlotly) {
      savePlot(address, "ConditionPlot", width, height)
  }
  plots <- vector("list", length(all_proteins))
  pb = utils::txtProgressBar(min = 0, max = length(all_proteins), style = 3)
  for (i in seq_along(all_proteins)) {
    single_protein = summarized[PROTEIN == all_proteins[i], ]
    single_protein = na.omit(single_protein)
    single_protein[, GROUP := factor(GROUP)]
    single_protein[, SUBJECT := factor(SUBJECT)]
    if (all(is.na(single_protein$ABUNDANCE))) {
      next()
    }
    sp_all = single_protein[, list(Mean = mean(ABUNDANCE, na.rm = TRUE),
                                   SD = sd(ABUNDANCE, na.rm = TRUE),
                                   numMeasurement = .N),
                            by = "GROUP"]
    if (interval == "CI") {
      sp_all[, ciw := qt(0.975, sp_all$numMeasurement) * sp_all$SD / sqrt(sp_all$numMeasurement)]
    } else {
      sp_all[, ciw := sp_all$SD]
    }
    if (sum(is.na(sp_all$ciw)) >= 1) {
      sp_all[is.na(sp_all$ciw), ciw := 0]
    }
    
    y.limup = ifelse(is.numeric(ylimUp), ylimUp, ceiling(max(sp_all$Mean + sp_all$ciw)))
    y.limdown = ifelse(is.numeric(ylimDown), ylimDown, floor(min(sp_all$Mean - sp_all$ciw)))
    
    sp_all = sp_all[order(GROUP), ]
    sp_all$Protein = all_proteins[i]
    if (save_plot) {
      results[[i]] = sp_all
    }
    con_plot = .makeConditionPlot(sp_all, scale, single_protein, y.limdown, 
                                  y.limup, x.axis.size, y.axis.size, 
                                  text.size, text.angle, legend.size, 
                                  dot.size.condition, yaxis.name)
    print(con_plot)
    plots[[i]] = con_plot
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  if (address != FALSE & !isPlotly) {
    dev.off()
  }
  
  if (save_plot) {
    result = data.table::rbindlist(results)
    data.table::setnames(result, "GROUP", "Condition")
    if (interval == "CI") {
      data.table::setnames(result, "ciw", "95% CI")
    } else if (interval == "SD") {
      data.table::setnames(result, "ciw", "SD")
    }
    .saveTable(result, address, "ConditionPlot_value")
  }
  if (isPlotly) {
      plots
  }
}

#' converter for plots from ggplot to plotly
#' @noRd
.convertGgplot2Plotly = function(plot, tips = "all") {
    converted_plot <- ggplotly(plot,tooltip = tips)
    converted_plot <- plotly::layout(
            converted_plot,
            width = 800,   # Set the width of the chart in pixels
            height = 600,  # Set the height of the chart in pixels
            title = list(
                font = list(
                    size = 18
                )
            ),
            xaxis = list(
                titlefont = list(
                    size = 15  # Set the font size for the x-axis label
                )
            ),
            legend = list(
                x = 0,     # Set the x position of the legend
                y = -0.25,    # Set the y position of the legend (negative value to move below the plot)
                orientation = "h",  # Horizontal orientation
                font = list(
                    size = 12  # Set the font size for legend item labels
                ),
                title = list(
                    font = list(
                        size = 12  # Set the font size for the legend title
                    )
                )
            )
        ) 
    converted_plot
}

.retainCensoredDataPoints = function(plot) {
    df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
    for (i in seq_along(plot$x$data)) {
        if (df$legend_entries[i] != "Detected data" && df$legend_entries[i] != "Censored missing data") {
            plot$x$data[[i]]$showlegend <- FALSE
        }
    }
    plot
}

.fixLegendPlotlyPlotsDataprocess = function(plot) {
    df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
    df$legend_group <- gsub("^\\((.*?),.*", "\\1", df$legend_entries)
    df$is_first <- !duplicated(df$legend_group)
    df$is_bool <- ifelse(grepl("TRUE|FALSE", df$legend_group), TRUE, FALSE)
    # df[nrow(df), "is_first"] <- FALSE 
    plot$x$data[[nrow(df)]]$showlegend <- FALSE # remove text legend
    for (i in df$id) {
        is_first <- df$is_first[[i]]
        is_bool <- df$is_bool[[i]]
        plot$x$data[[i]]$name <- df$legend_group[[i]]
        plot$x$data[[i]]$legendgroup <- plot$x$data[[i]]$name
        if (!is_first) plot$x$data[[i]]$showlegend <- FALSE
        if(is_bool) plot$x$data[[i]]$showlegend <- FALSE
    }
    plot

}

.fixCensoredPointsLegendProfilePlotsPlotly = function(plot) {
    df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
    first_false_index <- which(df$legend_entries == "FALSE")[1]
    first_true_index <- which(df$legend_entries == "TRUE")[1]

    # Update plot data for the first occurrence of "FALSE"
    if (!is.na(first_false_index)) {
        plot$x$data[[first_false_index]]$name <- "Detected data"
        plot$x$data[[first_false_index]]$showlegend <- TRUE
    }

    # Update plot data for the first occurrence of "TRUE"
    if (!is.na(first_true_index)) {
        plot$x$data[[first_true_index]]$name <- "Censored missing data"
        plot$x$data[[first_true_index]]$showlegend <- TRUE
    }
    plot
}

.fixLegendPlotlyPlotsVolcano = function(plot) {
    df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
    # Create a mapping
    color_mapping <- c("black" = "No regulation", "red" = "Up-regulated", "blue" = "Down-regulated")
    # Update the legend_entries column
    df$legend_group <- sapply(df$legend_entries, function(entry) {
        for (color in names(color_mapping)) {
            if (grepl(color, entry)) {
                entry <- gsub(color, color_mapping[color], entry)
                break
            }
        }
        entry <- gsub(",.+", "", entry)
        entry <- gsub("\\(|\\)", "", entry)  # Remove any remaining parentheses
        entry
    })
    for (i in df$id) {
        if(length(grep(df$legend_group[[i]], color_mapping)) == 0) { # keep only 3 legends
            plot$x$data[[i]]$showlegend <- FALSE
        }
        plot$x$data[[i]]$name <- df$legend_group[[i]]
        plot$x$data[[i]]$legendgroup <- plot$x$data[[i]]$name
    }
    plot <- plotly::layout(plot,legend=list(title=list(text="")))
    plot
}

.getPlotlyPlotHTML = function(plots, width, height) {
    doc <- htmltools::tagList(lapply(plots,function(x) htmltools::div(x, style = "float:left;width:100%;")))
    # Set a specific width for each plot
    plot_width <- 800
    plot_height <- 600

    # Create a div for each plot with style settings
    divs <- lapply(plots, function(x) {
        htmltools::div(x, style = paste0("width:", plot_width, "px; height:", plot_height, "px; margin: 10px;"))
    })

    # Combine the divs into a tagList
    doc <- htmltools::tagList(divs)
    doc
}

.savePlotlyPlotHTML = function(plots, address, file_name, width, height) {
    print("Saving plots as HTML")
    pb <- txtProgressBar(min = 0, max = 4, style = 3)
    
    setTxtProgressBar(pb, 1)
    file_name = getFileName(address, file_name, width, height)
    file_name = paste0(file_name,".html")
    
    setTxtProgressBar(pb, 2)
    doc <- .getPlotlyPlotHTML(plots, width, height)
    
    setTxtProgressBar(pb, 3)
    htmltools::save_html(html = doc, file = file_name) # works but lib same folder
    
    setTxtProgressBar(pb, 4)
    zip(paste0(gsub("\\.html$", "", file_name),".zip"), c(file_name, "lib"))
    unlink(file_name)
    unlink("lib",recursive = T)
    
    close(pb)
}
