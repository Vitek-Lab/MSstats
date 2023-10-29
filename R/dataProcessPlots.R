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
#' @param address the name of folder that will store the results. 
#' @param isPlotly This parameter is for MSstatsShiny application for plotly 
#' render, this cannot be used for saving PDF files as plotly do not have 
#' suppprt for PDFs currently. address and isPlotly cannot be set as TRUE at the
#' same time.
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
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly style add_trace plot_ly subplot
#' 
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
  dot.size.condition = 3, width = 10, height = 10, which.Protein = "all",
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
  
  if(isPlotly & address != FALSE) {
      stop("Both isPlotly and address cannot be set at the same time as plotly 
           plots cannot be saved to a PDF, Please set isPlotly to FALSE
           to generate ggplot Plots to a PDF")
  }
  
  if (type == "PROFILEPLOT") {
      plot <- .plotProfile(processed, summarized, featureName, ylimUp, ylimDown,
                           x.axis.size, y.axis.size, text.size, text.angle, legend.size, 
                           dot.size.profile, width, height, which.Protein, originalPlot, 
                           summaryPlot, remove_uninformative_feature_outlier, address, isPlotly)
      # return(plot)
  }
    
  if (type == "QCPLOT") {
      plot <- .plotQC(processed, featureName, ylimUp, ylimDown, x.axis.size, y.axis.size, 
                      text.size, text.angle, legend.size, dot.size.profile, width, height,
                      which.Protein, address, isPlotly)
      
  }
    
  if (type == "CONDITIONPLOT")
    plot <- .plotCondition(processed, summarized, ylimUp, ylimDown, scale, interval,
                   x.axis.size, y.axis.size, text.size, text.angle, legend.size, 
                   dot.size.profile, dot.size.condition, width, height,
                   which.Protein, save_condition_plot_result, address, isPlotly)
  
  if(isPlotly) {
      plotly_plot <- .convert.ggplot.plotly(plot)
      if(toupper(featureName) == "NA") {
          plotly_plot <- plotly_plot %>% style(plotly_plot, showlegend = FALSE)
      }
      plotly_plot
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
  
  all_proteins = levels(processed$PROTEIN)
  if (originalPlot) {
    savePlot(address, "ProfilePlot", width, height)
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
      
      # profile_plot_plotly = .makeProfilePlotPlotly(single_protein, profile_plot, featureName)
      setTxtProgressBar(pb, i)
      print(profile_plot)
      if(isPlotly & address == FALSE) {
          return(profile_plot)
      }
    }
    close(pb)
    
    if (address != FALSE) {
      dev.off()
    } 
  }
  
  if (summaryPlot) {
    protein_by_run = expand.grid(Protein = unique(summarized$Protein), 
                                 RUN = unique(summarized$RUN))
    summarized = merge(summarized, protein_by_run, by = c("Protein", "RUN"),
                       all.x = TRUE, all.y = TRUE)
    savePlot(address, "ProfilePlot_wSummarization", width, height)
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
      single_protein$analysis = "Processed feature-level data"
      combined = rbind(single_protein[
        , 
        list(PROTEIN, PEPTIDE, TRANSITION, FEATURE, LABEL,
             RUN, ABUNDANCE, FRACTION, censored, analysis)], quant)
      combined$analysis = factor(combined$analysis)
      combined$FEATURE = factor(combined$FEATURE)
      combined$RUN = as.numeric(combined$RUN)
      profile_plot = .makeSummaryProfilePlot(
        combined, is_censored, y.limdown, y.limup, x.axis.size, y.axis.size, 
        text.size, text.angle, legend.size, dot.size.profile, cumGroupAxis, 
        yaxis.name, lineNameAxis, groupNametemp
      )
      # profile_plot_plotly = .makeProfilePlotPlotly(combined, profile_plot, featureName)
      print(profile_plot)
      setTxtProgressBar(pb, i)
      
      if(isPlotly & address == FALSE) {
          return(profile_plot)
      }
    }
    close(pb)
    
    if (address != FALSE) {
      dev.off()
    } 
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
  
  savePlot(address, "QCPlot", width, height)
  if (protein %in% c("all", "allonly")) {
    qc_plot = .makeQCPlot(processed, TRUE, y.limdown, y.limup, x.axis.size, 
                          y.axis.size, text.size, text.angle, legend.size, 
                          label.color, cumGroupAxis, groupName, lineNameAxis, 
                          yaxis.name)
    print(qc_plot)
    if(isPlotly & address == FALSE) {
        return(qc_plot)
    }
  } 
  
  if (protein != 'allonly') {
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
      setTxtProgressBar(pb, i)
      
      if(isPlotly & address == FALSE) {
          return(qc_plot)
      }
    } 
    close(pb)
  } 
  if (address != FALSE) {
    dev.off()
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
  savePlot(address, "ConditionPlot", width, height)
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
    setTxtProgressBar(pb, i)
    
    if(isPlotly & address == FALSE) {
        return(con_plot)
    }
  }
  close(pb)
  
  if (address != FALSE) {
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
}

#' converter for plots from ggplot to plotly
#' @noRd
.convert.ggplot.plotly = function(plot) {
    converted_plot <- ggplotly(plot)
    converted_plot <- converted_plot %>% 
        plotly::layout(
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
                    size = 9  # Set the font size for legend item labels
                ),
                title = list(
                    font = list(
                        size = 12  # Set the font size for the legend title
                    )
                )
            )
        )
    fix_plot = .fix.legend.plotly.plots(converted_plot)
    return(fix_plot)
}

.fix.legend.plotly.plots = function(plot) {

    df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
    
    df$legend_group <- gsub("^(.*?),.*", "\\1", df$legend_entries)
    df$is_first <- !duplicated(df$legend_group)
    df$is_bool <- ifelse(grepl("TRUE|FALSE", df$legend_group), TRUE, FALSE)
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
