#' 
#' #############################################
#' ## dataProcessPlots
#' #############################################
#' #' @export
#' #' @import ggplot2 
#' #' @importFrom graphics axis image legend mtext par plot.new title plot
#' #' @importFrom grDevices dev.off hcl pdf
#' 
#' dataProcessPlots = function(
#'     data = data, type = type, featureName = "Transition", ylimUp = FALSE,
#'     ylimDown = FALSE, scale = FALSE, interval = "CI", x.axis.size = 10, 
#'     y.axis.size = 10, text.size = 4, text.angle = 0, legend.size = 7, 
#'     dot.size.profile = 2, dot.size.condition = 3, width = 10, height = 10, 
#'     which.Protein = "all",  originalPlot = TRUE, summaryPlot = TRUE,
#'     save_condition_plot_result = FALSE, 
#'     remove_uninformative_feature_outlier = FALSE, address = ""
#' ) {
#'     
#' }
#' 
#' dataProcessPlots <- function(data=data,
#'                              type=type,
#'                              featureName="Transition",
#'                              ylimUp=FALSE,
#'                              ylimDown=FALSE,
#'                              scale=FALSE,
#'                              interval="CI",
#'                              x.axis.size=10,
#'                              y.axis.size=10,
#'                              text.size=4,
#'                              text.angle=0,
#'                              legend.size=7,
#'                              dot.size.profile=2,
#'                              dot.size.condition=3,
#'                              width=10,
#'                              height=10, 
#'                              which.Protein="all", 
#'                              originalPlot=TRUE, 
#'                              summaryPlot=TRUE,
#'                              save_condition_plot_result=FALSE, 
#'                              remove_uninformative_feature_outlier=FALSE,
#'                              address="") {
#'     
#'     datafeature <- data$ProcessedData
#'     datarun <- data$RunlevelData
#'     
#'     datafeature$PROTEIN <- factor(datafeature$PROTEIN)	
#'     datarun$Protein <- factor(datarun$Protein)	
#'     
#'     if (!is.element("SUBJECT_NESTED", colnames(datafeature))) {
#'         stop("Input for dataProcessPlots function should be processed by dataProcess function previously. Please use 'dataProcess' function first.")
#'     }
#'     checkmate::assertChoice(toupper(type), c("PROFILEPLOT", "QCPLOT", "CONDITIONPLOT"),
#'                             .var.name = "type")
#'     
#'     if (as.character(address) == "FALSE") { 
#'         if (which.Protein == 'all') {
#'             stop('** Cannnot generate all plots in a screen. Please set one protein at a time.')
#'         } else if (length(which.Protein) > 1) {
#'             stop('** Cannnot generate multiple plots in a screen. Please set one protein at a time.')
#'         }
#'     }
#'     
#' }
#' 
#' .plotProfile = function(feature_data, run_data, proteins, original_plot, summary_plot,
#'                         remove_uninformative_feature_outlier) {
#'     is_censored = is.element("censored", colnames(feature_data))
#'     if (is_censored) {
#'         feature_data$is_censored = factor(feature_data$is_censored, levels = c("FALSE", "TRUE"))
#'     }
#'     
#'     if (remove_uninformative_feature_outlier) {
#'         if (is.element("feature_quality", colnames(input))) {
#'             feature_data$ABUDANDANCE = ifelse(
#'                 feature_data$feature_quality == "Noninformative" | feature_data$is_outlier,
#'                 NA, feature_data$ABUNDANCE)
#'             msg = "** Filtered out uninformative feature and outliers in the profile plots."
#'         } else {
#'             msg = "** To remove uninformative features or outliers, please use \"featureSubset == \"highQuality\" option in \"dataProcess\" function."
#'         }
#'         getOption("MSstatsMsg")("INFO", msg)
#'     }
#'     
#'     if (proteins != "all") {
#'         if (is.character(proteins)) {
#'             missing_proteins = setdiff(proteins, unique(feature_data$PROTEIN))
#'             if (length(missing_proteins) > 0) {
#'                 stop(paste0("Please check protein name. Data set does not have this protein. - ", toString(setdiff(proteins, unique(input$PROTEIN)))))
#'             }
#'         }
#'         
#'         if (is.numeric(proteins)) {
#'             temp.name <- levels(feature_data$PROTEIN)[proteins]
#'             if (length(levels(feature_data$PROTEIN)) < max(proteins)) {
#'                 stop(paste0("Please check your selection of proteins. There are ", 
#'                             length(levels(feature_data$PROTEIN))," proteins in this dataset."))
#'             }
#'         }
#'         
#'         feature_data = feature_data[PROTEIN %in% proteins, ]
#'         feature_data$PROTEIN = factor(feature_data$PROTEIN)
#'         run_data = run_data[Protein %in% proteins]
#'         run_data$PROTEIN = factor(run_data$Protein)
#'     }
#'     
#'     ## assign upper or lower limit
#'     # MC, 2016/04/21, default upper limit is maximum log2(intensity) after normalization+3, then round-up
#'     
#'     y.limup <- ceiling(max(datafeature$ABUNDANCE, na.rm=TRUE) + 3)
#'     if (is.numeric(ylimUp)) {
#'         y.limup <- ylimUp 
#'     }
#'     y.limdown <- -1
#'     if (is.numeric(ylimDown)) {
#'         y.limdown <- ylimDown 
#'     }
#'     
#'     datafeature <- datafeature[with(datafeature, order(GROUP_ORIGINAL, SUBJECT_ORIGINAL, LABEL)), ]
#'     datafeature$RUN <- factor(datafeature$RUN, levels=unique(datafeature$RUN), labels=seq(1, length(unique(datafeature$RUN))))
#'     datafeature$RUN <- as.numeric(datafeature$RUN)
#'     tempGroupName <- unique(datafeature[, c("GROUP_ORIGINAL", "RUN")])
#'     groupAxis <- as.numeric(xtabs(~GROUP_ORIGINAL, tempGroupName))
#'     cumGroupAxis <- cumsum(groupAxis)
#'     lineNameAxis <- cumGroupAxis[-nlevels(datafeature$GROUP_ORIGINAL)]
#'     groupName <- data.frame(RUN=c(0, lineNameAxis) + groupAxis / 2 + 0.5, 
#'                             ABUNDANCE=rep(y.limup-1, length(groupAxis)), 
#'                             Name=levels(datafeature$GROUP_ORIGINAL))
#'     if (length(unique(datafeature$LABEL)) == 2) {
#'         datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Reference", "Endogenous"))	
#'     } else {
#'         if (unique(datafeature$LABEL) == "L") {
#'             datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Endogenous"))	
#'         }
#'         if (unique(datafeature$LABEL) == "H") {
#'             datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Reference"))
#'         }
#'     }
#'     
#'     ## need to fill in incomplete rows for Runlevel data
#'     haverun <- FALSE
#'     if (sum(is.element(colnames(datarun), "RUN")) != 0) {
#'         datamat <- dcast( Protein ~ RUN, data=datarun, value.var='LogIntensities', keep=TRUE) 
#'         datarun <- melt(datamat, id.vars=c('Protein'))
#'         colnames(datarun)[colnames(datarun) %in% c("variable", "value")] <- c('RUN', 'ABUNDANCE')
#'         haverun <- TRUE
#'     }
#'     if (any(is.element(colnames(datafeature), "SuggestToFilter"))) {
#'         datafeature$SuggestToFilter <- NULL
#'     }
#'     if (any(is.element(colnames(datafeature), "Filter.Repro"))) {
#'         datafeature$Filter.Repro <- NULL
#'     }
#'     if (any(is.element(colnames(datafeature), "feature_quality"))) {
#'         datafeature$feature_quality <- NULL
#'     }
#'     if (any(is.element(colnames(datafeature), "is_outlier"))) {
#'         datafeature$is_outlier <- NULL
#'     }
#'     
#'     temp <- datafeature[!is.na(datafeature[, "ABUNDANCE"]) & !is.na(datafeature[, "INTENSITY"]), ]
#'     temp <- temp[1, ]
#'     
#'     yaxis_name = .getYaxis(input)
#'     proteins = unique(feature_data$PROTEIN)
#'     
#'     if (originalPlot) {
#'         for (i in 1:nlevels(datafeature$PROTEIN)) {	
#'             single_protein = feature_data[PROTEIN == proteins[i]]
#'             single_protein[, FEATURE := factor(as.character(FEATURE))]
#'             single_protein[, SUBJECT := factor(SUBJECT)]
#'             single_protein[, GROUP_ORIGINAL := factor(GROUP_ORIGINAL)]
#'             single_protein[, SUBJECT_ORIGINAL := factor(SUBJECT_ORIGINAL)]
#'             single_protein[, PEPTIDE := factor(PEPTIDE)]
#'             
#'             if (all(is.na(single_protein$ABUNDANCE))) {
#'                 getOption("MSstatsMsg")("INFO", 
#'                                         message(paste0("Can't the Profile plot for ", 
#'                                                        proteins[i], 
#'                                                        "because all measurements are NAs.")))
#'                 next()
#'             }
#'             b = unique(single_protein[, list(PEPTIDE, FEATURE)])
#'             b = b[order(PEPTIDE, FEATURE)]
#'             temp1 <- xtabs(~b[, 1])
#'             ss <- NULL
#'             s <- NULL
#'             for (j in 1:length(temp1)) {
#'                 temp3 <- rep(j, temp1[j])
#'                 s <- c(s, temp3)
#'                 temp2 <- seq(1, temp1[j])
#'                 ss <- c(ss, temp2)	
#'             }
#'             groupNametemp = data.frame(groupName, 
#'                                        "FEATURE"=unique(single_protein$FEATURE)[1], 
#'                                        "PEPTIDE"=unique(single_protein$PEPTIDE)[1])
#'             dot_colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#'             check_length = length(unique(s)) %/% length(dot_colors)
#'             if (check_length > 0) {
#'                 dot_colors = rep(dot_colors, times = check_length + 1)
#'             }
#'             profile_plot = .makeProfilePlot(single_protein, is_censored, featureName)
#'             print(profile_plot)
#'             getOption("MSstatsMsg")(paste("Drew the Profile plot for ", unique(single_protein$PROTEIN), 
#'                                           "(", i, " of ", length(unique(datafeature$PROTEIN)), ")"))
#'         }
#'     }
#'     
#'     ## 2st plot for original plot : summary ##
#'     if (summaryPlot) {
#'         for (i in 1:nlevels(datafeature$PROTEIN)) {	
#'             
#'             single_protein <- datafeature[datafeature$PROTEIN == levels(datafeature$PROTEIN)[i], ]
#'             single_protein$FEATURE <- factor(as.character(single_protein$FEATURE))	
#'             single_protein$SUBJECT <- factor(single_protein$SUBJECT)	
#'             single_protein$GROUP_ORIGINAL <- factor(single_protein$GROUP_ORIGINAL)	
#'             single_protein$SUBJECT_ORIGINAL <- factor(single_protein$SUBJECT_ORIGINAL)
#'             single_protein$PEPTIDE <- factor(as.character(single_protein$PEPTIDE))
#'             
#'             ## if all measurements are NA,
#'             if (nrow(single_protein) == sum(is.na(single_protein$ABUNDANCE))) {
#'                 message(paste("Can't the Profile plot for ", unique(single_protein$PROTEIN), 
#'                               "(", i, " of ", length(unique(datafeature$PROTEIN)), 
#'                               ") because all measurements are NAs."))
#'                 next()
#'             }
#'             
#'             ## seq for peptide and transition
#'             b <- unique(single_protein[, c("PEPTIDE", "FEATURE")])
#'             b <- b[with(b, order(PEPTIDE, FEATURE)), ] ## add because if there are missing value, orders are different.
#'             
#'             temp1 <- xtabs(~b[, 1])
#'             ss <- NULL
#'             s <- NULL
#'             
#'             for(j in 1:length(temp1)) {
#'                 temp3 <- rep(j, temp1[j])
#'                 s <- c(s, temp3)
#'                 temp2 <- seq(1, temp1[j])
#'                 ss <- c(ss, temp2)	
#'             }
#'             
#'             ## for annotation of condition
#'             groupNametemp <- data.frame(groupName, 
#'                                         FEATURE=unique(single_protein$FEATURE)[1], 
#'                                         analysis="Run summary")
#'             
#'             if (haverun) {
#'                 subrun <- datarun[datarun$Protein == levels(datafeature$PROTEIN)[i], ]
#'                 if (nrow(subrun) != 0) {
#'                     quantrun <- single_protein[1, ]
#'                     quantrun[, 2:ncol(quantrun)] <- NA
#'                     quantrun <- quantrun[rep(seq_len(nrow(subrun))), ]
#'                     quantrun$PROTEIN <- subrun$Protein 
#'                     quantrun$PEPTIDE <- "Run summary"
#'                     quantrun$TRANSITION <- "Run summary" 
#'                     quantrun$FEATURE <- "Run summary" 
#'                     quantrun$LABEL <- "Endogenous"
#'                     quantrun$RUN <- subrun$RUN
#'                     quantrun$ABUNDANCE <- subrun$ABUNDANCE
#'                     quantrun$FRACTION <- 1
#'                 } else { # if there is only one Run measured across all runs, no Run information for linear with censored
#'                     quantrun <- datafeature[1, ]
#'                     quantrun[, 2:ncol(quantrun)] <- NA
#'                     
#'                     quantrun$PROTEIN <- levels(datafeature$PROTEIN)[i]
#'                     quantrun$PEPTIDE <- "Run summary"
#'                     quantrun$TRANSITION <- "Run summary" 
#'                     quantrun$FEATURE <- "Run summary" 
#'                     quantrun$LABEL <- "Endogenous"
#'                     quantrun$RUN <- unique(datafeature$RUN)[1]
#'                     quantrun$ABUNDANCE <- NA
#'                     quantrun$FRACTION <- 1
#'                 }
#'                 
#'                 if (any(is.element(colnames(single_protein), "censored"))) {
#'                     quantrun$censored <- FALSE
#'                 }
#'                 
#'                 quantrun$analysis <- "Run summary"
#'                 single_protein$analysis <- "Processed feature-level data"
#'                 
#'                 ## if 'Filter' column after feature selection, remove this column in order to match columns with run quantification
#'                 filter_column <- is.element(colnames(single_protein), "Filter")
#'                 if (any(filter_column)) {
#'                     single_protein<-single_protein[, !filter_column]
#'                 }
#'                 
#'                 final <- rbind(single_protein, quantrun)
#'                 final$analysis <- factor(final$analysis)
#'                 final$FEATURE <- factor(final$FEATURE)
#'                 final$RUN <- as.numeric(final$RUN)
#'                 
#'                 profile_plot = .makeSummaryProfilePlot(input, is_censored)
#'                 print(profile_plot)
#'                 getOption("MSstatsMsg")("INFO", paste("Drew the Profile plot with summarization for ", unique(single_protein$PROTEIN), 
#'                                                       "(", i, " of ", length(unique(datafeature$PROTEIN)), ")"))
#'             }
#'         } # end-loop for each protein
#'     }  
#' } 
#' 
#' .getYaxis = function(input) {
#'     source = head(input[!is.na(INTENSITY) & !is.na(ABUNDANCE)], 1)
#'     original = source$INTENSITY[1]
#'     log2_diff = abs(original - log(original, 2))
#'     lo10_diff = abs(original - log(original, 10))
#'     if (log2_diff < log10_diff) {
#'         "Log2-intensities"
#'     } else {
#'         "Log10-intensities"
#'     }
#' }
#' 
#' .plotQC = function(input, protein) {
#'     temp <- input[!is.na(input[,"ABUNDANCE"]) & !is.na(input[,"INTENSITY"]), ]
#'     temp <- temp[1, ]
#'     
#'     yaxis_name = .getYaxis(input)
#'     
#'     y.limup = ifelse(is.numeric(ylimUp), ceiling(max(input$ABUNDANCE, na.rm=TRUE) + 3), ylimUp)    
#'     y.limdown = ifelse(is.numeric(ylimDown), -1, ylimDown)
#'     
#'     input = input[order(GROUP_ORIGINAL, SUBJECT_ORIGINAL), ]
#'     input$RUN = factor(input$RUN, levels = unique(input$RUN),
#'                        labels = seq(1, data.table::uniqueN(input$RUN)))
#'     if (nlevels(input$LABEL) == 2) {
#'         input$LABEL = factor(input$LABEL, labels = c("Reference", "Endogenous"))	
#'         label_color = c("darkseagreen1", "lightblue")
#'     } else {
#'         if (unique(input$LABEL) == "L") {
#'             input$LABEL = factor(input$LABEL, labels=c("Endogenous"))
#'             label_color = c("lightblue")	
#'         } else {
#'             input$LABEL = factor(input$LABEL, labels=c("Reference"))
#'             label_color = c("darkseagreen1")
#'         }
#'     }
#' 
#'     tempGroupName <- unique(input[, c("GROUP_ORIGINAL", "RUN")])
#'     input <- input[with(input, order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL)), ]
#'     groupAxis <- as.numeric(xtabs(~GROUP_ORIGINAL, tempGroupName))
#'     cumGroupAxis <- cumsum(groupAxis)
#'     lineNameAxis <- cumGroupAxis[-nlevels(input$GROUP_ORIGINAL)]
#'     groupName <- data.frame("RUN"=c(0, lineNameAxis)+groupAxis / 2 + 0.5, 
#'                             "ABUNDANCE"=rep(y.limup-1, length(groupAxis)), 
#'                             "Name"=levels(input$GROUP_ORIGINAL))
#' 
#'     if (protein %in% c("all", "allonly")) {
#'         qc_plot = .makeQCPlot(input, TRUE)
#'         getOption("MSstatsMsg")("INFO", "Drew the Quality Contol plot(boxplot) for all proteins.")
#'     } else {
#'         protein = .getQCProtein(input, protein)
#'         single_protein = input[PROTEIN == protein]
#'         qc_plot = .makeQCPlot(single_protein, FALSE)
#'         print(qc_plot)
#'         getOption("MSstatsMsg")("INFO", paste("Drew the Quality Control plot (boxplot) for",
#'                                               protein))
#'     }
#' }
#' 
#' 
#' .plotCondition = function(input) {
#'     colnames(datarun)[colnames(datarun) == "Protein"] <- "PROTEIN"
#'     colnames(datarun)[colnames(datarun) == "LogIntensities"] <- "ABUNDANCE"
#'     if (which.Protein != "all") {
#'         if (is.character(which.Protein)) {
#'             temp.name <- which.Protein
#'             if (length(setdiff(temp.name, unique(datarun$PROTEIN))) > 0) {
#'                 stop(paste("Please check protein name. Dataset does not have this protein. -", toString(temp.name), sep=" "))
#'             }
#'         }
#'         if (is.numeric(which.Protein)) {
#'             temp.name <- levels(datarun$PROTEIN)[which.Protein]
#'             if (length(levels(datarun$PROTEIN))<max(which.Protein)) {
#'                 stop(paste("Please check your selection of proteins. There are ", 
#'                            length(levels(datarun$PROTEIN))," proteins in this dataset."))
#'             }
#'         }
#'         datarun <- datarun[which(datarun$PROTEIN %in% temp.name), ]
#'         datarun$PROTEIN <- factor(datarun$PROTEIN)
#'     }
#'     
#'     resultall <- NULL
#'     
#'     temp <- datafeature[!is.na(datafeature[, "ABUNDANCE"]) & !is.na(datafeature[, "INTENSITY"]), ]
#'     temp <- temp[1,]
#'     temptest <- abs(log2(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"]) < abs(log10(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"])
#'     yaxis_name = .getYaxis(input)
#'     
#'     for (i in 1:nlevels(datarun$PROTEIN)) {	
#'         suball <- NULL
#'         sub <- datarun[datarun$PROTEIN == levels(datarun$PROTEIN)[i], ]
#'         sub <- na.omit(sub)	
#'         sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
#'         sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)	
#'         if (nrow(sub) == sum(is.na(sub$ABUNDANCE))) {
#'             message(paste("Can't the Condition plot for ", unique(sub$PROTEIN), 
#'                           "(", i, " of ",length(unique(datarun$PROTEIN)), ") because all measurements are NAs."))
#'             next()
#'         }
#'         sub.mean <- aggregate(ABUNDANCE ~ GROUP_ORIGINAL, data=sub, mean, na.rm=TRUE)
#'         sub.sd <- aggregate(ABUNDANCE ~ GROUP_ORIGINAL, data=sub, sd)
#'         sub.len <- aggregate(ABUNDANCE ~ GROUP_ORIGINAL, data=sub, length)
#'         colnames(sub.mean)[colnames(sub.mean) == "ABUNDANCE"] <- "Mean"
#'         colnames(sub.sd)[colnames(sub.sd) == "ABUNDANCE"] <- "SD"
#'         colnames(sub.len)[colnames(sub.len) == "ABUNDANCE"] <- "numMeasurement"
#'         suball <- merge(sub.mean, sub.sd, by="GROUP_ORIGINAL")
#'         suball <- merge(suball, sub.len, by="GROUP_ORIGINAL")
#'         if (interval == "CI") {
#'             suball$ciw <- qt(0.975, suball$numMeasurement) * suball$SD / sqrt(suball$numMeasurement)
#'         }
#'         if (interval == "SD") {
#'             suball$ciw <- suball$SD
#'         }
#'         if (sum(is.na(suball$ciw)) >= 1) {
#'             suball$ciw[is.na(suball$ciw)] <- 0
#'         }
#'         
#'         y.limup <- ceiling(max(suball$Mean + suball$ciw))
#'         if (is.numeric(ylimUp)) {
#'             y.limup <- ylimUp 
#'         }
#'         y.limdown <- floor(min(suball$Mean - suball$ciw))
#'         if (is.numeric(ylimDown)) {
#'             y.limdown <- ylimDown 
#'         }
#'         
#'         ## re-order (1, 10, 2, 3, -> 1, 2, 3, ... , 10)
#'         suball <- suball[order(suball$GROUP_ORIGINAL), ]
#'         suball <- data.frame(Protein=unique(sub$PROTEIN), suball)
#'         resultall <- rbind(resultall, suball)
#'         if (!scale) {  ## scale: false
#'             tempsummary <- suball
#'             colnames(tempsummary)[colnames(tempsummary) == "GROUP_ORIGINAL"] <- "Label"
#'             ptemp <- ggplot(aes_string(x='Label', y='Mean'), data=tempsummary)+
#'                 geom_errorbar(aes(ymax = Mean + ciw, ymin= Mean - ciw), 
#'                               data=tempsummary, width=0.1, colour="red")+
#'                 geom_point(size = dot.size.condition, colour = "darkred")+
#'                 scale_x_discrete('Condition')+
#'                 scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
#'                 geom_hline(yintercept = 0, linetype = "twodash", colour = "darkgrey", size = 0.6)+
#'                 labs(title=unique(sub$PROTEIN))+
#'                 theme(
#'                     panel.background=element_rect(fill='white', colour="black"),
#'                     panel.grid.major.y = element_line(colour="grey95"),
#'                     panel.grid.minor.y = element_blank(),
#'                     axis.text.x=element_text(size=x.axis.size, colour="black", angle=text.angle),
#'                     axis.text.y=element_text(size=y.axis.size, colour="black"),
#'                     axis.ticks=element_line(colour="black"),
#'                     axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
#'                     axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
#'                     title=element_text(size=x.axis.size+8, vjust=1.5))
#'         } else {
#'             tempsummary <- suball
#'             colnames(tempsummary)[colnames(tempsummary) == "GROUP_ORIGINAL"] <- "Label"
#'             tempsummary$Label <- as.numeric(gsub("\\D", "", unique(tempsummary$Label)))
#'             ptemp <- ggplot(aes_string(x='Label', y='Mean'), data=tempsummary)+
#'                 geom_errorbar(aes(ymax = Mean + ciw, ymin = Mean - ciw), 
#'                               data=tempsummary, width=0.1, colour="red")+
#'                 geom_point(size=dot.size.condition, colour="darkred")+
#'                 scale_x_continuous('Condition', breaks=tempsummary$Label, labels=tempsummary$Label)+
#'                 scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
#'                 geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+
#'                 labs(title=unique(sub$PROTEIN))+
#'                 theme(
#'                     panel.background=element_rect(fill='white', colour="black"),
#'                     panel.grid.major.y = element_line(colour="grey95"),
#'                     panel.grid.minor.y = element_blank(),
#'                     axis.text.x=element_text(size=x.axis.size, colour="black", angle=text.angle),
#'                     axis.text.y=element_text(size=y.axis.size, colour="black"),
#'                     axis.ticks=element_line(colour="black"),
#'                     axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
#'                     axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
#'                     title=element_text(size=x.axis.size+8, vjust=1.5))
#'         }
#'         print(ptemp)
#'         message(paste("Drew the condition plot for ", unique(sub$PROTEIN), 
#'                       "(", i, " of ", length(unique(datarun$PROTEIN)), ")"))
#'     } # end-loop
#'     
#'     if (save_condition_plot_result) {
#'         colnames(resultall)[colnames(resultall) == "GROUP_ORIGINAL"] <- 'Condition'
#'         if (interval == "CI") {
#'             colnames(resultall)[colnames(resultall) == "ciw"] <- '95% CI'
#'         }
#'         if (interval == "SD") {
#'             colnames(resultall)[colnames(resultall) == "ciw"] <- 'SD'
#'         }
#'     }
#' }
#' }
