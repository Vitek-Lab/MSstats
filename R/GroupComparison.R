#' Perform group comparisons using lm/lme4
#'
#' This function handles the comparison of intensities of each protein using a
#' linear model of the log2(intensities).  Conceptually it seems very similar to
#' what one does when using limma.
#'
#' @param contrast.matrix  A contrast matrix for comparing
#' @param data the result from processedData()
#' @param verbose Print the logging data?
#' @param scopeOfBioReplication Define the model with respect to biological replicates.
#' @param scopeOfTechReplication Define the model with respect to technical replicates.
#' @return a comparison
#'
#' @export
groupComparison <- function(contrast.matrix=NULL, data=NULL, verbose=TRUE,
                            scopeOfBioReplication="expanded",
                            scopeOfTechReplication="restricted") {
  allfiles <- list.files()
  logging::loginfo("MSstats - starting groupComparison().")

  rawinput <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                "ProductCharge", "IsotopeLabelType", "Condition",
                "BioReplicate", "Run", "Intensity")

  if (length(setdiff(toupper(rawinput), toupper(colnames(data[["ProcessedData"]])))) == 0) {
    logging::logwarn(paste0("The required input data was not from dataProcess()."))
    stop("Please use 'dataProcess()' first. Then groupComparison().")
  }

  ## contrast. matrix
  if (ncol(contrast.matrix) != length(unique(data[["ProcessedData"]][["GROUP_ORIGINAL"]]))) {
    logging::loginfo("The number of columns and groups are not the same.")
    stop("Check the contrast matrix. The groups are different than the columns of the contrast matrix.")
  }

  ## check whether row.names of contrast.matrix.sub exists or not
  if (sum(is.null(row.names(contrast.matrix))) > 0) {
    logging::loginfo("The contrast matrix requires row names.")
    stop("No row.names of comparison exist.\n")
  }

  if (!(scopeOfBioReplication == "restricted" | scopeOfBioReplication == "expanded")) {
    logging::loginfo(paste0("The variable scopeOfBioReplication is: ", scopeOfBioReplicateion, " is invalid."))
    stop("'scopeOfBioReplication' must be one of 'restricted' or 'expanded'.")
  }

  labeled <- ifelse(length(unique(data[["ProcessedData"]][["LABEL"]])) == 1, FALSE, TRUE)
  logging::loginfo(paste0("Some important variables: labeled=", labeled,
                          " scopeOfBioReplication=", scopeOfBioReplication))
  repeated <- checkRepeated(data[["ProcessedData"]])

  if (repeated) {
    logging::loginfo("Time course design of experiment is acceptable.")
  } else {
    logging::loginfo("Case control design of experiment is acceptable.")
  }

  ## for final result report
  out <- NULL
  outsummary <- NULL
  outfitted <- NULL
  dataafterfit <- NULL
  ## check input for labeled
  ## start to analyze by protein ID
  logging::loginfo("Starting inference methods.")
  if (data[["SummaryMethod"]] == "logOfSum") {
    logging::loginfo("Use t-test for log sum summarization per subject.")
    if (repeated) {
      logging::logwarn("Only group comparisons are available for t-testing.")
      stop()
    }
  }

  ## need original group information
  rqall <- data[["RunlevelData"]]
  rqall[["Protein"]] <- factor(rqall[["Protein"]])
  processall <- data[["ProcessedData"]]
  origGroup <- unique(rqall[["GROUP_ORIGINAL"]])
  groupinfo <- levels(data[["ProcessedData"]][["GROUP_ORIGINAL"]])

  ## It seems to me this for loop can be redone using dplyr, but since I don't
  ## truly understand wtf it is doing yet, such efforts are doomed to miserable
  ## failure.
  ##rq_tbl <- as_tibble(rqall)
  ##logfcs <- rq_tbl %>%
  ##  group_by(Protein) %>%
  ##  dplyr::do(gather_protein_fc(., contrast.matrix, origGroup, groupinfo)

  pb <- txtProgressBar(max=nlevels(rqall[["Protein"]]), style = 3)
  for (i in 1:nlevels(rqall[["Protein"]])) {
    subset_idx <- rqall[["Protein"]] == levels(rqall[["Protein"]])[i]
    subset <- rqall[subset_idx, ]
    setTxtProgressBar(pb, i)
    colnames(subset)[colnames(subset) == "LogIntensities"] <- "ABUNDANCE"
    colnames(subset)[colnames(subset) == "Protein"] <- "PROTEIN"
    ## it is important to remove NA first, before we have the correct structure for each factor
    sub <- subset[!is.na(subset[["ABUNDANCE"]]), ]
    ## 1. for logsum, t-test

    ## Provide a null fitted result in case something goes wrong.
    fitted <- list(result=NULL, valueresid=NULL, valuefitted=NULL, fittedmodel=NULL)
    method <- data[["SummaryMethod"]]
    if (method == "logOfSum") {
      ## subject level
      sub[["GROUP"]] <- factor(sub[["GROUP"]])
      sub[["SUBJECT"]] <- factor(sub[["SUBJECT"]])
      sub[["GROUP_ORIGINAL"]] <- factor(sub[["GROUP_ORIGINAL"]])
      ## testing and inference in whole plot
      if (message.show) {
        message(paste("Testing a comparison for protein ",
                      unique(sub[["PROTEIN"]]), "(", i, " of ",
                      length(unique(rqall[["Protein"]])), ")"))
      }
      fitted <- try(ttest.logsum(contrast.matrix, sub, origGroup), silent=TRUE)
    } else {
      ## linear model
      sub[["GROUP"]] <- factor(sub[["GROUP"]])
      sub[["SUBJECT"]] <- factor(sub[["SUBJECT"]])
      sub[["GROUP_ORIGINAL"]] <- factor(sub[["GROUP_ORIGINAL"]], levels=groupinfo)
      sub[["SUBJECT_ORIGINAL"]] <- factor(sub[["SUBJECT_ORIGINAL"]])
      sub[["SUBJECT_NESTED"]] <- factor(sub[["SUBJECT_NESTED"]])
      sub[["RUN"]] <- factor(sub[["RUN"]])

      ## singleFeature <- checkSingleFeature(sub)
      singleSubject <- checkSingleSubject(sub)
      TechReplicate <- checkTechReplicate(sub) ## use for label-free model
      ## MissGroupByFeature <- checkMissGroupByFeature(sub)
      ## MissRunByFeature <- checkMissRunByFeature(sub)
      MissSubjectByGroup <- checkRunbyFeature(sub)
      UnequalSubject <- checkUnequalSubject(sub)
      ## testing and inference in whole plot
      if (isTRUE(verbose)) {
        logging::loginfo(paste0("Testing a comparison of protein: ",
                                unique(sub[["PROTEIN"]]), "(", i, " of ",
                                length(unique(rqall[["Protein"]])), ")"))
      }
      ## fit the model
      fitted <- try(fit.model.single(contrast.matrix, sub, origGroup,
                                     TechReplicate=TechReplicate,
                                     singleSubject=singleSubject,
                                     repeated=repeated))
    }

    ##
    null_result <- NULL
    if (class(fitted) == "try-error") {
      logging::logwarn(paste0("Cannot analyze ", levels(rqall[["Protein"]])[i], " for comparison."))
      tmpresult <- data.frame()
      for (k in 1:nrow(contrast.matrix)) {
        fitted[["result"]] <- rbind(tempresult[["result"]],
                                   data.frame(
                                     "Protein" = levels(rqall[["Protein"]])[i],
                                     "Label" = row.names(contrast.matrix)[k],
                                     "logFC" = NA,
                                     "SE" = NA,
                                     "Tvalue" = NA,
                                     "DF" = NA,
                                     "pvalue" = NA,
                                     "issue" = NA))
      }
    }
    fitted_data <- fitted[["result"]]

    ## need to add information about % missingness and %imputation
    if (is.element("NumImputedFeature", colnames(data$RunlevelData))) {
      subprocess <- processall[processall[["PROTEIN"]] ==
                               as.character(unique(sub[["PROTEIN"]])), ]
      missing_pct <- count.missing.percentage(contrast.matrix, fitted_data, sub, subprocess)
      fitted_data <- missing_pct
    }
    ## comparison result table
    fitted_data[["Label"]] <- as.character(fitted_data[["Label"]])
    out <- data.table::rbindlist(list(out, fitted_data), fill=TRUE)

    ## for checking model assumptions
    ## add residual and fitted after fitting the model
    if (method != "logOfSum" & class(fitted) == "try-error") {
      if (nrow(sub) != 0) {
        sub[["residuals"]] <- NA
        sub[["fitted"]] <- NA
      }
    } else if (method != "logOfSum" & any(is.na(fitted[["result"]][["pvalue"]]))) {
      ## even though there is no error for fit.model function, still output can have empty fitted and residuals.
      sub[["residuals"]] <- NA
      sub[["fitted"]] <- NA
    } else if (method != "logOfSum") {
      sub[["residuals"]] <- fitted_data[["valueresid"]]
      sub[["fitted"]] <- fitted_data[["valuefitted"]]
    }

    dataafterfit <- data.table::rbindlist(list(dataafterfit, sub), fill=TRUE)
    ## save fitted model
    outfitted <- c(outfitted, list(fitted[["fittedmodel"]]))
  } ## End of this absurdist for loop.
  close(pb)

  logging::loginfo("Comparisons for all proteins are finished.")
  ## finalize result
  ## need to FDR per comparison
  out.all <- NULL
  out[["Label"]] <- as.character(out[["Label"]])

  ## What a bizarre way to add adjusted p-values.
  for (i in 1:length(unique(out[["Label"]]))) {
    outsub <- out[out[["Label"]] == unique(out[["Label"]])[i], ]
    outsub <- data.frame(outsub, adj.pvalue=p.adjust(outsub[["pvalue"]], method="BH"))
    out.all <- data.table::rbindlist(list(out.all, outsub))
  }
  out.all[["Label"]] <- factor(out.all[["Label"]])
  logging::loginfo("Adjusted p-values are calculated.")

  ## Figure out if the column for logFC should be log2 or log10
  ## There are much easier ways of doing this.
  ## Perhaps just pass the state of logTrans from dataProcess() into this
  ## function or save that value as an element in the dataProcess list?
  ## Whatever, there are more important strange things to understand than this.
  temp <- data[["ProcessedData"]][!is.na(data[["ProcessedData"]][, "ABUNDANCE"]) &
                             !is.na(data[["ProcessedData"]][, "INTENSITY"]), ]
  temp <- temp[temp[["ABUNDANCE"]] > 2, ]

  if (abs(log2(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"]) <
      abs(log10(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"])) {
    colnames(out.all)[3] <- "log2FC"
  }
  if (abs(log2(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"]) >
      abs(log10(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"])) {
    colnames(out.all)[3] <- "log10FC"
  }

  ## change the format as data.frame
  out.all <- data.frame(out.all)
  ## change order of columns,
  if (any(is.element(colnames(out.all), "ImputationPercentage"))) {
    out.all <- out.all[, c(1:7, 11, 8, 9, 10)]
  } else if (any(is.element(colnames(out.all), "MissingPercentage"))) {
    out.all <- out.all[, c(1:7, 10, 8, 9)]
  } else {
    ## logOfSum output
    out.all <- out.all[, c(1:7, 9)]
  }

  ## if just one condition is completely missing, replace adjust pvalue as zero
  if (any(is.element(colnames(out.all), "issue"))){
    out.all[!is.na(out.all[["issue"]]) & out.all[["issue"]] == "oneConditionMissing", "adj.pvalue"] <- 0
  }

  logging::loginfo("Group comparison is complete.")
  finalout <- list(
    "ComparisonResult" = out.all,
    "ModelQC" = dataafterfit,
    "fittedmodel" = outfitted)
  return(finalout)
}


#' Plot some QC metrics given a statistical model
#'
#' The first thing I am going to do is rip out the pdf calls, if I want to put
#' my plots into a file, I can do it myself.
#'
#' @param data processedData()
#' @param type the type of plot
#' @param axis.size  How big for the axes?
#' @param dot.size  How big for the dots?
#' @param text.size  How big for the text?
#' @param legend.size  And the legend?
#' @param width  How wide?
#' @param height and high?
#' @param featureName  Include the peptides on the plot?
#' @param feature.QQPlot  How many qqs to make?
#' @param which.Protein  Include which proteins?
#' @param address  That's it, I'm getting my tinfoil hat!
#' @return plots!
#' @export
modelBasedQCPlots <- function(data, type, axis.size=10, dot.size=3, text.size=7, legend.size=7,
                              width=10, height=10, featureName=TRUE, feature.QQPlot="all",
                              which.Protein="all", address="") {
  if (length(setdiff(toupper(type),
                     c("QQPLOTS", "RESIDUALPLOTS"))) != 0) {
    stop("Input for type=", type,
         ". However,'type' should be one of 'QQplots', 'ResidualPlots' .")
  }
  data <- data[!is.na(data[["fitted"]]), ]
  data <- data[!is.na(data[["residuals"]]), ]
  data$PROTEIN <- factor(data[["PROTEIN"]])

  ## choose Proteins or not
  if (which.Protein != "all") {
    ## check which.Protein is name of Protein
    if (is.character(which.Protein)) {
      temp.name <- which.Protein
      ## message if name of Protein is wrong.
      if (length(setdiff(temp.name, unique(data[["PROTEIN"]]))) > 0) {
        stop("Please check protein name. Dataset does not have this protein. -",
             toString(temp.name))
      }
    }

    ## check which.Protein is order number of Protein
    if (is.numeric(which.Protein)) {
      temp.name <- levels(data[["PROTEIN"]])[which.Protein]
      ## message if name of Protein is wrong.
      if (length(levels(data[["PROTEIN"]])) < max(which.Protein))
        stop("Please check your selection of proteins. There are ",
             length(levels(data[["PROTEIN"]])), " proteins in this dataset.")
    }
    ## use only assigned proteins
    data <- data[which(data$PROTEIN %in% temp.name), ]
    data$PROTEIN <- factor(data$PROTEIN)
  }

#############################################
### normality, QQ plot
#############################################
  if (toupper(type)=="QQPLOTS") {
    ## test feature.QQPlot is something wrong.
    ## one QQplot per protein, all features together
    if (toupper(feature.QQPlot)=="ALL") {
      ## save the plots as pdf or not
      ## If there are the file with the same name, add next numbering at the end of file name
      if (address!=FALSE) {
        allfiles <- list.files()
        num <- 0
        filenaming <- paste0(address, "QQPlot_allFeatures")
        finalfile <- paste0(address, "QQPlot_allFeatures.pdf")
        while (is.element(finalfile, allfiles)) {
          num <- num + 1
          finalfile <- paste0(paste(filenaming, num, sep="-"), ".pdf")
        }
        pdf(finalfile, width=width, height=height)
      }

      for (i in 1:nlevels(data$PROTEIN)) {
        sub <- data[data$PROTEIN == levels(data$PROTEIN)[i], ]
        ## get slope and intercept for qline
        y <- quantile(sub$residuals[!is.na(sub$residuals)], c(0.25, 0.75))
        x <- qnorm(c(0.25, 0.75))
        slope <- diff(y) / diff(x)
        int <- y[1L] - slope * x[1L]

        ptemp <- ggplot(sub, aes_string(sample="residuals")) +
          geom_point(stat="qq", alpha=0.8, shape=1, size=dot.size) +
          scale_shape(solid=FALSE) +
          geom_abline(slope = slope, intercept = int, colour="red") +
          scale_y_continuous("Sample Quantiles") +
          scale_x_continuous("Theoretical Quantiles") +
          labs(title=paste("Normal Q-Q Plot (", unique(sub$PROTEIN), ")")) +
          theme(
            panel.background=element_rect(fill="white", colour="black"),
            panel.grid.major=element_line(colour="grey95"),
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(size=axis.size, colour="black"),
            axis.text.y=element_text(size=axis.size, colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=axis.size + 5, vjust=-0.4),
            axis.title.y=element_text(size=axis.size + 5, vjust=0.3),
            title=element_text(size=axis.size + 8, vjust=1.5),
            legend.position="none")
        print(ptemp)
        message("Drew the QQ plot for ",
                unique(sub$PROTEIN), "(",
                i, " of ", length(unique(data$PROTEIN)), ")")
      } ## end loop

      if (address!=FALSE) {
        dev.off()
      }
    } ## end QQplot by all feature


    ## panel for each feature,
    if (toupper(feature.QQPlot) == "BYFEATURE") {
      ## save the plots as pdf or not
      ## If there are the file with the same name, add next numbering at the end of file name
      if (address != FALSE) {
        allfiles <- list.files()
        num <- 0
        filenaming <- paste0(address, "QQPlot_byFeatures")
        finalfile <- paste0(address, "QQPlot_byFeatures.pdf")
        while (is.element(finalfile, allfiles)) {
          num <- num + 1
          finalfile <- paste0(paste(filenaming, num, sep="-"), ".pdf")
        }
        pdf(finalfile, width=width, height=height)
      }

      for (i in 1:nlevels(data$PROTEIN)) {
        sub <- data[data$PROTEIN == levels(data$PROTEIN)[i], ]
        ## label-free
        if (length(unique(sub$LABEL)) == 1) {
          ## need to update for qline per feature
          ptemp <- ggplot(sub, aes_string(sample="residuals", color="FEATURE")) +
            geom_point(stat="qq", alpha=0.8, size=dot.size) +
            facet_wrap(~ FEATURE) +
            scale_y_continuous("Sample Quantiles") +
            scale_x_continuous("Theoretical Quantiles") +
            labs(title=paste("Normal Q-Q Plot (", unique(sub$PROTEIN), ")")) +
            theme(
              panel.background=element_rect(fill="white", colour="black"),
              panel.grid.major=element_line(colour="grey95"),
              panel.grid.minor=element_blank(),
              axis.text.x=element_text(size=axis.size, colour="black"),
              axis.text.y=element_text(size=axis.size, colour="black"),
              axis.ticks=element_line(colour="black"),
              axis.title.x=element_text(size=axis.size + 5, vjust=-0.4),
              axis.title.y=element_text(size=axis.size + 5, vjust=0.3),
              title=element_text(size=axis.size + 8, vjust=1.5),
              strip.text.x=element_text(size=text.size),
              legend.position="none")
          print(ptemp)
          message("Drew the QQ plot for ",
                  unique(sub$PROTEIN), "(", i, " of ",
                  length(unique(data$PROTEIN)), ")")
        }

        ## label-based : seperate endogenous and reference
        if (length(unique(sub$LABEL)) == 2) {
          ## need to update for qline per feature
          ## endogenous intensities
          sub.l <- sub[sub$LABEL=="L", ]
          ptemp <- ggplot(sub.l, aes_string(sample="residuals", color="FEATURE")) +
            geom_point(stat="qq", alpha=0.8, size=dot.size) +
            facet_wrap(~ FEATURE) +
            scale_y_continuous("Sample Quantiles") +
            scale_x_continuous("Theoretical Quantiles") +
            labs(title=paste("Normal Q-Q Plot (",
                             unique(sub$PROTEIN), ") - Endogenous Intensities")) +
            theme(
              panel.background=element_rect(fill="white", colour="black"),
              panel.grid.major=element_line(colour="grey95"),
              panel.grid.minor=element_blank(),
              axis.text.x=element_text(size=axis.size, colour="black"),
              axis.text.y=element_text(size=axis.size, colour="black"),
              axis.ticks=element_line(colour="black"),
              axis.title.x=element_text(size=axis.size + 5, vjust=-0.4),
              axis.title.y=element_text(size=axis.size + 5, vjust=0.3),
              title=element_text(size=axis.size + 8, vjust=1.5),
              strip.text.x=element_text(size=text.size),
              legend.position="none")
          print(ptemp)

          ## reference intensities
          sub.h <- sub[sub$LABEL=="H", ]
          ptemp <- ggplot(sub.h, aes_string(sample="residuals", color="FEATURE")) +
            geom_point(stat="qq", alpha=0.8, size=dot.size) +
            facet_wrap(~ FEATURE) +
            scale_y_continuous("Sample Quantiles") +
            scale_x_continuous("Theoretical Quantiles") +
            labs(title=paste("Normal Q-Q Plot (", unique(sub$PROTEIN), ") - Reference Intensities")) +
            theme(
              panel.background=element_rect(fill="white", colour="black"),
              panel.grid.major=element_line(colour="grey95"),
              panel.grid.minor=element_blank(),
              axis.text.x=element_text(size=axis.size, colour="black"),
              axis.text.y=element_text(size=axis.size, colour="black"),
              axis.ticks=element_line(colour="black"),
              axis.title.x=element_text(size=axis.size + 5, vjust=-0.4),
              axis.title.y=element_text(size=axis.size + 5, vjust=0.3),
              title=element_text(size=axis.size + 8, vjust=1.5),
              strip.text.x=element_text(size=text.size),
              legend.position="none")
          print(ptemp)
          message("Drew the QQ plot for ",
                  unique(sub$PROTEIN), "(", i, " of ",
                  length(unique(data$PROTEIN)), ")")
        }
      } ## end loop

      if (address!=FALSE) {
        dev.off()
      }
    } ### end for QQ by FEATURE
  } ### end QQplots

#############################################
#### Residual plot
#############################################
  if (toupper(type) == "RESIDUALPLOTS") {
    y.limdown <- min(data$residuals, na.rm=TRUE)
    y.limup <- max(data$residuals, na.rm=TRUE)
    x.limdown <- min(data$fitted, na.rm=TRUE)
    x.limup <- max(data$fitted, na.rm=TRUE)

    ## save the plots as pdf or not
    ## If there are the file with the same name, add next numbering at the end of file name
    if (address != FALSE) {
      allfiles <- list.files()
      num <- 0
      filenaming <- paste0(address, "ResidualPlot")
      finalfile <- paste0(address, "ResidualPlot.pdf")
      while (is.element(finalfile, allfiles)) {
        num <- num + 1
        finalfile <- paste0(paste(filenaming, num, sep="-"), ".pdf")
      }
      pdf(finalfile, width=width, height=height)
    }

    for (i in 1:nlevels(data$PROTEIN)) {
      sub <- data[data$PROTEIN==levels(data$PROTEIN)[i], ]
      sub$PEPTIDE <- factor(sub$PEPTIDE)
      sub$FEATURE <- factor(sub$FEATURE)
      ptemp <- ggplot(aes_string(x="fitted", y="residuals", color="FEATURE", shape="LABEL"),
                      data=sub) +
        geom_point(size=dot.size, alpha=0.5) +
        geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6) +
        scale_y_continuous("Residuals", limits=c(y.limdown, y.limup)) +
        scale_x_continuous("Predicted Abundance", limits=c(x.limdown, x.limup)) +
        labs(title=levels(data$PROTEIN)[i])
      if (length(unique(sub$LABEL))==2) {
        ptemp <- ptemp +
          scale_shape_manual(values=c(2, 19),
                             name="",
                             labels=c("Reference", "Endogenous"))
      } else{
        ptemp <- ptemp +
          scale_shape_manual(values=c(19),
                             name="",
                             labels=c("Endogenous"))
      }

      if (featureName) {
        ptemp <- ptemp +
          theme(
            panel.background=element_rect(fill="white", colour="black"),
            legend.key=element_rect(fill="white", xocolour="white"),
            legend.text=element_text(size=legend.size),
            panel.grid.major=element_line(colour="grey95"),
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(size=axis.size, colour="black"),
            axis.text.y=element_text(size=axis.size, colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=axis.size + 5, vjust=-0.4),
            axis.title.y=element_text(size=axis.size + 5, vjust=0.3),
            title=element_text(size=axis.size + 8, vjust=1.5),
            legend.position=c("top")) +
          guides(color=guide_legend(ncol=2))
      } else {
        ptemp <- ptemp +
          theme(
            panel.background=element_rect(fill="white", colour="black"),
            panel.grid.major=element_line(colour="grey95"),
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(size=axis.size, colour="black"),
            axis.text.y=element_text(size=axis.size, colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=axis.size + 5, vjust=-0.4),
            axis.title.y=element_text(size=axis.size + 5, vjust=0.3),
            title=element_text(size=axis.size + 8, vjust=1.5),
            legend.position="none")
      }
      print(ptemp)
      message("Drew the residual plot for ",
              unique(sub$PROTEIN), "(", i, " of ",
              length(unique(data[["PROTEIN"]])), ")")
    }
    if (address!=FALSE) {
      dev.off()
    }
  } ## end residualplots
}


#' Check if measurements are missing for entire group.
#'
#' If there is stuff missing from the experimental design, then it might be
#' problematic.
#'
#' @param sub1 I dunno
#' @param contrast.matrix  The matrix to test.
#' @return a list of the test.
#'
#' ok, so why not use qr()?
checkGroupComparisonAgreement <- function(sub1, contrast.matrix) {
  tempSub <- as.numeric(as.character(unique(sub1[, c("GROUP")])))
  positionMiss <- setdiff(seq(1, length(contrast.matrix)), tempSub)
  contrast.matrix.sub1 <- contrast.matrix[tempSub]
  ## either one of the groups for the comparison of interest is not present
  retlist <- list(
    "sign" = length(setdiff(contrast.matrix[tempSub], 0)) < 2,
    "positionMiss" = positionMiss)
  return(retlist)
}

#############################################
## check repeated (case-control? or time-course?)
#############################################
checkRepeated <- function(data) {
  data.light <- data[data[["LABEL"]] == "L", ]
  subjectByGroup <- table(data.light[["SUBJECT_ORIGINAL"]], data.light[["GROUP_ORIGINAL"]])
  subjectAppearances <- apply(subjectByGroup, 1, function(x) sum(x > 0))
  crossedIndicator <- any(subjectAppearances > 1)
  return(crossedIndicator)
}

#############################################
## check single subject for both case-control and time-course?
#############################################
checkSingleSubject <- function(data) {
  temp <- unique(data[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
  temp[["GROUP_ORIGINAL"]] <- factor(temp[["GROUP_ORIGINAL"]])
  temp1 <- xtabs(~ GROUP_ORIGINAL, data=temp)
  singleSubject <- all(temp1 == "1")
  return(singleSubject)
}

#############################################
## check checkSingleFeature
#############################################
checkSingleFeature <- function(data) {
  sigleFeature <- length(unique(data[["FEATURE"]])) < 2
  return(sigleFeature)
}

#############################################
## check checkTechReplicate
#############################################
checkTechReplicate <- function(data) {
  temp <- unique(data[, c("SUBJECT_NESTED", "RUN")])
  temp[["SUBJECT_NESTED"]] <- factor(temp[["SUBJECT_NESTED"]])
  temp1 <- xtabs(~ SUBJECT_NESTED, data=temp)
  TechReplicate <- all(temp1 != "1")
  return(TechReplicate)
}

#############################################
## check checkRunByFeature
#############################################
## it might not be right
checkRunbyFeature <- function(data) {
  data.light <- data[data[["LABEL"]] == "L", ]
  RunByFeature <- table(data.light[["RUN"]], data.light[["FEATURE"]])
  emptyRow <- apply(RunByFeature, 1, sum)
  noRunFeature <- any(emptyRow == 0)
  return(noRunFeature)
}

#############################################
## check checkMissGroupByFeature
#############################################
checkMissGroupByFeature <- function(data) {
  temp <- unique(data[, c("GROUP", "FEATURE")])
  temp1 <- xtabs(~ GROUP, data=temp)
  return(any(temp1 != temp1[1]))
}

#############################################
## check checkMissRunByFeature
#############################################
checkMissRunByFeature <- function(data) {
  ## temp <- unique(data[,c("RUN","FEATURE")])
  temp <- unique(data[data[["LABEL"]] == "L", c("RUN", "FEATURE")])
  temp1 <- xtabs(~ RUN, data=temp)
  ##return(any(temp1!=temp1[1]))
  return(any(temp1 != length(unique(data[["FEATURE"]]))))
}

#############################################
## check checkMissFeature for label-free missingness
#############################################
checkMissFeature <- function(data) {
  dataByPeptide <- tapply(as.numeric(data[["ABUNDANCE"]]),
                          list(data[["FEATURE"]], data[["GROUP_ORIGINAL"]]),
                          function(x) sum(x > 0, na.rm = TRUE))
  missPeptideInd <- apply(dataByPeptide, 1, function(x) any(x == 0 | is.na(x)))
  missingPeptides <- names(missPeptideInd)[missPeptideInd == TRUE]
  return(missingPeptides)
}

#############################################
## check checkUnequalSubject
#############################################
checkUnequalSubject <- function(data) {
  ##temp <- unique(data[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
  temp <- unique(data[data[["LABEL"]] == "L", c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
  temp1 <- xtabs(~ GROUP_ORIGINAL, data=temp)
  return(any(temp1 != temp1[1]))
}


#' Use lm to attempt to get a delta log(intensity) (abundance) metric.
#'
#' @param contrast.matrix  Apparently manually generated contrast matrix.
#' @param sub  Data subset to pass to lm (why not use all of the data for this?)
#' @param origGroup Defined right before the for loop which feeds this function
#'   as unique(rqall[["GROUP_ORIGINAL"]])
#' @param TechReplicate  Boolean from checkTechReplicate()
#' @param singleSubject Boolean from checkSingleSubject()
#' @param repeated  Boolean asking if this is repeated data, I presume, but need
#'   to doublecheck.
#'
#' Come _on_ choose a friggin standard for naming functions/variables.
fit.model.single <- function(contrast.matrix, sub, origGroup, TechReplicate=TRUE,
                             singleSubject=TRUE, repeated=TRUE) {
  ## input : output of run quantification
  sub[["GROUP"]] <- factor(sub[["GROUP"]])
  sub[["SUBJECT"]] <- factor(sub[["SUBJECT"]])
  ## if there is only one condition between two conditions, make error message and next
  if (length(unique(sub[["GROUP_ORIGINAL"]])) == 1) {
    logging::loginfo("Fitting a single group.")
    ## each comparison
    allout <- NULL
    for (k in 1:nrow(contrast.matrix)) {
      ## choose each comparison
      contrast.matrix.sub <- matrix(contrast.matrix[k, ], nrow=1)
      rownames(contrast.matrix.sub) <- rownames(contrast.matrix)[k]
      if (any(levels(origGroup)[contrast.matrix.sub != 0] == unique(sub[["GROUP_ORIGINAL"]]))) {
        logging::logwarn(paste0("Results of protein ",
                                toString(unique(sub[["PROTEIN"]])),
                                " for comparison ", rownames(contrast.matrix.sub)))
        ## need to check Inf vs -Inf
        ##if( contrast.matrix.sub[levels(origGroup)[contrast.matrix.sub != 0] == unique(data2$GROUP_ORIGINAL)] > 0 ){
        if (contrast.matrix.sub[contrast.matrix.sub != 0 &
                                (levels(origGroup) == unique(sub[["GROUP_ORIGINAL"]]))] > 0) {
          out <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "Label" = row.names(contrast.matrix.sub),
            "logFC" = Inf,
            "SE" = NA,
            "Tvalue" = NA,
            "DF" = NA,
            "pvalue" = NA,
            "issue" = "oneConditionMissing")
        } else {
          out <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "Label" = row.names(contrast.matrix.sub),
            "logFC" = (-Inf),
            "SE" = NA,
            "Tvalue" = NA,
            "DF" = NA,
            "pvalue" = NA,
            "issue" = "oneConditionMissing")
        }
      } else {
        logging::logwarn(paste0("Results of protein ",
                                toString(unique(sub[["PROTEIN"]])),
                                " for comparison ",
                                rownames(contrast.matrix.sub),
                                " are NA because there are no measurements in both conditions."))
        out <- data.frame(
          "Protein" = unique(sub[["PROTEIN"]]),
          "Label" = row.names(contrast.matrix.sub),
          "logFC" = NA,
          "SE" = NA,
          "Tvalue" = NA,
          "DF" = NA,
          "pvalue" = NA,
          "issue" = "completeMissing")
      }
      allout <- rbind(allout, out)
    } ## end loop for comparion
    finalresid <- NULL
    finalfitted <- NULL
    fit.full <- NULL
  } else {
    ## when subject is fixed, it is ok using lm function.
    ## when single feature, consider technical replicates for time-course.
    ## case-control
    if (isTRUE(repeated)) {
      ## time-course
      if (singleSubject) {
        fit.full <- lm(ABUNDANCE ~ GROUP, data=sub)
      } else {
        ## no single subject
        if (!TechReplicate) {
          fit.full <- lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data=sub)
          df.full <- lm(ABUNDANCE ~ GROUP + SUBJECT, data=sub)$df.residual
        } else {
          fit.full <- lmer(ABUNDANCE ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT),
                           data=sub) ## SUBJECT==SUBJECT_NESTED here
          df.full <- lm(ABUNDANCE ~ GROUP + SUBJECT + GROUP:SUBJECT,
                        data=sub)$df.residual
        }
      }
    } else {
      if (!TechReplicate | singleSubject) {
        fit.full <- lm(ABUNDANCE ~ GROUP, data=sub)
      } else {
        fit.full <- lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data=sub)
        df.full <- lm(ABUNDANCE ~ GROUP + SUBJECT, data=sub)$df.residual
      }
    } ## time-course

    ## get parameter from model
    if (class(fit.full) == "lm") {
      Para <- getParameterFixed(fit.full)
    } else {
      Para <- getParameterRandom(fit.full, df.full)
    }

    ## each comparison
    allout <- NULL
    ## get condition IDs which are completely missing.
    emptycondition <- setdiff(levels(origGroup), unique(sub[["GROUP_ORIGINAL"]]))
    for (k in 1:nrow(contrast.matrix)) {
      ## choose each comparison
      contrast.matrix.sub <- matrix(contrast.matrix[k, ], nrow=1)
      row.names(contrast.matrix.sub) <- row.names(contrast.matrix)[k]
      colnames(contrast.matrix.sub) <- colnames(contrast.matrix)
      if (length(emptycondition) != 0) {
        ## if there are any completely missing in any condition,
        ## one by one comparison is simple. However, for linear combination of condition can be complicated
        ## get + and - condition separately
        count.pos <- levels(origGroup)[contrast.matrix.sub != 0 & contrast.matrix.sub > 0]
        count.neg <- levels(origGroup)[contrast.matrix.sub != 0 & contrast.matrix.sub < 0]
        ## then check whether any + or - part is completely missing
        count.diff.pos <- intersect(levels(origGroup)[contrast.matrix.sub != 0 &
                                                      contrast.matrix.sub > 0],
                                    emptycondition)
        count.diff.neg <- intersect(levels(origGroup)[contrast.matrix.sub != 0 &
                                                      contrast.matrix.sub < 0],
                                    emptycondition)

        ## positive side
        if (length(count.diff.pos) != 0) {
          flag.issue.pos <- TRUE ## TRUE : there are problematic conditions
        } else {
          flag.issue.pos <- FALSE ## FALSE : no issue about completely missing
        }
        ## negative side
        if (length(count.diff.neg) != 0) {
          flag.issue.neg <- TRUE ## TRUE : there are problematic conditions
        } else {
          flag.issue.neg <- FALSE ## FALSE : no issue about completely missing
        }

        ## message and output
        if (flag.issue.pos & flag.issue.neg) {
          ## both sides are completely missing
          logging::logwarn(paste0("Results of protein ",
                                  unique(sub[["PROTEIN"]]),
                                  " for comparison ", rownames(contrast.matrix.sub),
                                  " are NA because there are no measurements in both conditions."))
          out <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "Label" = row.names(contrast.matrix.sub),
            "logFC" = NA,
            "SE" = NA,
            "Tvalue" = NA,
            "DF" = NA,
            "pvalue" = NA,
            "issue" = "completeMissing")
        } else if (flag.issue.pos | flag.issue.neg) {
          if (flag.issue.pos) {
            issue.side <- count.diff.pos
            logging::logwarn(paste0("Results of protein ",
                                    unique(sub[["PROTEIN"]]),
                                    " for comparison ",
                                    rownames(contrast.matrix.sub),
                                    " are NA because there are measurements only in group ",
                                    toString(issue.side)))
            out <- data.frame(
              "Protein" = unique(sub[["PROTEIN"]]),
              "Label" = rownames(contrast.matrix.sub),
              "logFC" = (-Inf),
              "SE" = NA,
              "Tvalue" = NA,
              "DF" = NA,
              "pvalue" = NA,
              "issue" = "oneConditionMissing")
          }

          if (flag.issue.neg) {
            issue.side <- count.diff.neg
            logging::logwarn(paste0("Results of protein ",
                                    unique(sub[["PROTEIN"]]),
                                    " for comparison ",
                                    rownames(contrast.matrix.sub),
                                    " are NA because there are measurements only in group ",
                                    toString(issue.side)))
            out <- data.frame(
              "Protein" = unique(sub[["PROTEIN"]]),
              "Label" = rownames(contrast.matrix.sub),
              "logFC" = Inf,
              "SE" = NA,
              "Tvalue" = NA,
              "DF" = NA,
              "pvalue" = NA,
              "issue" = "oneConditionMissing")
          }
        } else { # then same as regular calculation
          contrast <- make.contrast.free.single(fit.full, contrast.matrix.sub, sub)
          out <- estimableFixedRandom(Para, contrast)
          ## any error for out, just NA
          if (is.null(out)) {
            out <- data.frame(
              "Protein" = unique(sub[["PROTEIN"]]),
              "Label" = row.names(contrast.matrix.sub),
              "logFC" = NA,
              "SE" = NA,
              "Tvalue" = NA,
              "DF" = NA,
              "pvalue" = NA,
              "issue" = NA)
          } else {
            out <- data.frame(
              "Protein" = unique(sub[["PROTEIN"]]),
              "Label" = row.names(contrast.matrix.sub),
              out,
              "issue" = NA)
          }
        }
      } else {
        contrast <- make.contrast.free.single(fit.full, contrast.matrix.sub, sub)
        out <- estimableFixedRandom(Para, contrast)
        ## any error for out, just NA
        if (is.null(out)) {
          out <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "Label" = row.names(contrast.matrix.sub),
            "logFC" = NA,
            "SE" = NA,
            "Tvalue" = NA,
            "DF" = NA,
            "pvalue" = NA,
            "issue" = NA)
        } else {
          out <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "Label" = row.names(contrast.matrix.sub),
            out,
            "issue" = NA)
        }
      }
      allout <- rbind(allout, out)
    } ## end loop for comparion

    if (class(fit.full) == "lm") { ## lm model
      finalresid <- fit.full$residuals
      finalfitted <- fit.full$fitted.values
    } else {  ## lmer model
      finalresid <- resid(fit.full)
      finalfitted <- fitted(fit.full)
    }
  } ## more than 2 conditions in the dataset

  finalout <- list(
    "result" = allout,
    "valueresid" = finalresid,
    "valuefitted" = finalfitted,
    "fittedmodel" = fit.full)

  return(finalout)
}

#############################################
ttest.logsum <- function(contrast.matrix, data, origGroup) {
  ## each comparison
  allout <- NULL
  for (k in 1:nrow(contrast.matrix)) {
    ## choose each comparison
    contrast.matrix.sub <- matrix(contrast.matrix[k, ], nrow=1)
    row.names(contrast.matrix.sub) <- row.names(contrast.matrix)[k]
    ##GroupComparisonAgreement <- checkGroupComparisonAgreement(data,contrast.matrix.sub)
    ##    if (GroupComparisonAgreement$sign==TRUE) {
    ##     message("*** error : results of Protein ", unique(data$PROTEIN), " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")
    ##     out <- data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
    ##    }else{

    ## get two groups from contrast.matrix
    datasub <- data[which(data$GROUP_ORIGINAL %in% origGroup[contrast.matrix.sub != 0]), ]
    ## t test
    sumresult <- try(t.test(datasub$ABUNDANCE ~ datasub$GROUP_ORIGINAL,
                            var.equal=TRUE), silent=TRUE)
    if (class(sumresult) == "try-error") {
      out <- data.frame(
        "Protein" = unique(data$PROTEIN),
        "Label" = row.names(contrast.matrix.sub),
        "logFC" = NA,
        "SE" = NA,
        "Tvalue" = NA,
        "DF" = NA,
        "pvalue" = NA,
        "issue" = NA)
    } else {
      out <- data.frame(
        "Protein" = unique(data$PROTEIN),
        "Label" = paste0(names(sumresult$estimate)[1], " - ", names(sumresult$estimate)[2]),
        "logFC" = sumresult$estimate[1] - sumresult$estimate[2],
        "SE" = (sumresult$estimate[1] - sumresult$estimate[2]) / sumresult$statistic,
        "Tvalue" = sumresult$statistic,
        "DF" = sumresult$parameter,
        "pvalue" = sumresult$p.value,
        "issue" = NA)
      rownames(out) <- NULL
    }
    ## }
    allout <- rbind(allout, out)
  } # end loop for comparion

  finalout <- list(
    "result" = allout,
    "valueresid" = NULL,
    "valuefitted" = NULL,
    "fittedmodel" = NULL)
  return(finalout)
}

#############################################
count.missing.percentage <- function(contrast.matrix, temptempresult, sub, subprocess) {
  totaln.cell <- aggregate(PROTEIN ~ GROUP_ORIGINAL,
                           data=subprocess[subprocess$LABEL == "L", ], length)
  ## just for count total measurement, use PROTEIN instead of ABUNDANCE in order to prevent to remove NA in ABUNDANCE
  colnames(totaln.cell)[colnames(totaln.cell) == "PROTEIN"] <- "totalN"
  totaln.measured <- aggregate(NumMeasuredFeature ~ GROUP_ORIGINAL, data=sub, sum, na.rm=TRUE)
  totaln.imputed <- aggregate(NumImputedFeature ~ GROUP_ORIGINAL, data=sub, sum, na.rm=TRUE)
  totaln <- merge(totaln.cell, totaln.measured, by="GROUP_ORIGINAL", all=TRUE)
  totaln <- merge(totaln, totaln.imputed, by="GROUP_ORIGINAL", all=TRUE)

  if (any(is.element(colnames(sub), "NumImputedFeature"))) {
    temptempresult$MissingPercentage <- NA
    temptempresult$ImputationPercentage <- NA
  } else {
    temptempresult$MissingPercentage <- NA
  }
  for (k in 1:nrow(contrast.matrix)) {
    ## choose each comparison
    contrast.matrix.sub <- matrix(contrast.matrix[k, ], nrow=1)
    row.names(contrast.matrix.sub) <- row.names(contrast.matrix)[k]
    condition.needed <- contrast.matrix.sub != 0
    MissingPercentage.new <- NA
    ImputationPercentage.new <- NA
    ## total # missing
    MissingPercentage.new <- 1 -
      (sum(totaln[condition.needed, "NumMeasuredFeature"], na.rm=TRUE) /
       sum(totaln[condition.needed, "totalN"], na.rm=TRUE))

    ## # imputed intensity
    if (any(is.element(colnames(sub), "NumImputedFeature"))) {
      ImputationPercentage.new <- sum(
        totaln[condition.needed, "NumImputedFeature"], na.rm=TRUE) /
        sum(totaln[condition.needed, "totalN"], na.rm=TRUE)
      temptempresult[temptempresult$Label == row.names(contrast.matrix.sub),
                     "MissingPercentage"] <- MissingPercentage.new
      temptempresult[temptempresult$Label == row.names(contrast.matrix.sub),
                     "ImputationPercentage"] <- ImputationPercentage.new
    } else {
      temptempresult[temptempresult$Label == row.names(contrast.matrix.sub),
                     "MissingPercentage"] <- MissingPercentage.new
    }
  } # end loop for multiple comparisons
  return(temptempresult)
}
