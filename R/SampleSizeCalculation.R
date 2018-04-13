#' Perform a statistical power calculation given a dataset, experimental design,
#' and some parameters defining the desired parameters for analysis.
#'
#' @param data The dataset to examine, presumably after dataProcess()
#' @param desiredFC  The fold-change used to define 'significant'
#' @param FDR  The false discovery rate used to define 'significant'
#' @param numSample  Use the number of samples in the data to inform the power
#'   calculation?
#' @param power  Floating point number to describe the desired power in the
#'   data.  This is strange, as I only see it used in a if(isTRUE()), but it is
#'   defined here as 0.9... hmmm I need to read this code.
#' @param scopeofBioReplication The following parameters were just defined right
#'   after the function definition, so I put them into it, I dunno what they do
#'   yet.
#' @param interference  Personally, I hate it when people interfere with me, but
#'   in this case, TRUE!
#' @param equalFeatureVar  This is a parameter from dataProcess, I know that much.
#' @return Dataframe describing the power in the data.
#' @export
designSampleSize <- function(data, desiredFC=2.0, FDR=0.05, numSample=TRUE,
                             power=0.9, labeled=FALSE,
                             scopeofBioReplication="expanded",
                             interference=TRUE, equalFeatureVar=TRUE) {

  print_string <- paste0("MSstats - designSampleSize: desired fold change = ",
                         desiredFC, ", FDR = ", FDR, ", Power = ", power, ".")
  logging::loginfo(print_string)

  ## for label-free experiment
  ##if (!labeled) {
  sigma.error <- NULL
  VarComponent <- data.frame(Protein=seq(1, length(data)), Error=NA, Subject=NA, GroupBySubject=NA)
  for (i in 1:length(data)) {
    ## note: when run is fixed, we can obtain the same variance of error for
    ## both case-control and time course studies.
    fit.full <- data[[i]]
    ## if fit.full==NA (class(fit.full)=="try-error)
    if (is.null(fit.full)) {
      ## !!!!! but if we have NULL for last protein?
      next
    } else {
      ## get variance component
      if (class(fit.full) != "mer") {
        VarComponent[i, "Error"] <- summary(fit.full)$sigma ^ 2
      } else {
        stddev <- c(sapply(VarCorr(fit.full),
                           function(el) attr(el, "stddev")), attr(VarCorr(fit.full), "sc"))
        VarComponent[i, "Error"] <- stddev[names(stddev) == ""] ^ 2
        if (sum(names(stddev) %in% "SUBJECT_NESTED.(Intercept)") > 0) {
          VarComponent[i, "Subject"] <- stddev[names(stddev) == "SUBJECT_NESTED.(Intercept)"] ^ 2
        }
        if (sum(names(stddev) %in% "SUBJECT.(Intercept)") > 0) {
          VarComponent[i, "Subject"] <- stddev[names(stddev) == "SUBJECT.(Intercept)"] ^ 2
        }
        if (sum(names(stddev) %in% "SUBJECT:GROUP.(Intercept)") > 0) {
          VarComponent[i, "GroupBySubject"] <-
            stddev[names(stddev) == "SUBJECT:GROUP.(Intercept)"] ^ 2
        }
      }
    }
  } ## end-loop
  ##        VarComponent[is.na(VarComponent)] <- 0
  ## for label-free DDA, there are lots of missingness and lots of zero SE. So, remove NA SE.
  median.sigma.error <- median(VarComponent[, "Error"], na.rm=TRUE)
  if (sum(!is.na(VarComponent[, "GroupBySubject"])) > 0) {
    median.sigma.subject <- median(VarComponent[, "GroupBySubject"], na.rm=TRUE)
  } else {
    if (sum(!is.na(VarComponent[, "Subject"])) > 0) {
      median.sigma.subject <- median(VarComponent[, "Subject"], na.rm=TRUE)
    } else {
      median.sigma.subject <- 0
    }
  }
  ##
  print_string <- "designSampleSize(): Calculated variance component."
  logging::loginfo(print_string)

  ## power calculation
  if (isTRUE(power)) {
    delta <- log2(seq(desiredFC[1], desiredFC[2], 0.025))
    desiredFC <- 2 ^ delta
    m0_m1 <- 99
    t <- delta / sqrt(2 * (median.sigma.error / numSample + median.sigma.subject / numSample))
    ##t <- delta/sqrt(2*median.sigma.error/numPep/numTran/numSample+median.sigma.subject/numSample)
    powerTemp <- seq(0, 1, 0.01)

    power <- NULL
    for (i in 1:length(t)) {
      diff <- qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
      min(abs(diff), na.rm=TRUE)
      power[i] <- powerTemp[order(abs(diff))][1]
    }

    CV <- round((2 * (median.sigma.error / numSample +
                      median.sigma.subject / numSample)) / desiredFC, 3)

    print_string <- "designSampleSize(): Power has been calculated."
    logging::loginfo(print_string)
    out <- data.frame(desiredFC, numSample, FDR, power=power, CV)
    return(out)
  }

  if (is.numeric(power)) {
    ## Large portion of proteins are not changing
    m0_m1 <- 99 ## it means m0/m1=99, m0/(m0+m1)=0.99
    alpha <- power * FDR / (1 + (1 - FDR) * m0_m1)

    ## Num Sample calculation
    if (isTRUE(numSample)) {
      delta <- log2(seq(desiredFC[1], desiredFC[2], 0.025))
      desiredFC <- 2 ^ delta
      z_alpha <- qnorm(1 - alpha / 2)
      z_beta <- qnorm(power)
      aa <- (delta / (z_alpha + z_beta)) ^ 2
      numSample <- round(2 * (median.sigma.error + median.sigma.subject) / aa, 0)
      CV <- round(2 *
                  (median.sigma.error / numSample +
                   median.sigma.subject / numSample) / desiredFC, 3)

      print_string <- "designSampleSize(): The number of samples has been calculated."
      logging::loginfo(print_string)
      out <- data.frame(desiredFC, numSample, FDR, power, CV)
      return(out)
    }
  } # when power is numeric
  ##} ## label-free
}

#' A couple plots generated when looking at sample sizes of a msstats dataset.
#'
#' @param data  The dataProcess()'d data structure.
#' @param ...  Extra arguments for ggplot2.
#' @return plots!
#' @export
designSampleSizePlots <- function(data, ...) {
  arglist <- list(...)
  if (length(unique(data[["numSample"]])) > 1) {
    index <- "numSample"
  }
  if (length(unique(data[["power"]])) > 1) {
    index <- "power"
  }
  if (length(unique(data[["numSample"]])) == 1 &
      length(unique(data[["power"]])) == 1) {
    index <- "numSample"
  }

  text.size <- 1.2
  if (!is.null(arglist[["text.size"]])) {
    text.size <- arglist[["text.size"]]
  }
  axis.size <- 1.3
  if (!is.null(arglist[["axis.size"]])) {
    axis.size <- arglist[["axis.size"]]
  }
  lab.size <- 1.7
  if (!is.null(arglist[["lab.size"]])) {
    lab.size <- arglist[["lab.size"]]
  }

  retlist <- list()
  if (index == "numSample") {
    with(data, {
      plot(desiredFC, numSample, lwd=2, xlab="", ylab="", cex.axis=axis.size, type="l", xaxt="n")})
    axis(1, at=seq(min(data[["desiredFC"]]), max(data[["desiredFC"]]), 0.05),
         labels=seq(min(data[["desiredFC"]]), max(data[["desiredFC"]]), 0.05),
         cex.axis=axis.size)
    axis(3, at=seq(min(data[["desiredFC"]]), max(data[["desiredFC"]]), 0.05),
         labels=data[["CV"]][which(data[["desiredFC"]] %in%
                              seq(min(data[["desiredFC"]]), max(data[["desiredFC"]]), 0.05))],
         cex.axis=axis.size)
    mtext("Coefficient of variation, CV", 3, line=2.5, cex=lab.size)
    mtext("Desired fold change", 1, line=3.5, cex=lab.size)
    mtext("Minimal number of biological replicates", 2, line=2.5, cex=lab.size)
    legend("topright", c(paste("FDR is", unique(data$FDR), sep=" "),
                         paste("Statistical power is", unique(data[["power"]]), sep=" ")),
           bty="n", cex=text.size)
    retlist[["num_samples"]] <- grDevices::recordPlot()
  }

  if (index == "power") {
    with(data, {
      plot(desiredFC, power, lwd=2, xlab="", ylab="", cex.axis=axis.size, type="l", xaxt="n")})
    axis(1, at=seq(min(data[["desiredFC"]]), max(data[["desiredFC"]]), 0.05),
         labels=seq(min(data[["desiredFC"]]), max(data[["desiredFC"]]), 0.05),
         cex.axis=axis.size)
    mtext("Desired fold change", 1, line=3.5, cex=lab.size)
    mtext("Power", 2, line=2.5, cex=lab.size)
    legend("bottomright", c(paste("Number of replicates is", unique(data[["numSample"]]), sep=" "),
                            paste("FDR is", unique(data[["FDR"]]), sep=" ")), bty="n",
           cex=text.size)
    retlist[["power"]] <- grDevices::recordPlot()
  }

  ## Return a list of the plots created so that one may save/print/whatever them later.
  return(retlist)
}
