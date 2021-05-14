#' Planning future experimental designs of Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiments in sample size calculation
#'
#' @description Calculate sample size for future experiments of a Selected Reaction Monitoring (SRM), 
#' Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment 
#' based on intensity-based linear model. Two options of the calculation: 
#' (1) number of biological replicates per condition, 
#' (2) power.
#' 
#' @param data 'FittedModel' in testing output from function groupComparison.
#' @param desiredFC the range of a desired fold change which includes the lower 
#' and upper values of the desired fold change.
#' @param FDR a pre-specified false discovery ratio (FDR) to control the overall 
#' false positive rate. Default is 0.05
#' @param numSample minimal number of biological replicates per condition. 
#' TRUE represents you require to calculate the sample size for this category, 
#' else you should input the exact number of biological replicates.
#' @param power a pre-specified statistical power which defined as the probability 
#' of detecting a true fold change. TRUE represent you require to calculate the power 
#' for this category, else you should input the average of power you expect. Default is 0.9
#' @inheritParams .documentFunction
#' 
#' @details The function fits the model and uses variance components to calculate 
#' sample size. The underlying model fitting with intensity-based linear model with 
#' technical MS run replication. Estimated sample size is rounded to 0 decimal.
#' The function can only obtain either one of the categories of the sample size 
#' calculation (numSample, numPep, numTran, power) at the same time.
#' 
#' @return data.frame - sample size calculation results including varibles:
#' desiredFC, numSample, FDR,  and power.
#' 
#' @importFrom stats median
#' 
#' @export
#'  
#' @author Meena Choi, Ching-Yun Chang, Olga Vitek. 
#' 
#' @examples
#' # Consider quantitative data (i.e. QuantData) from yeast study.
#' # A time course study with ten time points of interests and three biological replicates.
#' QuantData <- dataProcess(SRMRawData)
#' head(QuantData$FeatureLevelData)
#' ## based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
#' comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
#' comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
#' comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
#' comparison<-rbind(comparison1,comparison2, comparison3)
#' row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")
#' colnames(comparison)<-unique(QuantData$ProteinLevelData$GROUP)
#' 
#' testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=QuantData)
#' 
#' ## Calculate sample size for future experiments:
#' #(1) Minimal number of biological replicates per condition
#' designSampleSize(data=testResultMultiComparisons$FittedModel, numSample=TRUE,
#'                  desiredFC=c(1.25,1.75), FDR=0.05, power=0.8)
#' #(2) Power calculation
#' designSampleSize(data=testResultMultiComparisons$FittedModel, numSample=2,
#'                  desiredFC=c(1.25,1.75), FDR=0.05, power=TRUE)
#'                  
designSampleSize = function(
    data, desiredFC, FDR = 0.05, numSample = TRUE, power = 0.9,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, "MSstats_sampleSize_log_")
    getOption("MSstatsLog")("INFO", "** MSstats - designSampleSize function")
    getOption("MSstatsLog")("INFO", paste0("Desired fold change = ", 
                                           paste(desiredFC, collapse=" - ")))
    getOption("MSstatsLog")("INFO", paste0("FDR = ", FDR))
    getOption("MSstatsLog")("INFO", paste0("Power = ", power))
    
    var_component = .getVarComponent(data)
    ## for label-free DDA, there are lots of missingness and lots of zero SE. So, remove NA SE.
    median_sigma_error = median(var_component[["Error"]], na.rm = TRUE)
    median_sigma_subject = .getMedianSigmaSubject(var_component)
    getOption("MSstatsLog")("INFO", "Calculated variance component. - okay")

    ## power calculation
    if (isTRUE(power)) {
        delta = log2(seq(desiredFC[1], desiredFC[2], 0.025))
        desiredFC = 2 ^ delta
        power_output = .calculatePower(desiredFC, FDR, delta, median_sigma_error, 
                                       median_sigma_subject, numSample)        
        CV = round( (2 * (median_sigma_error / numSample + median_sigma_subject / numSample)) / desiredFC, 3)
        getOption("MSstatsLog")("INFO", "Power is calculated. - okay")
        sample_size = data.frame(desiredFC, numSample, FDR, 
                                 power = power_output, CV)
    }	
    
    if (is.numeric(power)) {
        delta = log2(seq(desiredFC[1], desiredFC[2], 0.025))
        desiredFC = 2 ^ delta
        ## Large portion of proteins are not changing
        m0_m1 = 99 ## it means m0/m1=99, m0/(m0+m1)=0.99
        alpha = power * FDR / (1 + (1 - FDR) * m0_m1)
        if (isTRUE(numSample)) {
            numSample = .getNumSample(desiredFC, power, alpha, delta,
                                      median_sigma_error, median_sigma_subject)
            CV = round(2 * (median_sigma_error / numSample + median_sigma_subject / numSample) / desiredFC, 3)
            getOption("MSstatsLog")("INFO", "The number of sample is calculated. - okay")
            sample_size = data.frame(desiredFC, numSample, FDR, power, CV)
        }
    } 
    sample_size
}


#' Get variances from models fitted by the groupComparison function
#' @param fitted_models FittedModels element of groupComparison output
#' @keywords internal
.getVarComponent = function(fitted_models) {
    Protein = NULL
    
    result = data.table::rbindlist(
        lapply(fitted_models, function(fit) {
            if (!is.null(fit)) {
                if (!is(fit, "lmerMod")) {
                    error = summary(fit)$sigma^2
                    subject <- NA
                    group_subject <- NA
                } else {
                    stddev = c(sapply(lme4::VarCorr(fit), function(el) attr(el, "stddev")), 
                               attr(lme4::VarCorr(fit), "sc"))
                    error = stddev[names(stddev) == ""]^2
                    if (any(names(stddev) %in% "SUBJECT.(Intercept)")) {
                        subject = stddev["SUBJECT.(Intercept)"]^2
                    } else {
                        subject = NA
                    }
                    if (any(names(stddev) %in% "SUBJECT:GROUP.(Intercept)")) {
                        group_subject = stddev["SUBJECT:GROUP.(Intercept)"]^2
                    } else {
                        group_subject = NA
                    }
                }
                list(Error = error,
                     Subject = subject,
                     GroupBySubject = group_subject)
            } else {
                NULL
            }
        })
    )
    result[, Protein := 1:.N]
    result
}


#' Get median per subject or group by subject
#' @param var_component data.frame, output of .getVarComponent
#' @importFrom stats median
#' @keywords internal
.getMedianSigmaSubject = function(var_component) {
    if (sum(!is.na(var_component[, "GroupBySubject"])) > 0) {
        median_sigma_subject = median(var_component[["GroupBySubject"]], na.rm=TRUE)
    } else {
        if (sum(!is.na(var_component[, "Subject"])) > 0) {
            median_sigma_subject = median(var_component[["Subject"]], na.rm=TRUE)
        } else {
            median_sigma_subject = 0
        }
    }
    median_sigma_subject
}


#' Power calculation
#' @inheritParams designSampleSize
#' @param delta difference between means (?)
#' @param median_sigma_error median of error standard deviation
#' @param median_sigma_subject median standard deviation per subject
#' @importFrom stats qnorm
#' @keywords internal
.calculatePower = function(desiredFC, FDR, delta, median_sigma_error, 
                           median_sigma_subject, numSample) {
    m0_m1 = 99
    t = delta / sqrt(2 * (median_sigma_error / numSample + median_sigma_subject / numSample))
    powerTemp = seq(0, 1, 0.01)
    power = numeric(length(t))
    for (i in seq_along(t)) {
        diff = qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
        min(abs(diff), na.rm = TRUE)
        power[i] = powerTemp[order(abs(diff))][1]
    }
    power
}


#' Get sample size
#' @inheritParams designSampleSize
#' @inheritParams .calculatePower
#' @param alpha significance level
#' @param delta difference between means (?)
#' @importFrom stats qnorm
#' @keywords internal
.getNumSample = function(desiredFC, power, alpha, delta, median_sigma_error, 
                         median_sigma_subject){
    z_alpha = qnorm(1 - alpha / 2)
    z_beta = qnorm(power)
    aa = (delta / (z_alpha + z_beta)) ^ 2
    numSample = round(2 * (median_sigma_error + median_sigma_subject) / aa, 0)
    numSample
}


#' Visualization for sample size calculation
#' 
#' @description To illustrate the relationship of desired fold change and the calculated 
#' minimal number sample size which are (1) number of biological replicates per condition, 
#' (2) number of peptides per protein, 
#' (3) number of transitions per peptide, and 
#' (4) power. The input is the result from function (\code{\link{designSampleSize}}.
#' 
#' @param data output from function designSampleSize.
#' 
#' @details Data in the example is based on the results of sample size calculation from function \code{\link{designSampleSize}}
#' 
#' @return Plot for estimated sample size with assigned variable.
#' 
#' @export
#' 
#' @author Meena Choi, Ching-Yun Chang, Olga Vitek. 
#' @examples
#' # Based on the results of sample size calculation from function designSampleSize,
#' # we generate a series of sample size plots for number of biological replicates, or peptides, 
#' # or transitions or power plot.
#' QuantData<-dataProcess(SRMRawData)
#' head(QuantData$ProcessedData)
#' ## based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
#' comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
#' comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
#' comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
#' comparison<-rbind(comparison1,comparison2, comparison3)
#' row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")
#' colnames(comparison)<-unique(QuantData$ProteinLevelData$GROUP)
#' 
#' testResultMultiComparisons<-groupComparison(contrast.matrix=comparison, data=QuantData)
#' 
#' # plot the calculated sample sizes for future experiments:
#' # (1) Minimal number of biological replicates per condition
#' result.sample<-designSampleSize(data=testResultMultiComparisons$FittedModel, numSample=TRUE,
#'                                 desiredFC=c(1.25,1.75), FDR=0.05, power=0.8)
#' designSampleSizePlots(data=result.sample)
#' # (2) Power
#' result.power<-designSampleSize(data=testResultMultiComparisons$FittedModel, numSample=2,
#'                                desiredFC=c(1.25,1.75), FDR=0.05, power=TRUE)
#' designSampleSizePlots(data=result.power)
#' 
designSampleSizePlots = function(data) {
    if (length(unique(data$numSample)) > 1) {
        index = "numSample"
    } else if (length(unique(data$power)) > 1) {
        index = "power"
    } else if (length(unique(data$numSample)) == 1 & length(unique(data$power)) == 1) {
        index = "numSample"
    } else {
        stop ("Invalid input")
    }
    
    text.size = 1.2	
    axis.size = 1.3	
    lab.size = 1.7
    
    if (index == "numSample") {
        plot(data$desiredFC, data$numSample, 
             lwd=2, xlab="", ylab="", 
             cex.axis=axis.size, type="l", xaxt="n")
        axis(1, at=seq(min(data$desiredFC), max(data$desiredFC), 0.05), 
             labels=seq(min(data$desiredFC), max(data$desiredFC), 0.05),
             cex.axis=axis.size)
        mtext("Desired fold change", 1, line=3.5, cex=lab.size)
        mtext("Minimal number of biological replicates", 2, line=2.5, cex=lab.size)
        legend("topright",
               c(paste("FDR is", unique(data$FDR)),
                 paste("Statistical power is", unique(data$power))),
               bty="n",
               cex=text.size)
    }
    
    if (index == "power") {
        plot(data$desiredFC, data$power, 
             lwd=2, xlab="", ylab="", 
             cex.axis=axis.size, type="l", xaxt="n")
        axis(1, at=seq(min(data$desiredFC), max(data$desiredFC), 0.05), 
             labels=seq(min(data$desiredFC), max(data$desiredFC), 0.05), 
             cex.axis=axis.size)
        mtext("Desired fold change", 1, line=3.5, cex=lab.size)
        mtext("Power", 2, line=2.5, cex=lab.size)
        legend("bottomright", 
               c(paste("Number of replicates is", unique(data$numSample)), 
                 paste("FDR is", unique(data$FDR))), 
               bty="n", 
               cex=text.size)
    }
}
