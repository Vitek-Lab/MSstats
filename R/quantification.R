#' Protein sample quantification or group quantification
#'
#' @description Model-based quantification for each condition or for each biological
#' sample per protein in a targeted Selected Reaction Monitoring (SRM),
#' Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition
#' (DIA or SWATH-MS) experiment. Quantification takes the processed data set
#' by \code{\link{dataProcess}} as input and automatically generate the quantification
#' results (data.frame) in a long or matrix format.
#'
#' @param data name of the (processed) data set.
#' @param type choice of quantification. "Sample" or "Group" for protein sample
#' quantification or group quantification.
#' @param format choice of returned format. "long" for long format which has
#' the columns named Protein, Condition, LogIntensities (and BioReplicate if it is
#' subject quantification), NumFeature for number of transitions for a protein,
#' and NumPeaks for number of observed peak intensities for a protein.
#' "matrix" for data matrix format which has the rows for Protein and the columns,
#' which are Groups(or Conditions) for group quantification or the combinations
#' of BioReplicate and Condition (labeled by "BioReplicate"_"Condition")
#' for sample quantification. Default is "matrix"
#' @inheritParams .documentFunction
#'
#' @details
#' \itemize{
#' \item{Sample quantification : individual biological sample quantification for each protein. The label of each biological sample is a combination of the corresponding group and the sample ID. If there are no technical replicates or experimental replicates per sample, sample quantification is the same as run summarization from dataProcess. If there are technical replicates or experimental replicates, sample quantification is median among run quantification corresponding MS runs.}
#' \item{Group quantification : quantification for individual group or individual condition per protein. It is median among sample quantification.}
#' \item{The quantification for endogenous samples is based on run summarization from subplot model, with TMP robust estimation.}
#' }
#'
#' @return data.frame as described in details.
#'
#' @importFrom stats median
#'
#' @export
#'
#' @examples
#' # Consider quantitative data (i.e. QuantData) from a yeast study with ten time points of
#' # interests, three biological replicates, and no technical replicates which is
#' # a time-course experiment.
#' # Sample quantification shows model-based estimation of protein abundance in each biological
#' # replicate within each time point.
#' # Group quantification shows model-based estimation of protein abundance in each time point.
#' QuantData<-dataProcess(SRMRawData, use_log_file = FALSE)
#' head(QuantData$FeatureLevelData)
#' # Sample quantification
#' sampleQuant<-quantification(QuantData, use_log_file = FALSE)
#' head(sampleQuant)
#' # Group quantification
#' groupQuant<-quantification(QuantData, type="Group", use_log_file = FALSE)
#' head(groupQuant)
#'
quantification = function(
    data, type = "Sample", format="matrix", use_log_file = TRUE, append = FALSE,
    verbose = TRUE, log_file_path = NULL
) {
    LogIntensities = NULL
    
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose,
                                        log_file_path, base = "MSstats_quant_log_")
    getOption("MSstatsLog")("INFO", "MSstats - quantification function")
    checkmate::assertChoice(toupper(type), c("SAMPLE", "GROUP"),
                            .var.name = "type")
    checkmate::assertChoice(toupper(format), c("MATRIX", "LONG"),
                            .var.name = "format")
    getOption("MSstatsLog")("INFO", paste0("type of quantification = ", type))
    getOption("MSstatsLog")("INFO", paste0("format of output = ", format))

    if (toupper(type) == "SAMPLE") {
        datarun = data.table::as.data.table(data$ProteinLevelData)
        datarun$Protein = factor(datarun$Protein)
        datarun = datarun[!is.na(LogIntensities), ]
        datam = data.table::dcast(Protein ~ GROUP + SUBJECT, 
                                  data = datarun,
                                  value.var = 'LogIntensities',
                                  fun.aggregate = median)
        if (format == "long") {
            data_l = melt(datam, id.vars=c('Protein'))
            colnames(data_l)[colnames(data_l) %in% c("variable", "value")] = c('Group_Subject', 'LogIntensity')
        }
        getOption("MSstatsLog")("INFO", "Finish sample quantificiation - okay.")
        if (format == "long") {
            return(data_l)
        }
        if (format == "matrix") {
            return(datam)
        }
    }

    if (toupper(type) == "GROUP") {
        datarun = data.table::as.data.table(data$ProteinLevelData)
        datarun$Protein = factor(datarun$Protein)
        datarun = datarun[!is.na(LogIntensities), ]
        datam = data.table::dcast(Protein + GROUP ~ SUBJECT,
                      data = datarun,
                      value.var = 'LogIntensities',
                      fun.aggregate = median)
        datam2 = data.table::melt(datam, id.vars = c('Protein', "GROUP"))
        data.table::setnames(datam2, c("variable", "value"),
                             c("Subject", "LogIntensity"))
        datam3 = data.table::dcast(Protein ~ GROUP,
                       data = datam2,
                       value.var = 'LogIntensity',
                       fun.aggregate = function(x) median(x, na.rm=TRUE))
        if (format == "long") {
            data_l = melt(datam3, id.vars=c('Protein'))
            data.table::setnames(data_l, c("variable", "value"),
                                 c("Group", "LogIntensity"))
        }
        getOption("MSstatsLog")("INFO", "Finish group quantificiation - okay.")
        if (format == "long") {
            return(data_l)
        }
        if (format == "matrix") {
            return(datam3)
        }
    }
}
