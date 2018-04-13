#' My voice is my password, quantify me.
#'
#' No seriously, this function performs quantifications on a MSMS data set.
#'
#' @param data  Input data from processData()
#' @param type  Quantify based on which column in the experimental design?
#' @param format  Return data as a ... matrix?
#' @return the quantified data!
#' @export
quantification <- function(data, type="Sample", format="matrix") {
  logging::loginfo("Starting quantification().")
  ## data format
  ## rawinput <-
  ## c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge",
                                        # "IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
  ## if (length(setdiff(toupper(rawinput),toupper(colnames(data))))==0) {
  ## }
  type <- toupper(type)
  format <- toupper(format)
  ## check other options
  if (!(type == "SAMPLE" | type == "GROUP")) {
    print_string <- "The required input 'type' is wrong or missing, setting it to 'SAMPLE'."
    type <- "SAMPLE"
    logging::logwarn(print_string)
  }
  if (!(format == "MATRIX" | format == "LONG")) {
    print_string <- "The required input 'format' is wrong or missing, setting it to 'MATRIX'."
    type <- "MATRIX"
    logging::logwarn(print_string)
  }

  ## all input
  logging::loginfo(paste0("The type of quantification: ", type, "."))
  logging::loginfo(paste0("The format of the output: ", format, "."))

#################################
### sample quantification
#################################
  if (type == "SAMPLE") {
    datarun <- data[["RunlevelData"]]
    datarun[["Protein"]] <- factor(datarun[["Protein"]])
    datarun <- datarun[!is.na(datarun[["LogIntensities"]]), ]
    datam <- reshape2::dcast(Protein ~ GROUP_ORIGINAL + SUBJECT_ORIGINAL,
                             data=datarun,
                             value.var="LogIntensities",
                             fun.aggregate=median)
    if (format=="LONG") {
      data_l <- reshape2::melt(datam, id.vars=c("Protein"))
      colnames(data_l)[colnames(data_l) %in%
                       c("variable", "value")] <- c("Group_Subject", "LogIntensity")
    }

    logging::loginfo("Finished sample quantification.")
    if (format == "LONG") {
      return(data_l)
    }
    if (format == "MATRIX") {
      return(datam)
    }
  } ## end sample quantification

#################################
### Group quantification
#################################
  if (type == "GROUP") {
    datarun <- data[["RunlevelData"]]
    datarun[["Protein"]] <- factor(datarun[["Protein"]])
    datarun <- datarun[!is.na(datarun[["LogIntensities"]]), ]
    datam <- reshape2::dcast(Protein + GROUP_ORIGINAL ~ SUBJECT_ORIGINAL,
                             data=datarun,
                             value.var="LogIntensities",
                             fun.aggregate=median)
    datam2 <- reshape2::melt(datam, id.vars=c("Protein", "GROUP_ORIGINAL"))
    colnames(datam2)[colnames(datam2) %in% c("variable", "value")] <- c("Subject", "LogIntensity")
    datam3 <- reshape2::dcast(Protein ~ GROUP_ORIGINAL,
                              data=datam2,
                              value.var="LogIntensity",
                              fun.aggregate=function(x) {
                                median(x, na.rm=TRUE)
                              })
    rm(datam)
    rm(datam2)

    if (format == "LONG") {
      data_l <- melt(datam3, id.vars=c("Protein"))
      colnames(data_l)[colnames(data_l) %in% c("variable", "value")] <- c("Group", "LogIntensity")
    }

    logging::loginfo("Finished group quantification.")
    if (format == "LONG") {
      return(data_l)
    }
    if (format == "MATRIX") {
      return(datam3)
    }
  } ## end group quantification
}
