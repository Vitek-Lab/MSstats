#######################################
## MSnSet -> input for MSstats
#######################################

#' Transform a MSNset into a MSstat data set
#'
#' As someone who is new to proteomics, I know basically nothing about the MSNset or MSnSet or
#' whatever it is.  I guess this is my chance to learn...
#'
#' @param ProteinName  A protein name!
#' @param PeptideSequence  And the sequence of a peptide!
#' @param PrecursorCharge  The charge of the observed precursor.
#' @param FragmentIon  What is the state of the observed fragment?
#' @param ProductCharge  The charge of the observed product.
#' @param IsotopeLabelType  How were they labeled?
#' @param Bioreplicate  What is the replicate state?
#' @param Run  Which run is this?
#' @param Condition  What is the condition of this sample?
#' @param data  The data set to convert.
#'
#' It seems to me that data should be the first argument and then fill in the
#'   rest of the arguments with default values.
#'
#' @import Rcpp
#' @importFrom MSnbase pData fData MSnSet
#' @export
transformMSnSetToMSstats <- function(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                                     ProductCharge, IsotopeLabelType, Bioreplicate,
                                     Run, Condition, data) {

  if (!inherits(data, "MSnSet")) {
    stop("Only MSnSet class can be converted to input format for MSstats.")
  }
  ## if the data is an expression set, extract the different data components and reshape to "long" format

  ## extract the three data components
  ## Hmmm will Biobase pass this correctly to MSnbase::pData?  I dunno.  I think
  ## so.
  ## If the results look weird, look HERE!
  sampleData <- Biobase::pData(data)
  featureData <- Biobase::fData(data)
  expressionMatrix <- as.data.frame(Biobase::exprs(data))
  ## extract the variable names
  sample.variables <- dimnames(sampleData)[[2]]
  feature.variables <- dimnames(featureData)[[2]] ## ProteinAccession, PeptideSequence, charge
  ## create a patient variable by which to merge with intensities
  sampleData[["pData_rownames"]] <- rownames(sampleData)
  ## create a feature variable by which to merge with intensities
  featureData[["fData_rownames"]] <- rownames(featureData)
  ## transform the matrix of intensities so that each row is a sample and feature combination
  long.abundances <- stats::reshape(
                              expressionMatrix, idvar="fData_rownames",
                              ids=rownames(expressionMatrix), times=colnames(expressionMatrix),
                              timevar="pData_rownames", varying=list(dimnames(sampleData)[[1]]),
                              v.names="ABUNDANCE", direction="long")
  rownames(long.abundances) <- NULL

  ## merge with featureData and patientData to get the feature -> protein and the bio.rep -> group mappings
  long.abundances[["ABUNDANCE"]] <- as.numeric(as.character(long.abundances[["ABUNDANCE"]]))
  long.abundances.2 <- merge(long.abundances, featureData, by="fData_rownames")
  final.data <- merge(long.abundances.2, sampleData, by="pData_rownames")
  ## extract required information
  ## default : Protein="ProteinAccession",PeptideSequence="PeptideSequence",
  ## PrecursorCharge="charge", FragmentIon, ProductCharge, IsotopeLabelType="mz",
  ## Bioreplicate=NA, Run=NA,
  if (missing(ProteinName)) {
    ProteinName <- "ProteinAccession"
  }

  if (missing(PeptideSequence)) {
    PeptideSequence <- "PeptideSequence"
  }
  if (missing(PrecursorCharge)) {
    PrecursorCharge <- "charge"
  }
  if (missing(FragmentIon)) {
    final.data[["FragmentIon"]] <- NA
    FragmentIon <- "FragmentIon"
  }
  if (missing(ProductCharge)) {
    final.data[["ProductCharge"]] <- NA
    ProductCharge <- "ProductCharge"
  }
  if (missing(IsotopeLabelType)) {
    IsotopeLabelType <- "mz"
  }
  if (missing(Bioreplicate)) {
    Bioreplicate <- "pData_rownames"
  }
  if (missing(Run)) {
    Run <- "fileIdx"
  }
  if (missing(Condition)) {
    stop("The condition arguments must be specified.")
  }

  ## if there are any missing variable name, warn it and stop
  check.name <- c(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
                  IsotopeLabelType, Bioreplicate, Run, Condition)
  diff.name <- setdiff(check.name, colnames(final.data))
  if (length(diff.name) > 0) {
    stop(paste("Please check the variable name. The provided variable name",
               paste(diff.name, collapse=","), "is not present in the data set.", sep=" "))
  }
  ## define the group column
  if (length(Condition) > 1) {
    group.variable <- final.data[, which(colnames(final.data) == Condition[1])]
    for (i in 2:length(Condition)) {
      var <- final.data[, which(colnames(final.data) == Condition[i])]
      group.variable <- paste(group.variable, var, sep = ".")
    }
    final.data[["group.column"]] <- group.variable
  } else {
    final.data[["group.column"]] <- final.data[, which(colnames(final.data) == Condition)]
  }

  ## combine into a data frame
  final.data.2 <- final.data[, c(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                                 ProductCharge, IsotopeLabelType,
                                 "group.column", Bioreplicate, Run, "ABUNDANCE")]
  colnames(final.data.2) <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                              "ProductCharge", "IsotopeLabelType", "Condition", "BioReplicate",
                              "Run", "Intensity")
  rownames(final.data.2) <- NULL
  return(final.data.2)
}

#' Generate a MSnSet from a MSstats data set.
#'
#' @param data  Input data from MSstats.
#' @export
transformMSstatsToMSnSet <- function(data) {
  data <- droplevels(data)
  colnames(data) <- toupper(colnames(data))
  ## the expression matrix
  xx <- with(data, by(INTENSITY, list(ISOTOPELABELTYPE, BIOREPLICATE, CONDITION, RUN), c))
  xx <- do.call(cbind, xx)
  ## sample annotation
  pd <- data[!duplicated(data[, c("ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN")]),
             c("ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN")]
  ## feature annotation
  fd <- data[!duplicated(data[, c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE",
                                  "FRAGMENTION", "PRODUCTCHARGE")]),
             c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", "FRAGMENTION", "PRODUCTCHARGE")]
  ## need to make as MSnSet class
  e <- MSnbase::MSnSet(xx, fd, pd)
  MSnbase::sampleNames(e) <- paste(
             e[["ISOTOPELABELTYPE"]],
             e[["CONDITION"]],
             e[["BIOREPLICATE"]],
             e[["RUN"]],
             sep = ".")
  MSnbase::featureNames(e) <- paste(
             MSnbase::fData(e)[["PEPTIDESEQUENCE"]],
             MSnbase::fData(e)[["PRECURSORCHARGE"]],
             MSnbase::fData(e)[["FRAGMENTION"]],
             MSnbase::fData(e)[["PRODUCTCHARGE"]],
             sep = ".")
  if (methods::validObject(e)){
    return(e)
  } else {
    stop("Did not receive a valid object.")
  }
}
