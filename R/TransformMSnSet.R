
#######################################
## MSnSet -> input for MSstats
#######################################

transformMSnSetToMSstats  <-  function(ProteinName,
										PeptideSequence, 
										PrecursorCharge, 
										FragmentIon, 
										ProductCharge, 
										IsotopeLabelType, 
										Bioreplicate,
										Run, 
										Condition, 
										data) {
  
  	if (!inherits(data, "MSnSet")) {
    	stop("Only MSnSet class can be converted to input format for MSstats.")
  	}
  	## if the data is an expression set, extract the different data components and reshape to "long" format 
  
  	## extract the three data components
  	sampleData  <-  pData(data)
  	featureData  <-  fData(data)
  	expressionMatrix  <-  as.data.frame(exprs(data))
  
  	## extract the variable names 
  	sample.variables  <-  dimnames(sampleData)[[2]]
  	feature.variables  <-  dimnames(featureData)[[2]] ## ProteinAccession, PeptideSequence, charge
  
  	## create a patient variable by which to merge with intensities 
  	sampleData$pData_rownames  <-  rownames(sampleData)
  
  	## creat a feature variable by which to merge with intensities
  	featureData$fData_rownames  <-  rownames(featureData)
  
  	## transform the matrix of intensities so that each row is a sample and feature combination
  	long.abundances  <-  reshape(expressionMatrix, idvar = "fData_rownames", ids = rownames(expressionMatrix), times = colnames(expressionMatrix), timevar = "pData_rownames", varying = list(dimnames(sampleData)[[1]]), v.names = "ABUNDANCE", direction = "long")
  	rownames(long.abundances)  <-  NULL
  
  	## merge with featureData and patientData to get the feature -> protein and the bio.rep -> group mappings
  	long.abundances$ABUNDANCE  <-  as.numeric(as.character(long.abundances$ABUNDANCE))
  	long.abundances.2  <-  merge(long.abundances, featureData, by = "fData_rownames")
  	final.data  <-  merge(long.abundances.2, sampleData, by = "pData_rownames")
  
  	## extract required information
  	## default : Protein="ProteinAccession",PeptideSequence="PeptideSequence", PrecursorCharge="charge", FragmentIon, ProductCharge ,IsotopeLabelType="mz", Bioreplicate=NA,Run=NA,
  
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
    	final.data$FragmentIon <- NA
    	FragmentIon <- "FragmentIon"
  	}
  	
  	if (missing(ProductCharge)) {
    	final.data$ProductCharge <- NA
    	ProductCharge <- "ProductCharge"
  	} 
  	
  	if (missing(IsotopeLabelType)) {
  		IsotopeLabelType <- "mz"
  	}
  	
  	if (missing(Bioreplicate)) {
  		Bioreplicate <- "pData_rownames"
  	}
  	
  	if (missing(Run)) {
  		Run <- "file"
  	}
  	
  	if (missing(Condition)) {
  		stop("The condition arguments must be specified.")
  	}
  
  
  	## if there are any missing variable name, warn it and stop
  	check.name <- c(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Bioreplicate, Run, Condition)
  
  	diff.name <- setdiff(check.name, colnames(final.data))
  	if (length(diff.name) > 0){
    	stop(paste("Please check the variable name. The provided variable name", paste(diff.name, collapse=","), "is not present in the data set.", sep=" "))
  	}
  
  	## define the group column 
  
  	if (length(Condition) > 1) {
    	group.variable  <-  final.data[, which(colnames(final.data) == Condition[1])]
    
    	for(i in 2:length(Condition)) {
      		var  <-  final.data[, which(colnames(final.data) == Condition[i])]
      		group.variable  <-  paste(group.variable, var, sep = ".")
    	}
    	final.data$group.column  <-  group.variable
  	} else {
  		final.data$group.column  <-  final.data[, which(colnames(final.data) == Condition)]
  	}
  
  	## combine into a data frame
  	final.data.2  <-  final.data[, c(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, "group.column", Bioreplicate, Run, "ABUNDANCE")]
  
  	colnames(final.data.2) <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition", "BioReplicate", "Run", "Intensity")
  	rownames(final.data.2) <- NULL
  
  	return(final.data.2)
}


#######################################
#### MSnSet -> input for MSstats
#######################################

transformMSstatsToMSnSet <- function(data) {
  	
  	data <- droplevels(data)
  
  	colnames(data) <- toupper(colnames(data))
  
  	## the expression matrix
  	xx <- with(data, by(INTENSITY, list(ISOTOPELABELTYPE, BIOREPLICATE, CONDITION,RUN), c))
  	xx <- do.call(cbind, xx)
  
  	## sample annotation
  	pd <- data[!duplicated(data[, c("ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN")]), c("ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN")]
  
  	## feature annotation
  	fd <- data[!duplicated(data[, c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", "FRAGMENTION", "PRODUCTCHARGE")]), c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", "FRAGMENTION", "PRODUCTCHARGE")]
  
  	## need to make as MSnSet class
  	e  <-  MSnSet(xx, fd, pd)
  
  	sampleNames(e) <- paste(e$ISOTOPELABELTYPE, e$CONDITION, e$BIOREPLICATE, e$RUN, sep = ".")
  	featureNames(e) <- paste(fData(e)$PEPTIDESEQUENCE, fData(e)$PRECURSORCHARGE, fData(e)$FRAGMENTION, fData(e)$PRODUCTCHARGE, sep = ".")
  
  	if (validObject(e)){
  		return(e)
  	}
    
}


