
## Converter for Progenesis output
## Thanks Ulrich Omasits : use R scripts by Ulrich Omasits, 2015, version 2.1

## output from Progenesis : wide format 

#' @export
ProgenesistoMSstatsFormat <- function(input, 
                                      annotation,
                                      useUniquePeptide=TRUE,
                                      summaryforMultipleRows=max,
                                      fewMeasurements="remove",
                                      removeProtein_with1Peptide=FALSE){
  
  ##############################
  ### 1. use only 'use in quantitation = true'
  ##############################
  
  input <- input[input$Use.in.quantitation == TRUE,]

  
  ##############################
  ### 2. modify column names and remove some columnts
  ##############################

  input$Sequence <- paste(input$Sequence, input$Modifications, sep="")
  colnames(input)[colnames(input) == 'Accession'] <- 'Protein'
  
  input <- input[, which(colnames(input) %in% c('Protein', 'Sequence', 'Charge', as.character(annotation$Run)))]
  
  ################################################
  ## 3. remove peptides which are used in more than one protein
  ## we assume to use unique peptide
  ################################################
  if(useUniquePeptide){
    
    pepcount <- unique(input[, c("Protein","Sequence")]) 
    pepcount$Sequence <- factor(pepcount$Sequence)
    
    ## count how many proteins are assigned for each peptide
    structure <- aggregate(Protein ~., data=pepcount, length)
    remove_peptide <- structure[structure$Proteins!=1, ]
    
    ## remove the peptides which are used in more than one protein
    if(length(remove_peptide$Proteins != 1) != 0){
      input <- input[-which(input$Sequence %in% remove_peptide$Sequence), ]
    }
    
    message('** Peptides, that are used in more than one proteins, are removed.')
    
  }
  
  
  ##############################
  ### 4. remove multiple measurements per feature and run
  ##############################
  input_remove <- melt(input, id=c('Protein', 'Sequence', 'Charge'))

  colnames(input_remove) <- c("ProteinName","PeptideSequence","PrecursorCharge","Run","Intensity")

  ## maximum or sum up abundances among intensities for identical features within one run
  input <- dcast( ProteinName + PeptideSequence + PrecursorCharge ~ Run, data=input_remove, 
                   value.var='Intensity', 
                   fun.aggregate=summaryforMultipleRows, fill=NA_real_) 
  
  ## reformat for long format
  input <- melt(input, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge'))
  colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
  
  ##############################
  ### 5. add annotation
  ##############################
  input <- merge(input, annotation, by="Run")
  
  ## add other required information
  input$FragmentIon <- NA
  input$ProductCharge <- NA
  input$IsotopeLabelType <- "L"
  
  input <- input[, c(2,3,4,8,9,10,6,7,1,5)]
  
  message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  
  
  ##############################
  ###  6. remove features which has 1 or 2 measurements across runs
  ##############################
  if(fewMeasurements=="remove"){
    
    ## it is the same across experiments. # measurement per feature. 
    input <- .remove_feature_with_few_progenesis(input)
  }
  

  ##############################
  ###  7. remove proteins with only one peptide and charge per protein
  ##############################
	
	if(removeProtein_with1Peptide){
	  ######## remove protein which has only one peptide
	  input$feature <- paste(input$PeptideSequence, input$PrecursorCharge, input$FragmentIon, input$ProductCharge, sep="_")
	  
	  tmp <- unique(input[, c("ProteinName", 'feature')])
	  tmp$ProteinName <- factor(tmp$ProteinName)
	  count <- xtabs( ~ ProteinName, data=tmp)
    lengthtotalprotein <- length(count)
    
	  removepro <- names(count[count <= 1])
	  
	  if (length(removepro) > 0) {
	    
	    input <- input[-which(input$ProteinName %in% removepro), ]
	    message(paste("** ", length(removepro), ' proteins, which have only one feature in a protein, are removed among ', lengthtotalprotein, ' proteins.', sep=""))
	  }
	  
	  input <- input[, -which(colnames(input) %in% c('feature'))]
	}
  
  input$ProteinName <- input$ProteinName
  
	return(input)
}


.remove_feature_with_few_progenesis <- function(x){
  
  xtmp <- x[!is.na(x$Intensity) & x$Intensity > 0, ]
  xtmp$feature <- paste(xtmp$PeptideSequence, xtmp$PrecursorCharge, sep="_")
  count_measure <- xtabs( ~feature, xtmp)
  
  remove_feature_name <- count_measure[count_measure < 3]
  
  x$feature <- paste(x$PeptideSequence, x$PrecursorCharge, sep="_")
  
  if( length(x) > 0 ){
    x <- x[-which(x$feature %in% names(remove_feature_name)), ]
  }

  x <- x[, -which(colnames(x) %in% c('feature'))]
  
  return(x)
}


