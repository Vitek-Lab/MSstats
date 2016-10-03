
## Converter for Spectronaut output

## output from Spectronaut : long format 
## columns : Condition, BioReplicate, Run, ProteinName, FragmentIon, PeptideSequence            
##           ProductCharge, PrecursorCharge, IsotopeLabelType, Intensity, 
##           F.ExcludedFromQuantification

#' @export
SpectronauttoMSstatsFormat <- function(input, 
                                       summaryforMultipleRows=max){
  
  ##############################
  ### 1. use only 'F.ExcludedFromQuantification != True'
  ##############################
  
  input <- input[input$F.ExcludedFromQuantification !="True",]


  ##############################
  ### 2. remove multiple measurements per feature and run
  ##############################
  
  annotation <- unique(input[, c("Condition","BioReplicate","Run", "IsotopeLabelType")])
  
  ## maximum or sum up abundances among intensities for identical features within one run
  input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon + ProductCharge ~ Run, data=input, 
                   value.var='Intensity', 
                   fun.aggregate=summaryforMultipleRows, fill=NA_real_) 
  
  ## reformat for long format
  input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon', 'ProductCharge'))
  colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
  
  message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  
  ##############################
  ### 3. add annotation
  ##############################
  input <- merge(input, annotation, by="Run")
  
  input <- input[, c(2,3,4,5,6,10,8,9,1,7)]
  
  input$ProteinName <- input$ProteinName
  
	return(input)
}


