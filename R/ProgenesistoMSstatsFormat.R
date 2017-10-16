
## Converter for Progenesis output
## Thanks Ulrich Omasits : use R scripts by Ulrich Omasits, 2015, version 2.1

## output from Progenesis : wide format 

#' @export
ProgenesistoMSstatsFormat <- function(input, 
                                      annotation,
                                      useUniquePeptide=TRUE,
                                      summaryforMultipleRows=max,
                                      fewMeasurements="remove",
                                      removeOxidationMpeptides=FALSE,
                                      removeProtein_with1Peptide=FALSE){

    ##############################
    ### 1. use only 'use in quantitation = true'
    ##############################
    #position.Normedpeak <- which(colnames(input) == 'Normalized.abundance')
    #position.Rawpeak <- which(colnames(input) == 'Raw.abundance')
    #position.Spectra <- which(colnames(input) == 'Spectral.counts')
    
    if( is.element('Spectral.counts', colnames(input)) & 
        is.element('Raw.abundance', colnames(input)) & 
        is.element('Normalized.abundance', colnames(input)) ){
        
        input <- input[,
                       c(1:(which(colnames(input) == 'Normalized.abundance')-1),
                         which(colnames(input) == 'Raw.abundance') : (which(colnames(input) == 'Spectral.counts')-1))]
        
    } else if ( is.element('Raw.abundance', colnames(input)) & 
                is.element('Normalized.abundance', colnames(input)) ){
       
         input <- input[,
                       c(1:(which(colnames(input) == 'Normalized.abundance')-1),
                         which(colnames(input) == 'Raw.abundance') : ncol(input))]
         
    } else if ( is.element('Raw.abundance', colnames(input)) ) {
        
        input <- input[,
                       c(1:(which(colnames(input) == 'Raw.abundance')-1),
                         which(colnames(input) == 'Raw.abundance') : ncol(input))]
    }
    
    
    input <- input[-1,]
    colnames(input) <- input[1, ]
    input <- input[-1,]
  
    ##############################
    ### 2. use only 'use in quantitation = true'
    ##############################
    ## there are space in column name
    colnames(input)[grep('quantitation', colnames(input))] <- 'Use.in.quantitation'
  
    ## value for use in quantitation is True vs False
    if( length( grep('True', unique(input$Use.in.quantitation)) ) > 0 ){ 
        input <- input[input$Use.in.quantitation == 'True', ]
    } else if( length( grep('TRUE', unique(input$Use.in.quantitation)) ) > 0){
        input <- input[input$Use.in.quantitation == TRUE, ]
    }

    ##############################
    ### 2. modify column names and remove some columnts
    ##############################
  
    ## get modified sequence
    input$ModifiedSequence <- paste(input$Sequence, input$Modifications, sep="")
  
    ## use 'Accession' for Protein ID
    colnames(input)[colnames(input) == 'Accession'] <- 'Protein'
  
    ## get subset of datasets
    input <- input[, which(colnames(input) %in% c('Protein', 'Sequence', 'Charge', 'ModifiedSequence',
                                                as.character(annotation$Run)))]
    ## remove completely duplicated rows
    input <- input[!duplicated(input), ]
  
    ################################################
    ### 3. remove the peptides including oxidation (M) sequence
    if (removeOxidationMpeptides) {
        remove_m_sequence <- unique(input[grep("Oxidation", input$ModifiedSequence), "ModifiedSequence"])
    
        if(length(remove_m_sequence) > 0){
            input <- input[-which(input$ModifiedSequence %in% remove_m_sequence), ]
        }
    
        message('Peptides including oxidation(M) in the sequence are removed.')
    
    }
 
    ################################################
    ## 4. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    if(useUniquePeptide){
    
        pepcount <- unique(input[, c("Protein","Sequence")]) 
        pepcount$Sequence <- factor(pepcount$Sequence)
    
        ## count how many proteins are assigned for each peptide
        structure <- aggregate(Protein ~., data=pepcount, length)
        remove_peptide <- structure[structure$Protein!=1, ]
    
        ## remove the peptides which are used in more than one protein
        if(length(remove_peptide$Protein != 1) != 0){
            input <- input[-which(input$Sequence %in% remove_peptide$Sequence), ]
        }
    
        message('** Peptides, that are used in more than one proteins, are removed.')
    
    }
  
    ##############################
    ### 4. remove multiple measurements per feature and run
    ##############################
    input <- input[, -which(colnames(input) %in% 'Sequence')]
    input_remove <- melt(input, id=c('Protein', 'ModifiedSequence', 'Charge'))

    colnames(input_remove) <- c("ProteinName","PeptideModifiedSequence","PrecursorCharge","Run","Intensity")
    input_remove$Intensity <- as.double(input_remove$Intensity)
  
    ## maximum or sum up abundances among intensities for identical features within one run
    input <- dcast( ProteinName + PeptideModifiedSequence + PrecursorCharge ~ Run, data=input_remove, 
                   value.var='Intensity', 
                   fun.aggregate=summaryforMultipleRows, 
                   fill=NA_real_) 
  
    ## reformat for long format
    input <- melt(input, id=c('ProteinName', 'PeptideModifiedSequence', 'PrecursorCharge'))
    colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
  
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  
  
    ##############################
    ### 5. add annotation
    ##############################
    input <- merge(input, annotation, by="Run", all=TRUE)
  
    ## add other required information
    input$FragmentIon <- NA
    input$ProductCharge <- NA
    input$IsotopeLabelType <- "L"
    
    input.final <- data.frame(ProteinName = input$ProteinName,
                              PeptideModifiedSequence = input$PeptideModifiedSequence,
                              PrecursorCharge = input$PrecursorCharge,
                              FragmentIon = input$FragmentIon,
                              ProductCharge = input$ProductCharge,
                              IsotopeLabelType = input$IsotopeLabelType,
                              Condition = input$Condition,
                              BioReplicate = input$BioReplicate,
                              Run = input$Run,
                              Intensity = input$Intensity)
    
    if( any(is.element(colnames(input), 'Fraction')) ) {
        input.final <- data.frame(input.final,
                                  Fraction = input$Fraction)
    }
    
    input <- input.final
    rm(input.final)
  
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
	    input$feature <- paste(input$PeptideModifiedSequence,
	                           input$PrecursorCharge,
	                           input$FragmentIon,
	                           input$ProductCharge,
	                           sep="_")
	  
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
    xtmp$feature <- paste(xtmp$PeptideModifiedSequence, xtmp$PrecursorCharge, sep="_")
    count_measure <- xtabs( ~feature, xtmp)
  
    remove_feature_name <- count_measure[count_measure < 3]
  
    x$feature <- paste(x$PeptideModifiedSequence, x$PrecursorCharge, sep="_")
  
    if( length(remove_feature_name) > 0 ){
        x <- x[-which(x$feature %in% names(remove_feature_name)), ]
    }

    x <- x[, -which(colnames(x) %in% c('feature'))]
  
    return(x)
}


