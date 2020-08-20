
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
    ## 0. check input
    ##############################

    ## there are space in column name

    if (!is.element('Modifications', input[2, ])) {
        input[2, ][grep('modifications', input[2 ,])] <- 'Modifications'
    }

    if (!is.element('Protein', input[2,]) & is.element('Accession', input[2,])) {
        ## use 'Accession' for Protein ID
        input[2, ][input[2, ] == 'Accession'] <- 'Protein'
    }
    
    required.column <- c('Protein', 'Sequence', 'Charge', 'Modifications')
    
    if (length(grep('quantitation', input[2, ])) > 0){
        input[2, ][grep('quantitation', input[2, ])] <- 'Use.in.quantitation'
        required.column <- c(required.column, 'Use.in.quantitation')
    }
    
    if (!all(required.column %in% input[2, ])) {
        
        missedInput <- which(!(required.column %in% input[2, ]))

        stop(paste("**", toString(required.column[missedInput]), 
                   "is not provided in input. Please check the input."))
    }
    
    ## check annotation
    required.annotation <- c('Run', 'BioReplicate', 'Condition')
    
    if (!all(required.annotation %in% colnames(annotation))) {
        
        missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))
        
        stop(paste("**", toString(required.annotation[missedAnnotation]), 
                   "is not provided in Annotation. Please check the annotation."))
    }
    
    ## check annotation information
    ## get annotation
    annotinfo <- unique(annotation[, c("Run", "Condition", 'BioReplicate')])	
    
    ## Each Run should has unique information about condition and bioreplicate
    check.annot <- xtabs(~Run, annotinfo)
    if ( any(check.annot > 1) ) {
        stop('** Please check annotation. Each MS run can\'t have multiple conditions or BioReplicates.' )
    }

    ## get abundance information
    if (is.element('Raw.abundance', colnames(input)) & 
        is.element('Normalized.abundance', colnames(input))) {
        
        start.column <- which(colnames(input) == 'Raw.abundance')
        check.numRun <- which(colnames(input) == 'Normalized.abundance')
        
        if (start.column-check.numRun != nrow(annotation)) {
            stop(message('** Please check annotation file. The numbers of MS runs in annotation and output are not matched.'))
        }
        
        raw.abundance.column <- c(start.column:(start.column + nrow(annotation)-1))
        
        input <- input[, c(which(input[2, ] %in% required.column),
                           raw.abundance.column)]
         
    } else if (is.element('Raw.abundance', colnames(input))) {
        
        start.column <- which(colnames(input) == 'Raw.abundance')
        raw.abundance.column <- c(start.column:(start.column + nrow(annotation)-1))
        
        input <- input[, c(which(input[2, ] %in% required.column),
                           raw.abundance.column)]
    }
    
    input <- input[-1, ]
    colnames(input) <- input[1, ]
    input <- input[-1, ]
  
    ##############################
    ## 1. use only 'use in quantitation = true'
    ##############################

    if (is.element('Use.in.quantitation', colnames(input))) {
        ## value for use in quantitation is True vs False
        if (length( grep('True', unique(input$Use.in.quantitation))) > 0) { 
            input <- input[input$Use.in.quantitation == 'True', ]
        } else if (length(grep('TRUE', unique(input$Use.in.quantitation))) > 0) {
            input <- input[input$Use.in.quantitation == TRUE, ]
        }
        
        input <- input[, -which(colnames(input) %in% c('Use.in.quantitation'))]
    }
    
    ##############################
    ## 2. modify column names and remove some columnts
    ##############################
    
    input <- input[!is.na(input$Protein) & input$Protein != '', ]
    input <- input[!is.na(input$Sequence) & input$Sequence != '', ]
    
    ## get modified sequence
    input$ModifiedSequence <- paste(input$Sequence, input$Modifications, sep="")
  
    ## remove completely duplicated rows
    input <- input[!duplicated(input), ]
  
    ################################################
    ## 3. remove the peptides including oxidation (M) sequence
    if (removeOxidationMpeptides) {
        remove_m_sequence <- unique(input[grep("Oxidation", input$ModifiedSequence), "ModifiedSequence"])
    
        if (length(remove_m_sequence) > 0) {
            input <- input[-which(input$ModifiedSequence %in% remove_m_sequence), ]
        }
    
        message('Peptides including oxidation(M) in the sequence are removed.')
    }
 
    ################################################
    ## 4. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    if (useUniquePeptide) {
    
        pepcount <- unique(input[, c("Protein", "Sequence")]) 
        pepcount$Sequence <- factor(pepcount$Sequence)
    
        ## count how many proteins are assigned for each peptide
        structure <- aggregate(Protein ~ ., data=pepcount, length)
        remove_peptide <- structure[structure$Protein!=1, ]
    
        ## remove the peptides which are used in more than one protein
        if (length(remove_peptide$Protein != 1) != 0) {
            input <- input[-which(input$Sequence %in% remove_peptide$Sequence), ]
        }
    
        message('** Peptides, that are used in more than one proteins, are removed.')
    }
  
    ##############################
    ## 5. remove multiple measurements per feature and run
    ##############################
    input <- input[, -which(colnames(input) %in% c('Sequence', 'Modifications'))]
    input_remove <- melt(input, id=c('Protein', 'ModifiedSequence', 'Charge'))

    colnames(input_remove) <- c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "Run", "Intensity")
    input_remove$Intensity <- as.double(input_remove$Intensity)
  
    ## maximum or sum up abundances among intensities for identical features within one run
    input <- dcast(ProteinName + PeptideModifiedSequence + PrecursorCharge ~ Run, data=input_remove, 
                   value.var='Intensity', 
                   fun.aggregate=summaryforMultipleRows, na.rm=T,
                   fill=NA_real_) 
  
    ## reformat for long format
    input <- melt(input, id=c('ProteinName', 'PeptideModifiedSequence', 'PrecursorCharge'))
    colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
  
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  
    ##############################
    ## 6. add annotation
    ##############################
    input <- merge(input, annotation, by="Run", all=TRUE)
  
    ## add other required information
    input$FragmentIon <- NA
    input$ProductCharge <- NA
    input$IsotopeLabelType <- "L"
    
    input.final <- data.frame("ProteinName" = input$ProteinName,
                              "PeptideModifiedSequence" = input$PeptideModifiedSequence,
                              "PrecursorCharge" = input$PrecursorCharge,
                              "FragmentIon" = input$FragmentIon,
                              "ProductCharge" = input$ProductCharge,
                              "IsotopeLabelType" = input$IsotopeLabelType,
                              "Condition" = input$Condition,
                              "BioReplicate" = input$BioReplicate,
                              "Run" = input$Run,
                              "Intensity" = input$Intensity)
    
    if (any(is.element(colnames(input), 'Fraction'))) {
        input.final <- data.frame(input.final,
                                  "Fraction" = input$Fraction)
    }
    
    input <- input.final
    rm(input.final)
  
    ##############################
    ##  7. remove features which has 1 or 2 measurements across runs
    ##############################
    if (fewMeasurements == "remove") {
    
        ## it is the same across experiments. # measurement per feature. 
        input <- .remove_feature_with_few_progenesis(input)
    }
  

    ##############################
    ##  8. remove proteins with only one peptide and charge per protein
    ##############################
	
	if (removeProtein_with1Peptide) {
	    ##remove protein which has only one peptide
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
	        message(paste0("** ", length(removepro), 
	                       ' proteins, which have only one feature in a protein, are removed among ', 
	                       lengthtotalprotein, ' proteins.'))
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
  
    if (length(remove_feature_name) > 0) {
        x <- x[-which(x$feature %in% names(remove_feature_name)), ]
    }

    x <- x[, -which(colnames(x) %in% c('feature'))]
  
    return(x)
}


