
## Converter for MSFragger+Philosopher

## output from MSFragger+Philosopher : MSstats report

#' @export
PhilosophertoMSstatsFormat <- function(input,
                              annotation,
                              which.sequence = 'Peptide.Sequence',
                              useIsUniqueColumn = FALSE,
                              useUniquePeptide=TRUE,
                              summaryforMultipleRows=max,
                              fewMeasurements="remove",
                              removeOxidationMpeptides=FALSE,
                              removeProtein_with1Peptide=FALSE){

    ## check annotation
    required.annotation <- c('Condition', 'BioReplicate', 'Run')
    
    if (!all(required.annotation %in% colnames(annotation))) {
        
        missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))
        
        stop(paste("**", toString(required.annotation[missedAnnotation]), 
                   "is not provided in Annotation. Please check the annotation file.",
                   "'Run' will be matched with 'Spectrum.File' "))
    }
    
    ## check annotation information
    ## get annotation
    annotinfo <- unique(annotation[, c("Run", "Condition", 'BioReplicate')])	

    ## Each Run should has unique information about condition and bioreplicate
    check.annot <- xtabs(~Run, annotinfo)
    if ( any(check.annot > 1) ) {
        stop('** Please check annotation. Each MS run can\'t have multiple conditions or BioReplicates.' )
    }
    
    ################################################
    ## 0.1. which intensity : Intensity column, no any other quantified information
    ################################################

    ################################################
    ## 0.2. which protein id : Protein.Accessions
    ################################################

    ################################################
    ## 0.3. which sequence : Peptide.Sequence or combination with Modified.Peptide.Sequence
    ################################################
    ## default : Peptide.Sequence
    which.seq <- NULL
    
    if (which.sequence == 'Peptide.Sequence') {
        which.seq <- 'Peptide.Sequence'
    } else if (which.sequence == 'Modified.Peptide.Sequence') {
        which.seq <- 'Modified.Peptide.Sequence'
        
        ## Modified.peptide Sequence column has empty for non-modified.peptide.sequence.
        ## Therefore, we need to fill in with peptide.sequence
        input$Modified.Peptide.Sequence <- as.character(input$Modified.Peptide.Sequence)
        input$Peptide.Sequence <- as.character(input$Peptide.Sequence)
        
        input[is.na(input$Modified.Peptide.Sequence), 'Modified.Peptide.Sequence'] <- input[is.na(input$Modified.Peptide.Sequence), 'Peptide.Sequence']
    } 
    
    if (is.null(which.sequence)) {
        stop('** Please select which columns should be used for peptide sequence, between twp options (Peptide.Sequence or Modified.Peptide.Sequence).')
    }
    
    if (!is.element(which.seq, colnames(input))) {
        stop('** Please select which columns should be used for peptide sequence, between twp options (Peptide.Sequence or Modified.Peptide.Sequence).')
    }
    
    ################################################
    ## 1. get subset of columns
    ################################################
  
    input <- input[, which(colnames(input) %in% c("Protein.Accessions", 
                                                which.seq, "Charge",
                                                "Spectrum.File", "Intensity", "Is.Unique"))]
    
    colnames(input)[colnames(input) == "Protein.Accessions"] <- 'ProteinName'
    
    colnames(input)[colnames(input) == which.seq] <- 'PeptideSequence'

    colnames(input)[colnames(input) == 'Spectrum.File'] <- 'Run'
    
    colnames(input)[colnames(input) == "Intensity"] <- 'Intensity'

    ################################################
    ## 2. remove the peptides including oxidation (M) sequence
    ################################################
    
    ## need to check from felipe
    if (removeOxidationMpeptides) {
        remove_m_sequence <- unique(input[grep("Oxidation", input$Modifications), "Modifications"])
    
        if (length(remove_m_sequence) > 0) {
            input <- input[-which(input$Modifications %in% remove_m_sequence), ]
        }
    
        message('Peptides including oxidation(M) in the Modifications are removed.')
    }
 
    ################################################
    ## 3. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    
    if (useIsUniqueColumn) {
        
        ## remove rows with Is.Unique is false
        input <- input[input$Is.Unique != 'false', ]
        
        message('** Rows with Is.Unique = false are removed.')
    }
    
    if (useUniquePeptide) {
        ## If useIsUniqueColumn = T and Is.Unique column is correct, we dont need this part.
        ## But, double check
        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")]) 
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)
        
        ## count how many proteins are assigned for each peptide
        structure <- aggregate(ProteinName ~ ., data=pepcount, length)
        remove_peptide <- structure[structure$ProteinName != 1, ]
        
        ## remove the peptides which are used in more than one protein
        if (length(remove_peptide$Proteins != 1) != 0) {
            input <- input[-which(input$Sequence %in% remove_peptide$Sequence), ]
            
            message('** Peptides, that are used in more than one proteins, are removed.')
        }
    }
    
    ##############################
    ## 4. remove featuares with all na or zero
    ## some rows have all zero values across all MS runs. They should be removed.
    ##############################
    
    input$fea <- paste(input$PeptideSequence,
                       input$Charge,
                       sep="_")
    
    inputtmp <- input[!is.na(input$Intensity) & input$Intensity > 1, ]
    
    count <- inputtmp %>% group_by(fea) %>% summarise(length=length(Intensity))
    
    ## get feature with all NA or zeros
    getfea <- count[count$length > 0, 'fea']
    
    if (nrow(getfea) > 0) {
        
        nfea.remove <- length(unique(input$fea)) - nrow(getfea)
        input <- input[which(input$fea %in% getfea$fea), ]
        
        message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
        
    }
    
    rm(inputtmp)
    input <- input[, -which(colnames(input) %in% c('fea'))]
    

    ##############################
    ## 5. remove multiple measurements per feature and run
    ##############################
    ## maximum or sum up abundances among intensities for identical features within one run
    input_sub <- dcast( ProteinName + PeptideSequence + Charge ~ Run, data=input, 
                        value.var='Intensity', 
                        fun.aggregate=summaryforMultipleRows, 
                        fill=NA_real_) 
    ## need to check whether Modification is needed or not.
    # input_sub <- dcast( ProteinName + PeptideSequence + Modifications + Charge ~ Run, data=input, 
    #                  value.var='Intensity', 
    #                  fun.aggregate=summaryforMultipleRows, 
    #                     fill=NA_real_) 
 
    ## reformat for long format
    #input_sub <- melt(input_sub, id=c('ProteinName', 'PeptideSequence', 'Modifications', 'Charge'))
    input_sub <- melt(input_sub, id=c('ProteinName', 'PeptideSequence', 'Charge'))
    
    colnames(input_sub)[which(colnames(input_sub) %in% c('variable','value'))] <- c("Run", "Intensity")
  
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  
    input <- input_sub

    
    ##############################
    ##  6. remove features which has 1 or 2 measurements across runs
    ##############################
    if (fewMeasurements=="remove") {
        
        ## it is the same across experiments. # measurement per feature. 
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        xtmp$feature <- paste(xtmp$PeptideSequence, xtmp$Charge, sep="_")
        count_measure <- xtabs( ~feature, xtmp)
    
        remove_feature_name <- count_measure[count_measure < 3]
    
        input$feature <- paste(input$PeptideSequence, input$PrecursorCharge, sep="_")
    
        if (length(remove_feature_name) > 0) {
            input <- input[-which(input$feature %in% names(remove_feature_name)), ]
        }
    
        input <- input[, -which(colnames(input) %in% c('feature'))]
    }
  
    ##############################
    ##  8. remove proteins with only one peptide and charge per protein
    ##############################
	
	if (removeProtein_with1Peptide) {
	  
	    ## remove protein which has only one peptide
	    input$feature <- paste(input$PeptideSequence, 
	                           input$Charge, 
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
    
    ##############################
    ## 9. add annotation
    ##############################
    
    noruninfo <- setdiff(unique(input$Run), unique(annotation$Run))
    if (length(noruninfo) > 0) {
        stop(paste('** Annotation for Run :', 
                   paste(noruninfo, collapse = ', '), 
                   ' are needed. Please update them in annotation file.') )
    }
    
    input <- merge(input, annotation, by="Run", all=TRUE)
    
    ## add other required information
    # need to update
    #input$PeptideModifiedSequence <- paste(input$PeptideSequence, input$Modifications, sep="_")
    
    input.final <- data.frame("ProteinName" = input$ProteinName,
                              "PeptideSequence" = input$PeptideSequence,
                              "PrecursorCharge" = input$Charge,
                              "FragmentIon" = NA,
                              "ProductCharge" = NA,
                              "IsotopeLabelType" = "L",
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
  
	return(input)
}


