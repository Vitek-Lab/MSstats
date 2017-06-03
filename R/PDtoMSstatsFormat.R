
## Converter for Proteome discoverer output

## output from Proteome discoverer : PSM sheet

#' @export
PDtoMSstatsFormat <- function(input,
                              annotation,
                              useNumProteinsColumn=FALSE,
                              useUniquePeptide=TRUE,
                              summaryforMultipleRows=max,
                              fewMeasurements="remove",
                              removeOxidationMpeptides=FALSE,
                              removeProtein_with1Peptide=FALSE,
                              which.quantification = 'Precursor.Area'){

    ################################################
    ## 0. which intensity : Precursor.Area vs. Intensity
    ################################################
    ## 2017.01.11 : use 'Precursor.Area' instead of 'Intensity'
    ## default : Precursor.Area
    if(which.quantification == 'Intensity'){
        which.quant <- 'Intensity'
    } else {
        which.quant <- 'Precursor.Area'
    }
    
    ################################################
    ## 1. get subset of columns
    ################################################
  
    input <- input[, which(colnames(input) %in% c("Protein.Group.Accessions", "X..Proteins",
                                                "Sequence", "Modifications", "Charge",
                                                "Spectrum.File", which.quant))]
    
    colnames(input)[colnames(input) == 'Protein.Group.Accessions'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'X..Proteins'] <- 'numProtein'
    colnames(input)[colnames(input) == 'Sequence'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'Spectrum.File'] <- 'Run'
    colnames(input)[colnames(input) == 'Precursor.Area'] <- 'Intensity'
  
  
    ################################################
    ## 2. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
  
    if( useNumProteinsColumn ){
        
        ## remove rows with #proteins is not 1
        input <- input[input$numProtein == '1', ]
        
        message('** Rows with #Proteins, which are not equal to 1, are removed.')
        
    }
    
    if( useUniquePeptide ){
    

        ## double check
        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")]) 
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)
        
        ## count how many proteins are assigned for each peptide
        structure <- aggregate(ProteinName ~., data=pepcount, length)
        remove_peptide <- structure[structure$ProteinName != 1, ]
    
        ## remove the peptides which are used in more than one protein
        if( length(remove_peptide$Proteins != 1) != 0 ){
            input <- input[-which(input$Sequence %in% remove_peptide$Sequence), ]
            
            message('** Peptides, that are used in more than one proteins, are removed.')
            
        }
    }
  

    ################################################
    ### 3. remove the peptides including oxidation (M) sequence
    ################################################
  
    if (removeOxidationMpeptides) {
        remove_m_sequence <- unique(input[grep("Oxidation", input$Modifications), "Modifications"])
    
        if(length(remove_m_sequence) > 0){
            input <- input[-which(input$Modifications %in% remove_m_sequence), ]
        }
    
        message('Peptides including oxidation(M) in the Modifications are removed.')
    
    }
 

    ##############################
    ### 4. remove multiple measurements per feature and run
    ##############################
    ## maximum or sum up abundances among intensities for identical features within one run
    input_sub <- dcast( ProteinName + PeptideSequence + Modifications + Charge ~ Run, data=input, 
                        value.var='Intensity', 
                        fun.aggregate=summaryforMultipleRows, 
                        fill=NA_real_) 
 
    ## reformat for long format
    input_sub <- melt(input_sub, id=c('ProteinName', 'PeptideSequence', 'Modifications', 'Charge'))
    colnames(input_sub)[which(colnames(input_sub) %in% c('variable','value'))] <- c("Run", "Intensity")
  
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  
    input <- input_sub

    ##############################
    ### 5. add annotation
    ##############################
    
    noruninfo <- setdiff(unique(input$Run), unique(annotation$Run))
    if ( length(noruninfo) > 0 ) {
        stop( paste('** Annotation for Run :', 
                    paste(noruninfo, collapse = ', '), 
                    ' are needed. Please update them in annotation file.') )
    }
    
    input <- merge(input, annotation, by="Run", all=TRUE)
  
    ## add other required information
    input$FragmentIon <- NA
    input$ProductCharge <- NA
    input$IsotopeLabelType <- "L"
  
    input$PeptideModifiedSequence <- paste(input$PeptideSequence, input$Modifications, sep="_")
  
    input.final <- data.frame(ProteinName = input$ProteinName,
                        PeptideModifiedSequence = input$PeptideModifiedSequence,
                        PrecursorCharge = input$Charge,
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
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        xtmp$feature <- paste(xtmp$PeptideModifiedSequence, xtmp$PrecursorCharge, sep="_")
        count_measure <- xtabs( ~feature, xtmp)
    
        remove_feature_name <- count_measure[count_measure < 3]
    
        input$feature <- paste(input$PeptideModifiedSequence, input$PrecursorCharge, sep="_")
    
        if( length(remove_feature_name) > 0 ){
            input <- input[-which(input$feature %in% names(remove_feature_name)), ]
        }
    
        input <- input[, -which(colnames(input) %in% c('feature'))]
    
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


