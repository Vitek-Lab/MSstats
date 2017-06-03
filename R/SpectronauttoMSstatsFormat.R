
## Converter for Spectronaut output

## output from Spectronaut : long format 
## columns : Condition, BioReplicate, Run, ProteinName, FragmentIon, PeptideSequence            
##           ProductCharge, PrecursorCharge, IsotopeLabelType, Intensity, 
##           F.ExcludedFromQuantification

#' @export
SpectronauttoMSstatsFormat <- function(input, 
                                       intensity = 'PeakArea',
                                       filter_with_Qvalue = TRUE,
                                       qvalue_cutoff = 0.01,
                                       useUniquePeptide = TRUE,
                                       fewMeasurements="remove",
                                       removeProtein_with1Feature = FALSE,
                                       summaryforMultipleRows=max){

    ## Check correct option or input
    requiredinput.general <- c('F.FrgLossType', 'F.ExcludedFromQuantification',
                               'PG.ProteinGroups', 'EG.ModifiedSequence', 'FG.Charge',
                               'F.FrgIon', 'R.Condition', 
                               'R.FileName', 'R.Replicate', 'EG.Qvalue')
    
    requiredinput.int <- c('F.PeakArea', 'F.NormalizedPeakArea')
    
    requiredinput.charge <- c('F.Charge', 'F.FrgZ')
    
    ################################
    ## general input
    ################################
    
    if ( !all( requiredinput.general %in% colnames(input) ) ) {
        
        misssing.col <- requiredinput.general[!requiredinput.general %in% colnames(input)]
        
        stop(paste("** Please check the required input. The required input needs '", 
                   paste(missing.col, collapse = ", "), "'", sep=""))
    }
    
    ## intensity columns
    if ( sum( requiredinput.int %in% colnames(input) ) == 0 ) {
        
        stop( paste("** Please check the required input. The required input needs at least one of '", 
                   paste(requiredinput.int, collapse = "' or '"), "'", sep="") )
    }
    
    ## general input
    if ( sum( requiredinput.charge %in% colnames(input) ) == 0 ) {
        
        stop( paste("** Please check the required input. The required input needs at least one of '", 
                    paste(requiredinput.charge, collapse = "' or '"), "'", sep="") )
    }

    
    ##############################
    ### 1. loss type : use only 'no loss'
     ##############################
    input <- input[input$F.FrgLossType == 'noloss', ]
  
    ##############################
    ### 2. use only 'F.ExcludedFromQuantification == False' : XIC quality
    ##############################
  
    input <- input[input$F.ExcludedFromQuantification == 'False', ]
  
    ##############################
    ### 3. get useful subset of column
    ##############################
    if(is.element('F.Charge', colnames(input))){
        f.charge <- 'F.Charge'
    } else if(is.element('F.FrgZ', colnames(input))) {
        f.charge <- 'F.FrgZ'
    } else {
        f.charge <- NULL
    }
    
    subsetcolumn <- c('PG.ProteinGroups', 'EG.ModifiedSequence', 'FG.Charge',
                    'F.FrgIon', f.charge,
                    'R.Condition', 'R.FileName', 'R.Replicate',
                    'EG.Qvalue')
  
    if ( intensity == 'NormalizedPeakArea' ){
        ## use normalized peak area by SN
        input <- input[, c(subsetcolumn, 'F.NormalizedPeakArea')]
    } else {
        ## use original peak area without any normalization
        input <- input[, c(subsetcolumn, 'F.PeakArea')]
    }
  
    colnames(input)[colnames(input) == 'PG.ProteinGroups'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'EG.ModifiedSequence'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'FG.Charge'] <- 'PrecursorCharge'
    colnames(input)[colnames(input) == 'F.FrgIon'] <- 'FragmentIon'
    colnames(input)[colnames(input) == f.charge] <- 'ProductCharge'
    colnames(input)[colnames(input) == 'R.Condition'] <- 'Condition'
    colnames(input)[colnames(input) == 'R.FileName'] <- 'Run'
    colnames(input)[colnames(input) == 'R.Replicate'] <- 'BioReplicate'
    colnames(input)[colnames(input) == 'F.PeakArea'] <- 'Intensity'
    colnames(input)[colnames(input) == 'F.NormalizedPeakArea'] <- 'Intensity'
    colnames(input)[colnames(input) == 'EG.Qvalue'] <- 'Qvalue'
  
    ##############################
    ### 4. filter by Qvalue
    ##############################

    if( filter_with_Qvalue ){
    
        if( !is.element(c('Qvalue'), colnames(input)) ){
      
            stop('** EG.Qvalue column is needed in order to filter out by Qvalue. Please add EG.Qvalue column in the input.')
      
        } else {
      
            ## when qvalue > qvalue_cutoff, replace with zero for intensity
            input[!is.na(input$Qvalue) & input$Qvalue > qvalue_cutoff, "Intensity"] <- 0
      
            message(paste('** Intensities with great than ', qvalue_cutoff, ' in EG.Qvalue are replaced with zero.', sep=''))
        }
    }
  
    ##############################
    ### 5. remove featuares with all na or zero
    ### some rows have all zero values across all MS runs. They should be removed.
    ##############################
  
    input$fea <- paste(input$PeptideSequence,
                        input$PrecursorCharge,
                        input$FragmentIon,
                        input$ProductCharge,
                        sep="_")
  
    input1 <- input[is.na(input$Intensity), ]
    input2 <- input[!is.na(input$Intensity) & input$Intensity <= 1, ]
    inputtmp <- rbind(input1, input2)
  
    count <- aggregate(Intensity ~ fea, data=inputtmp, FUN=length)
  
    ## get feature with all NA or zeros
    getfea <- count[count$Intensity == length(unique(input$Run)), 'fea']
  
    if( length(getfea) > 0 ){
        input <- input[-which(input$fea %in% getfea), ]
    }
  

    ################################################
    ## 6. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    if( useUniquePeptide ){
    
        pepcount <- unique(input[, c("ProteinName","PeptideSequence")]) ## Protein.group.IDs or Sequence
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)
    
        ## count how many proteins are assigned for each peptide
        structure <- aggregate(ProteinName ~., data=pepcount, length)
        remove_peptide <- structure[structure$ProteinName != 1, ]
    
        ## remove the peptides which are used in more than one protein
        if( length(remove_peptide$ProteinName != 1 ) != 0 ){
            input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]
      
            message('** Peptides, that are used in more than one proteins, are removed.')
        } else {
            message('** All peptides are unique peptides in proteins.')
        }
    
        rm(structure)
        rm(remove_peptide)
    }
  
    ##############################
    ###  7. remove features which has 1 or 2 measurements across runs
    ##############################
    if( fewMeasurements=="remove" ){
    
        ## it is the same across experiments. # measurement per feature. 
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        count_measure <- xtabs( ~fea, xtmp)
    
        remove_feature_name <- count_measure[count_measure < 3]
    
        if( length(remove_feature_name) > 0 ){
            input <- input[-which(input$fea %in% names(remove_feature_name)), ]
        }
    
    }
  
    ##############################
    ### 8. remove proteins with only one peptide and charge per protein
    ##############################
  
    if( removeProtein_with1Feature ){
        ######## remove protein which has only one peptide
        tmp <- unique(input[, c("ProteinName", 'fea')])
        tmp$ProteinName <- factor(tmp$ProteinName)
        count <- xtabs( ~ ProteinName, data=tmp)
        lengthtotalprotein <- length(count)
    
        removepro <- names(count[count <= 1])
    
        if (length(removepro) > 0) {
      
            input <- input[-which(input$ProteinName %in% removepro), ]
            message(paste("*** ", length(removepro), ' proteins, which have only one feature in a protein, are removed among ', lengthtotalprotein, ' proteins.', sep=""))
        }
    }
  

    ##############################
    ### 9. remove multiple measurements per feature and run
    ##############################
  
    count <- aggregate(Intensity ~ fea, data=input, FUN=length)
  
    ## if any feature has more number of total MS runs, 
    if( any(unique(count$Intensity) > length(unique(input$Run))) ){
    
        annotation <- unique(input[, c("Condition","BioReplicate","Run")])
    
        ## maximum or sum up abundances among intensities for identical features within one run
        input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon + ProductCharge ~ Run, data=input, 
                          value.var='Intensity', 
                          fun.aggregate=summaryforMultipleRows, fill=NA_real_) 
    
        ## reformat for long format
        input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon', 'ProductCharge'))
        colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    
        ### add annotation
        input <- merge(input, annotation, by="Run")
    
        ## reorder columns
        input <- input[, c(2,3,4,5,6,8,1,9,7)]
    
        message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
    
    } else {
        ## remove column, named as 'fea'
        input <- input[, -which(colnames(input) %in% c('fea', 'Qvalue'))]
        message('** No multiple measurements in a feature and a run.')
    }
  
    input$IsotopeLabelType <- 'L'
    input$ProteinName <- input$ProteinName
  
	return(input)
}


