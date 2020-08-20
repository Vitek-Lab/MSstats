
## Converter for Spectronaut output

## output from Spectronaut : long format 
## columns : Condition, BioReplicate, Run, ProteinName, FragmentIon, PeptideSequence            
##           ProductCharge, PrecursorCharge, IsotopeLabelType, Intensity, 
##           F.ExcludedFromQuantification

#' @export
SpectronauttoMSstatsFormat <- function(input, 
                                       annotation = NULL,
                                       intensity = 'PeakArea',
                                       filter_with_Qvalue = TRUE,
                                       qvalue_cutoff = 0.01,
                                       useUniquePeptide = TRUE,
                                       fewMeasurements="remove",
                                       removeProtein_with1Feature = FALSE,
                                       summaryforMultipleRows=max){

    if (is.null(fewMeasurements)) {
        stop('** Please select \'remove\' or \'keep\' for \'fewMeasurements\'.')
    }
    
    if (!is.element(fewMeasurements, c('remove', 'keep'))) {
        stop('** Please select \'remove\' or \'keep\' for \'fewMeasurements\'.')
    }
    
    ## Check correct option or input
    requiredinput.general <- c('F.FrgLossType', 'F.ExcludedFromQuantification',
                               'PG.ProteinGroups', 'EG.ModifiedSequence', 'FG.Charge',
                               'F.FrgIon', 
                               'R.FileName', 'EG.Qvalue')
    
    requiredinput.int <- c('F.PeakArea', 'F.NormalizedPeakArea')
    
    requiredinput.charge <- c('F.Charge', 'F.FrgZ')
    
    ################################
    ## general input
    ################################
    
    if (!all( requiredinput.general %in% colnames(input))) {
        
        misssing.col <- requiredinput.general[!requiredinput.general %in% colnames(input)]
        
        stop(paste0("** Please check the required input. The required input needs : ", 
                    toString(missing.col)))
    }
    
    ## intensity columns
    if (sum( requiredinput.int %in% colnames(input) ) == 0) {
        
        stop(paste0("** Please check the required input. The required input needs at least one of ", 
                   toString(requiredinput.int)))
    }
    
    ## general input
    if (sum( requiredinput.charge %in% colnames(input) ) == 0) {
        
        stop(paste0("** Please check the required input. The required input needs at least one of '", 
                    toString(requiredinput.charge)))
    }
    
    ## get annotation
    if (is.null(annotation)) {
        annotinfo <- unique(input[, c("R.FileName", "R.Condition", "R.Replicate")])	
        colnames(annotinfo) <- c('Run', 'Condition', 'BioReplicate')
    } else {
        ## check annotation
        required.annotation <- c('Condition', 'BioReplicate', 'Run')
        
        if (!all(required.annotation %in% colnames(annotation))) {
            
            missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))
            
            stop(paste("**", toString(required.annotation[missedAnnotation]), 
                       "is not provided in Annotation. Please check the annotation file.",
                       "'Run' will be matched with 'R.FileName' "))
        } else {
            annotinfo <- annotation
        }
    }
    
    ## check annotation information
    ## Each Run should has unique information about condition and bioreplicate
    check.annot <- xtabs(~Run, annotinfo)
    if ( any(check.annot > 1) ) {
        stop('** Please check annotation. Each MS run can\'t have multiple conditions or BioReplicates.' )
    }
    
    
    ##############################
    ## 1. loss type : use only 'no loss'
    ##############################
    input <- input[input$F.FrgLossType == 'noloss', ]
  
    ##############################
    ## 2. use only 'F.ExcludedFromQuantification == False' : XIC quality
    ##############################
    if (is.logical(input$F.ExcludedFromQuantification)) {
        input$F.ExcludedFromQuantification <- factor(input$F.ExcludedFromQuantification)
        input$F.ExcludedFromQuantification <- factor(input$F.ExcludedFromQuantification,
                                                    labels = c('False', 'True'))
    }
    
    if (!all(unique(input$F.ExcludedFromQuantification) %in% c('False', 'True'))) {
        stop( paste("** Please check the column called F.ExcludedFromQuantification. Only False or True are allowed in this column."))
    }
    
    input <- input[input$F.ExcludedFromQuantification == 'False', ]
  
    ##############################
    ## 3. get useful subset of column
    ##############################
    if (is.element('F.Charge', colnames(input))) {
        f.charge <- 'F.Charge'
    } else if (is.element('F.FrgZ', colnames(input))) {
        f.charge <- 'F.FrgZ'
    } else {
        f.charge <- NULL
    }
    
    if (is.element('PG.Qvalue', colnames(input))) {
        pg.qvalue <- 'PG.Qvalue'
    } else if(is.element('PG.Qvalue', colnames(input))) {
        pg.qvalue <- 'PG.Qvalue'
    } else {
        pg.qvalue <- NULL
    }
    
    subsetcolumn <- c('PG.ProteinGroups', 'EG.ModifiedSequence', 'FG.Charge',
                    'F.FrgIon', f.charge,
                    'R.FileName', 
                    'EG.Qvalue', pg.qvalue)
  
    if (intensity == 'NormalizedPeakArea') {
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
    colnames(input)[colnames(input) == 'R.FileName'] <- 'Run'
    colnames(input)[colnames(input) == 'F.PeakArea'] <- 'Intensity'
    colnames(input)[colnames(input) == 'F.NormalizedPeakArea'] <- 'Intensity'
    colnames(input)[colnames(input) == 'EG.Qvalue'] <- 'Qvalue'
  
    ##############################
    ## 4. filter by Qvalue
    ##############################

    ## protein FDR
    if (is.element('PG.Qvalue', colnames(input))) {
        input[!is.na(input$PG.Qvalue) & input$PG.Qvalue > 0.01, "Intensity"] <- NA
        message('** Intensities with great than 0.01 in PG.Qvalue are replaced with NA.')
        
        input <- input[, -which(colnames(input) %in% 'PG.Qvalue')]
        
    }
    
    ## precursor qvalue
    if (filter_with_Qvalue) {
    
        if (!is.element(c('Qvalue'), colnames(input))) {
      
            stop('** EG.Qvalue column is needed in order to filter out by Qvalue. Please add EG.Qvalue column in the input.')
      
        } else {
      
            ## when qvalue > qvalue_cutoff, replace with zero for intensity
            input[!is.na(input$Qvalue) & input$Qvalue > qvalue_cutoff, "Intensity"] <- 0
      
            message(paste0('** Intensities with great than ', qvalue_cutoff, ' in EG.Qvalue are replaced with zero.'))
        }
    }
  
    ##############################
    ## 5. remove featuares with all na or zero
    ## some rows have all zero values across all MS runs. They should be removed.
    ##############################
  
    input$fea <- paste(input$PeptideSequence,
                        input$PrecursorCharge,
                        input$FragmentIon,
                        input$ProductCharge,
                        sep="_")
  
    inputtmp <- input[!is.na(input$Intensity) & input$Intensity > 1, ]
    
    count <- inputtmp %>% group_by(fea) %>% summarise(length=length(Intensity))
    
    ## get feature with all NA or zeros
    getfea <- count[count$length > 0, 'fea']
    
    if (nrow(getfea) > 0) {
        
        nfea.remove <- length(unique(input$fea))-nrow(getfea)
        input <- input[which(input$fea %in% getfea$fea), ]

        message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
        
    }
    
    rm(inputtmp)
    
    ################################################
    ## 6. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    if (useUniquePeptide) {
    
        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")]) ## Protein.group.IDs or Sequence
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)
    
        ## count how many proteins are assigned for each peptide
        structure <- pepcount %>% group_by(PeptideSequence) %>% summarise(length=length(ProteinName))
        
        remove_peptide <- structure[structure$length != 1, ]
        
        ## remove the peptides which are used in more than one protein
        if(nrow(remove_peptide) != 0){
            input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]
            
            message('** Peptides, that are used in more than one proteins, are removed.')
        } else {
            message('** All peptides are unique peptides in proteins.')
        }
    
        rm(structure)
        rm(remove_peptide)
    }
  
    ##############################
    ##  7. remove features which has 1 or 2 measurements across runs
    ##############################
    if (fewMeasurements == "remove"){
    
        ## it is the same across experiments. # measurement per feature. 
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        count_measure <- xtabs( ~fea, xtmp)
    
        remove_feature_name <- count_measure[count_measure < 3]
    
        if (length(remove_feature_name) > 0) {
            
            input <- input[-which(input$fea %in% names(remove_feature_name)), ]
            
            message(paste0('** ', length(remove_feature_name), ' features have 1 or 2 intensities across runs and are removed.'))
        }
    }
  
    ##############################
    ## 8. remove proteins with only one peptide and charge per protein
    ##############################
  
    if (removeProtein_with1Feature) {
        ## remove protein which has only one peptide
        tmp <- unique(input[, c("ProteinName", 'fea')])
        tmp$ProteinName <- factor(tmp$ProteinName)
        count <- xtabs( ~ ProteinName, data=tmp)
        lengthtotalprotein <- length(count)
    
        removepro <- names(count[count <= 1])
    
        if (length(removepro) > 0) {
      
            input <- input[-which(input$ProteinName %in% removepro), ]
            
            message(paste0("** ", length(removepro), 
                           ' proteins, which have only one feature in a protein, are removed among ', 
                           lengthtotalprotein, ' proteins.'))
        } else {
            message("** All proteins have at least two features.")
        }
    }

    ##############################
    ## 9. remove multiple measurements per feature and run
    ##############################
  
    count <- aggregate(Intensity ~ fea, data=input, FUN=length)
  
    ## if any feature has more number of total MS runs, 
    if (any(unique(count$Intensity) > length(unique(input$Run)))) {
    
        ## maximum or sum up abundances among intensities for identical features within one run
        input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon + ProductCharge ~ Run, data=input, 
                          value.var='Intensity', 
                          fun.aggregate=summaryforMultipleRows, na.rm=T,
                          fill=NA_real_) 
    
        ## reformat for long format
        input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon', 'ProductCharge'))
        colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    
        message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
    
    } else {
        ## remove column, named as 'fea'
        input <- input[, -which(colnames(input) %in% c('fea', 'Qvalue'))]
        message('** No multiple measurements in a feature and a run.')
    }
    
    ##############################
    ## 10. merge annotation
    ##############################
    input <- merge(input, annotinfo, all=TRUE)

    input.final <- data.frame("ProteinName" = input$ProteinName,
                              "PeptideSequence" = input$PeptideSequence,
                              "PrecursorCharge" = input$PrecursorCharge,
                              "FragmentIon" = input$FragmentIon,
                              "ProductCharge" = input$ProductCharge,
                              "IsotopeLabelType" = "L",
                              "Condition" = input$Condition,
                              "BioReplicate" = input$BioReplicate,
                              "Run" = input$Run,
                              "Intensity" = input$Intensity)
    
    input <- input.final
    input$ProteinName <- factor(input$ProteinName)
    
    rm(input.final)
  
	return(input)
}


