
## Converter for openswath output

#' @export

OpenSWATHtoMSstatsFormat <- function(input,
                                     annotation = NULL,
                                     filter_with_mscore = TRUE,
                                     mscore_cutoff = 0.01,
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
    
    if (is.null(annotation)) {
        stop('** Please prepare \'annotation\' as one of input.')
    } else {
        
        ## check annotation
        required.annotation <- c('Condition', 'BioReplicate', 'Run')
        
        if (!all(required.annotation %in% colnames(annotation))) {
            
            missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))
            
            stop(paste("**", toString(required.annotation[missedAnnotation]), 
                       "is not provided in Annotation. Please check the annotation file."))
        } else {
            annotinfo <- annotation
        }
    }
  
    ## check annotation information
    ## get annotation
    ## Each Run should has unique information about condition and bioreplicate
    check.annot <- xtabs(~Run, annotinfo)
  
    if ( any(check.annot > 1) ) {
        stop('** Please check annotation. Each MS run can\'t have multiple conditions or BioReplicates.' )
    }

    ## Check correct option or input
    requiredinput.general <- c('ProteinName', 'FullPeptideName', 'Charge', 'Sequence',
                               'decoy', 'm_score',
                               'aggr_Fragment_Annotation', 'aggr_Peak_Area', 
                               'filename')
    
    ################################
    ## 1. check general input and use only required columns.
    ################################
    
    if (!all(requiredinput.general %in% colnames(input))) {
        
        misssing.col <- requiredinput.general[!requiredinput.general %in% colnames(input)]
        
        stop(paste0("** Please check the required input. The required input needs : ", 
                   toString(missing.col)))
    } else {
        input <- input[, colnames(input) %in% requiredinput.general]
    }
    
    ##############################
    ## 2. remove the decoys
    ##############################
    if (length(unique(input$decoy)) == 2) {
        
        ## decoy=1 means the decoys.
        input <- input[input$decoy == 0, ]
        input <- input[, -which(colnames(input) %in% 'decoy')]
        
        message("** Remove the decoys.")
    }
  
    ##############################
    ## 3. filter by mscore
    ##############################
    
    ## mscore
    if (filter_with_mscore) {
    
        if (!is.element(c('m_score'), colnames(input))) {
      
            stop('** m_score column is needed in order to filter out by m_scoe. Please add m_score column in the input.')
      
        } else {
      
            ## when mscore > mscore_cutoff, replace with zero for intensity
            input <- input[!is.na(input$m_score) & input$m_score <= mscore_cutoff, ] 
      
            message(paste0('** Features with great than ', mscore_cutoff, ' in m_score are removed.'))
        }
        
        input <- input[, -which(colnames(input) %in% 'm_score')]
    }
    
    ##############################
    ## 4. Make required long format - disaggregate : one row for each transition
    ##############################
    
    ## The columns "aggr_Fragment_Annotation" : separate by ';' and "aggr_Peak_Area" : separate by ';' 
    ## are disaggregated into the new columns "FragmentIon" and "Intensity". 
    
    ## 2020.10.05 : change factor to character
    input$aggr_Fragment_Annotation <- as.character(input$aggr_Fragment_Annotation)
    input$aggr_Peak_Area <- as.character(input$aggr_Peak_Area)
    
    input <- input %>%
        separate_rows(aggr_Fragment_Annotation, aggr_Peak_Area, sep="[;]")
    
    colnames(input)[colnames(input) == 'aggr_Fragment_Annotation'] <- 'FragmentIon'
    colnames(input)[colnames(input) == 'aggr_Peak_Area'] <- 'Intensity'
    
    ## use FullPeptideName as peptide sequence
    ##   Sequence       FullPeptideName
    ##  TAEICEHLKR  TAEIC(UniMod:4)EHLKR
    ##     SCTILIK     SC(UniMod:4)TILIK
    
    ## FullPeptideName -> PeptideSequence, 
    ## Charge -> PrecursorCharge, 
    
    colnames(input)[colnames(input) == 'FullPeptideName'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'Charge'] <- 'PrecursorCharge'
    colnames(input)[colnames(input) == 'filename'] <- 'Run'
    
    input <- input[, -which(colnames(input) %in% 'Sequence')]
    
    ## Unimod Identifier should be replaced from ":" to "_".
    input$PeptideSequence <- gsub(':', '_', input$PeptideSequence)
    input$FragmentIon <- gsub(':', '_', input$FragmentIon)
    
    ## class of intensity is character, change it as numeric
    input$Intensity <- as.numeric(input$Intensity)	
	
    ## there are incomplete rows, which potentially NA
    ## if negative and 0 values should be replaced with NA
    input[input$Intensity < 1, "Intensity"] <- 0
    
    ##############################
    ## 5. remove featuares with all na or zero
    ## some rows have all zero values across all MS runs. They should be removed.
    ##############################
  
    input$fea <- paste(input$PeptideSequence,
                        input$PrecursorCharge,
                        input$FragmentIon,
                        #input$ProductCharge,
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
        if (nrow(remove_peptide) != 0) {
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
    if (fewMeasurements=="remove") {
    
        ## it is the same across experiments. # measurement per feature. 
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        count_measure <- xtabs( ~fea, xtmp)
    
        remove_feature_name <- count_measure[count_measure < 3]
    
        if (length(remove_feature_name) > 0) {
            
            input <- input[-which(input$fea %in% names(remove_feature_name)), ]
            
            message(paste0('** ', length(remove_feature_name), 
                           ' features have 1 or 2 intensities across runs and are removed.'))
            
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
        input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon ~ Run, data=input, 
                          value.var='Intensity', 
                          fun.aggregate=summaryforMultipleRows, na.rm=T,
                          fill=0) 
    
        ## reformat for long format
        input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon'))
        colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    
        message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
    
    } else {
        
        ## still need to fill incomplete rows
        input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon ~ Run, data=input, 
                          value.var='Intensity', 
                          fill=0) 
        
        ## reformat for long format
        input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon'))
        colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
        
        message('** No multiple measurements in a feature and a run.')
    }
    
    
    ##############################
    ## 10. merge annotation
    ##############################
    input <- merge(input, annotinfo, by='Run', all=TRUE)
    
    ## fill in extra columns
    input.final <- data.frame("ProteinName" = input$ProteinName,
                              "PeptideSequence" = input$PeptideSequence,
                              "PrecursorCharge" = input$PrecursorCharge,
                              "FragmentIon" = input$FragmentIon,
                              "ProductCharge" = NA,
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


