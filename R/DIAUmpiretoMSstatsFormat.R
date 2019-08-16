
## raw.frag : FragSummary
## raw.pep : PeptideSummary with selected_fragmnet
## raw.pro : ProteinSummary with selected_peptide

## useUniquePeptide : remove peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
## summaryforMultipleRows : max or sum - when there are multiple measurements for certain feature and certain fun, use highest or sum of all.
## fewMeasurements : if 1 or 2 measurements across runs per feature, 'remove' will remove those featuares. It can affected for unequal variance analysis.

#' @export
#' @importFrom dplyr %>% left_join semi_join group_by summarise
#' @importFrom tidyr separate_rows
#' @import stringr
DIAUmpiretoMSstatsFormat <- function(raw.frag, raw.pep, raw.pro, 
                                annotation,
                                useSelectedFrag = TRUE,
                                useSelectedPep = TRUE,
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
            ### annotation.txt : Condition, BioReplicate, Run, 
            annot <- annotation
        }
    }
    
    ## Each Run should has unique information about condition and bioreplicate
    check.annot <- xtabs(~Run, annot)
    if ( any(check.annot > 1) ) {
        stop('** Please check annotation. Each MS run can\'t have multiple conditions or BioReplicates.' )
    }
    
    ########################
    ## 1.  get selected frag from DIA-Umpire
    ########################
    
    ## need to check whether 'Selected_fragments' column is available or not.
    if (!any(is.element(colnames(raw.pep), 'Selected_fragments'))) {
        stop('** Selected_fragments column is required. Please check it.')
    }
    
    ## Selected_fragments
    raw.pep2 <- raw.pep[, c('Peptide.Key', 'Proteins', 'Selected_fragments')]
    
    ## remove empty protein name
    raw.pep2 <- raw.pep2[raw.pep2$Proteins != '', ]
    #raw.pep2$Proteins <- gsub(';', '', raw.pep2$Proteins) ## be careful for multiple protein ids
    raw.pep2$Selected_fragments <- as.character(raw.pep2$Selected_fragments)
    
    ## get infor for selected fragment
    raw.pep3 <- raw.pep2 %>%
        separate_rows(Selected_fragments, sep="[|]")
    ## raw.pep3 has selected fragment information
    
    ########################
    ## 2. get selected peptides from DIA-Umpire
    ########################
    
    if (!any(is.element(colnames(raw.pro), 'Selected_peptides'))) {
        stop('** Selected_peptides column is required. Please check it.')
    }
    
    ## Selected_peptide
    raw.pro2 <- raw.pro[, c('Protein.Key', 'Selected_peptides')]
    
    ## remove empty protein name
    raw.pro2 <- raw.pro2[raw.pro2$Protein.Key != '', ]
    raw.pro2$Selected_peptides <- as.character(raw.pro2$Selected_peptides)
    
    ## get info for selected peptides
    raw.pro3 <- raw.pro2 %>%
        separate_rows(Selected_peptides, sep="[|]")

    ## top 6 peptides are selected, raw.pro3 has selected peptide information
    
    ########################
    ## 3. subtract the rows with selected peptides and fragments
    ########################
    colnames(raw.pro3) <- c('Protein', 'Peptide')
    colnames(raw.pep3) <- c('Peptide', 'Protein', 'Fragment')
    
    ## some proteins have no selected peptides
    ## ex : raw.pro3[raw.pro3$Peptide == '', ]
    ## !! if empty for selected peptide or fragments: 
    ## all zero or NA if there is no peptide information in protein
    ## remove all, none of peptide are selected
    ## conclusion, if 'selected_xxx' is empty, no selected.

    raw.pro3 <- raw.pro3[raw.pro3$Peptide != '', ]
    raw.pro3$Peptide <- as.character(raw.pro3$Peptide)
    
    raw.pep3 <- raw.pep3[raw.pep3$Fragment != '', ]
    raw.pep3$Peptide <- as.character(raw.pep3$Peptide)
    
    if (useSelectedFrag & useSelectedPep) {
        
        raw.pep3$Protein <- gsub(';', '', raw.pep3$Protein)
        
        ## remove the shared peptides
        count <- xtabs(~Peptide, raw.pro3)
        notuni.pep <- names(count)[count > 1]
        
        if (length(notuni.pep) > 0){
            ## remove not-unique peptide in proteien quantification
            raw.pro3 <- raw.pro3[-which(raw.pro3$Peptide %in% notuni.pep), ]
            raw.pep3 <- raw.pep3[-which(raw.pep3$Peptide %in% notuni.pep), ]
        } 
        
        ## get selected peptides only in raw.pro3
        raw.pep3 <- raw.pep3[which(raw.pep3$Peptide %in% unique(raw.pro3$Peptide)), ]
        unique(nchar(as.character(unique(raw.pep3$Protein))))
        
        ## The number of peptides are the same in raw.pep3 or raw.pro3.
        ## but the number of proteins are different.
        ## usd these selected peptides, and use protein id in raw.pro3. because pro3 assign one protein id for each peptides
        
        ## now raw.pep3 can be used for final list
        ## but need to assign protein id from raw.pro3. It should be updated in raw.frag.
        ## raw.frag can have multiple protein ids.
        raw.pep3 <- raw.pep3[, -which(colnames(raw.pep3) %in% 'Protein')]
        
        ## make the final list of selected proteins and peptides
        raw.final <- left_join(raw.pep3, raw.pro3, by=c('Peptide'))
        raw.final$Protein <- as.character(raw.final$Protein)
        
        message('** Get the selected fragments and peptides.')
        
        #length(unique(raw.final$Peptide)) # 12975
        #length(unique(raw.final$Protein)) # 3627
        #length(unique(raw.pep3$Peptide))
        #length(unique(raw.pro3$Protein))
        
    } else if (useSelectedFrag & !useSelectedPep) { # only use selected fragment, and use all peptides
        
        ## can use raw.pep3 only.
        ## but doublecheck protein id : protein id in raw.frag is always one id.
        ## need to find shared peptides
        count.pro <- str_count(raw.pep3$Protein, '\\;')

        ## remove not-unique peptide in proteien quantification
        raw.pep3 <- raw.pep3[count.pro < 2, ]
        
        raw.pep3$Protein <- gsub(';', '', raw.pep3$Protein)
        
        raw.final <- raw.pep3
        
        message('** Got the selected fragments.')
        
    } else if (!useSelectedFrag & !useSelectedPep) {
        
        stop('** MSstats recommends to use at least selected fragments.')
    }
    

    ############################
    ## 4. subtract fragment data with selected peptide and fragment
    ############################
    
    ## 0.1  check protein id
    raw.frag$Protein <- as.character(raw.frag$Protein)
    
    ## 0.2  reduce column
    inputlevel <- 'Intensity'
    
    ## columns which include 'inputlevel
    int.column <- colnames(raw.frag)[grep(inputlevel, colnames(raw.frag))]
    deft.column <- c('Fragment.Key', 'Protein','Peptide','Fragment')
    
    ## get subset of data first
    raw <- raw.frag[, which(colnames(raw.frag) %in% c(int.column, deft.column))]
    
    ## 0.3 update fragment info
    raw$Fragment <- as.character(raw$Fragment)
    
    if (length(grep('\\+', raw$Fragment[1])) < 1) { ## if Fragment column does not have charge
        
        raw$Fragment.Key <- as.character(raw$Fragment.Key)
        raw$ion <- substr(raw$Fragment.Key, start=nchar(raw$Fragment.Key), stop = nchar(raw$Fragment.Key))
        raw$Fragment <- paste(raw$Fragment, raw$ion, sep="_")
        
        raw <- raw[, -which(colnames(raw) %in% c('ion'))]
        
    } else {
        raw$Fragment <- gsub('\\+', '_', raw$Fragment)
    }
    
    ## remove extra column
    raw <- raw[, -which(colnames(raw) %in% c('Fragment.Key'))]
    
    ## 0.4 get selected fragment and peptide only
    raw$Peptide <- as.character(raw$Peptide)
    raw$Fragment <- as.character(raw$Fragment)
    raw$Fragment <- as.character(raw$Fragment)
    
    raw.selected <- semi_join(raw, raw.final, by=c('Protein', 'Peptide', 'Fragment'))
    
    message('** Extract the data from selected fragments and/or peptides.')
    
    ############################
    ## 5. reformat
    ############################
    ## change column names
    colnames(raw.selected)[colnames(raw.selected) == 'Protein'] <- 'ProteinName'
    colnames(raw.selected)[colnames(raw.selected) == 'Peptide'] <- 'PeptideSequence'
    colnames(raw.selected)[colnames(raw.selected) == 'Fragment'] <- 'FragmentIon'
    
    ## make long format
    raw.l <- melt(raw.selected, id.vars=c('ProteinName', 'PeptideSequence', 'FragmentIon'),
                  variable.name='Run', value.name='Intensity')
    
    ## remove suffix in the Run
    raw.l$Run <- gsub(paste("_", inputlevel, sep=""), '', raw.l$Run)
    
    ## fill in other column
    raw.l$PrecursorCharge <- NA
    raw.l$ProductCharge <- NA
    raw.l$IsotopeLabelType <- 'light'
    
    ## 4. add annotation
    input <- merge(raw.l, annot, by='Run')

    rm(raw.l)
    rm(raw.selected)
    
    ##############################
    ## 6. remove featuares with all na or zero
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
        
        nfea.remove <- length(unique(input$fea)) - nrow(getfea)
        input <- input[which(input$fea %in% getfea$fea), ]
        
        message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
        
    }
    
    rm(inputtmp)
    
    ##############################
    ## 7. remove features which has 1 or 2 measurements across runs
    ##############################
    if (fewMeasurements=="remove") {
        
        ## it is the same across experiments. # measurement per feature. 
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        count_measure <- xtabs(~ fea, xtmp)
        
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
            
            message("All proteins have at least two features.")
            
        }
    }
    
    ##############################
    ## 9. remove multiple measurements per feature and run
    ##############################
    
    count <- aggregate(Intensity ~ fea, data=input, FUN=length)
    
    ## if any feature has more number of total MS runs, 
    if (any(unique(count$Intensity) > length(unique(input$Run)))){
        
        annotation <- unique(input[, c("Condition", "BioReplicate", "Run")])
        
        ## maximum or sum up abundances among intensities for identical features within one run
        input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon + ProductCharge ~ Run, data=input, 
                          value.var='Intensity', 
                          fun.aggregate=summaryforMultipleRows, 
                          fill=NA_real_) 
        
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
        input <- input[, -which(colnames(input) %in% c('fea'))]
        message('** No multiple measurements in a feature and a run.')
    }
    
    input$IsotopeLabelType <- 'L'
    input$ProteinName <- factor(input$ProteinName)

	return(input)
}


