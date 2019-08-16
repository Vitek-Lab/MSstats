
## Pre-processing for Skyline output
## columns from Skyline : ProteinName, PeptideSequence, PeptideModifiedSequence,
##                        PrecursorCharge, PrecursorMz, FragmentIon, ProductCharge, ProductMz, IsotopeLabelType,   
##                        Condition, BioReplicate, FileName, Area, StandardType, Truncated, DetectionQValue 

#' @export
SkylinetoMSstatsFormat <- function(input, 
                                   annotation = NULL,
                                   removeiRT = TRUE, 
                                   filter_with_Qvalue = TRUE,
                                   qvalue_cutoff = 0.01,
                                   useUniquePeptide = TRUE,
                                   fewMeasurements="remove",
                                   removeOxidationMpeptides = FALSE,
                                   removeProtein_with1Feature = FALSE){
  
    if (is.null(fewMeasurements)) {
        stop('** Please select \'remove\' or \'keep\' for \'fewMeasurements\'.')
    }
    
    if (!is.element(fewMeasurements, c('remove', 'keep'))) {
        stop('** Please select \'remove\' or \'keep\' for \'fewMeasurements\'.')
    }
    
    ##############################
    ## 1. Rename column names
    ##############################
  
    ## replace '.' between words with no spzce
    colnames(input) <- gsub('\\.', '', colnames(input))

    ## use PeptideModifiedSequence for PeptideSequence,
    ## if there are both peptidesequence and peptidemodifiedsequence,
    ## use peptidemodifiedsequence.
    if (any(is.element(colnames(input), 'PeptideSequence')) & any(is.element(colnames(input), 'PeptideModifiedSequence'))){
        input <- input[, -which(colnames(input) %in% 'PeptideSequence')]
    }
    
    if (any(is.element(colnames(input), 'PeptideModifiedSequence'))){
        colnames(input)[colnames(input) == 'PeptideModifiedSequence'] <- 'PeptideSequence'
    }
    
    ## replace 'FileName' with Run
    colnames(input)[colnames(input) == 'FileName'] <- 'Run'
    ## replace 'Area' with Intensity
    colnames(input)[colnames(input) == 'Area'] <- 'Intensity'
    
    
    ##############################
    ## 1.1. check annotation information
    ##############################
    ## get annotation
    if (is.null(annotation)) {
        annotinfo <- unique(input[, c("Run", "Condition", 'BioReplicate')])	
    } else {
        annotinfo <- annotation
    }
    
    ## Each Run should has unique information about condition and bioreplicate
    check.annot <- xtabs(~Run, annotinfo)
    if ( any(check.annot > 1) ) {
        stop('** Please check annotation. Each MS run can\'t have multiple conditions or BioReplicates.' )
    }
    
    ##############################
    ## 2. Remove decoy protein name
    ##############################
  
    proname <- unique(input$ProteinName)
    decoy1 <- proname[grep('DECOY', proname)] 
    decoy2 <- proname[grep('Decoys', proname)] 
    decoy <- c(as.character(decoy1), decoy2)
  
    if (length(decoy) > 0) {
    
        input <- input[-which(input$ProteinName %in% decoy), ]
    
        message('** Proteins, which names include DECOY, are removed.')
    }

    ##############################
    ## 3. Remove iRT proteins
    ##############################
  
    if (removeiRT) {
        irt <- unique(input$StandardType)
  
        if (sum(is.element(irt, 'iRT')) > 0) {
    
            input <- input[-which(input$StandardType %in% c("iRT")), ]
            message('** iRT proteins/peptides are removed.')
        }
    }

    ################################################
    ## 4. remove the peptides including M sequence
    ################################################

    if (removeOxidationMpeptides) {
        remove_oxim_sequence <- unique(input[grep("+16", input$PeptideSequence), "PeptideSequence"])
    
        if (length(remove_oxim_sequence) > 0) {
            input <- input[-which(input$PeptideSequence %in% remove_oxim_sequence), ]
        }
    
        message('** Peptides including M[+16] in the sequence are removed.')
    }
  
    ################################################
    ## 5. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    if (useUniquePeptide) {
    
        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")]) ## Protein.group.IDs or Sequence
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)
    
        ## count how many proteins are assigned for each peptide
        structure <- pepcount %>% group_by(PeptideSequence) %>% summarise(length=length(ProteinName))
        
        remove_peptide <- structure[structure$length != 1, ]
    
        ## remove the peptides which are used in more than one protein
        if (nrow(remove_peptide) != 0){
            input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]
        }
    
        message('** Peptides, that are used in more than one proteins, are removed.')
    }
  
    ##############################
    ## 6. class of intensity is factor, change it as numeric
    ##############################

    input$Intensity <- as.numeric(as.character(input$Intensity))
  
    ##############################
    ##  7. remove truncated peaks with NA
    ##############################
  
    if (is.element('True', input$Truncated)) {
        if (sum(!is.na(input$Truncated) & input$Truncated == 'True') > 0) {
      
            input[!is.na(input$Truncated) & input$Truncated == "True", "Intensity"] <- NA
            message('** Truncated peaks are replaced with NA.')
        }
    }
  
    if (is.element(TRUE, input$Truncated)) {
    
        if (sum(!is.na(input$Truncated) & input$Truncated) > 0) {
    
            input[!is.na(input$Truncated) & input$Truncated, "Intensity"] <- NA
            message('** Truncated peaks are replaced with NA.')
        }
    }
  
    ##############################
    ##  8. DDA : Sum for isotopic peaks per peptide and precursor charge for DDA
    ##############################

    DDA <- FALSE
  
    ## check whether the dataset for DDA or not
    input$FragmentIon <- factor(input$FragmentIon)
    checkDDA <- setdiff(c('precursor', 'precursor [M+1]', 'precursor [M+2]'), levels(input$FragmentIon))
    
    ## need to check mixed in fragmention column.
    any.fragment <- setdiff( levels(input$FragmentIon), c('precursor', 'precursor [M+1]', 'precursor [M+2]'))
    any.precursor3 <- intersect( levels(input$FragmentIon), c('precursor', 'precursor [M+1]', 'precursor [M+2]'))
    
    ## if there are fragment ion and also have any 'precursor', it is the issue.
    if (length(any.fragment) > 0 & length(any.precursor3) > 0) {
        stop("** Please check precursors information. If your experiment is DIA, please remove the precursors. If your experiments is DDA, please check the precursor information.")
    }
    
    if (length(checkDDA) < 3) {
    
        DDA <- TRUE
        ## add the column for unique peptide and precursor
        input$pepprecursor <- paste(input$PeptideSequence, input$PrecursorCharge, sep="_")
        input <- input[!is.na(input$Intensity), ]
        
        ## keep StandardType information
        standard.info <- input[, c('ProteinName', 'PeptideSequence', 
                                   'PrecursorCharge', 'StandardType')]
        standard.info <- standard.info[!duplicated(standard.info), ]
        
        ## sum of mooisotopic peaks
        ## zero is kept as zero, missing rows replace with NA (fill option)
        data_w <- dcast(pepprecursor ~ Run, data=input, 
                        value.var='Intensity', 
                        fun.aggregate=function(x) sum(x, na.rm=TRUE), 
                        fill=NA_real_) 
    
        ## make long format
        newdata <- melt(data_w, id.vars=c('pepprecursor'))
        colnames(newdata)[colnames(newdata) %in% c("variable","value")] <- c('Run','Intensity')
    
        ## assignn protein name
        uniinfo <- unique(input[, c("ProteinName", "PeptideSequence", "PrecursorCharge", "pepprecursor")])	
    
        ## get annotation
        if (is.null(annotation)) {
            annotinfo <- unique(input[, c("Run", "Condition", 'BioReplicate')])	
        } else {
            annotinfo <- annotation
        }
    	
        input <- merge(newdata, uniinfo, by="pepprecursor")
    
        ## assign the annotation
        ## merge it by Run
        input <- merge(input, annotinfo, by="Run")
    
        ## add other required information
        input$FragmentIon <- "sum"
        input$ProductCharge <- NA
        input$IsotopeLabelType <- "L"
        
        ## merge standard type information
        input <- merge(input, standard.info, all=TRUE)
        
        input.final <- data.frame("ProteinName" = input$ProteinName,
                                  "PeptideSequence" = input$PeptideSequence,
                                  "PrecursorCharge" = input$PrecursorCharge,
                                  "FragmentIon" = input$FragmentIon,
                                  "ProductCharge" = input$ProductCharge,
                                  "IsotopeLabelType" = input$IsotopeLabelType,
                                  "Condition" = input$Condition,
                                  "BioReplicate" = input$BioReplicate,
                                  "Run" = input$Run,
                                  "Intensity" = input$Intensity,
                                  "StandardType" = input$StandardType)
        
        if (any(is.element(colnames(input), 'Fraction'))) {
            input.final <- data.frame(input.final,
                                      "Fraction" = input$Fraction)
        }
        
        input <- input.final
        rm(input.final)
    
        message('** For DDA datasets, three isotopic peaks per feature and run are summed.')
    
    }
  
    
    ##############################
    ## 9. DIA : filter by Qvalue
    ##############################
    
    if (!DDA & filter_with_Qvalue){
        
        if (!is.element(c('DetectionQValue'), colnames(input))) {
            
            stop('** DetectionQValue column is needed in order to filter out by Qvalue. Please add DectionQValue column in the input.')
            
        } else {
            
            ## make Q value as numeric
            input$DetectionQValue <- as.numeric(as.character(input$DetectionQValue))
            
            ## when qvalue > qvalue_cutoff, replace with zero for intensity
            input[!is.na(input$DetectionQValue) & input$DetectionQValue > qvalue_cutoff, "Intensity"] <- 0
            
            message(paste0('** Intensities with great than ', qvalue_cutoff, 
                           ' in DetectionQValue are replaced with zero.'))
        }
    }
    
    ##############################
    ## 10. remove featuares with all na or zero
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
    
    if (nrow(getfea) > 0){
        
        nfea.remove <- length(unique(input$fea))-nrow(getfea)
        input <- input[which(input$fea %in% getfea$fea), ]
        
        message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
    }
    
    rm(inputtmp)
    
    ##############################
    ##  11. remove features which has 1 or 2 measurements across runs
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
    ##  12. remove proteins with only one feature per protein
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
    
    input <- input[, -which(colnames(input) %in% c('fea'))]
  
    ##############################
    ##  13. if annotation is missing,
    ##############################
    missing.annotation <- any( is.na(input$Condition) | is.na(input$BioReplicate) )
    
    if (missing.annotation & is.null(annotation)) {
        stop('** Please check annotation for Condition and BioReplicat column. There is missing information.')	
    } else if (missing.annotation & !is.null(annotation)) {
        annotinfo <- annotation
        
        input <- input[, -which(colnames(input) %in% c('Condition', 'BioReplicate'))]
        
        ## assign the annotation
        ## merge it by Run
        input <- merge(input, annotinfo, by="Run")
        
        input.final <- data.frame("ProteinName" = input$ProteinName,
                                  "PeptideModifiedSequence" = input$PeptideSequence,
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
                                      Fraction = input$Fraction)
        }
        
        if (any(is.element(colnames(input), 'DetectionQValue'))) {
            input.final <- data.frame(input.final,
                                      "DetectionQValue" = input$DetectionQValue)
        }
        
        input <- input.final
        rm(input.final)
    }
  
    input$ProteinName <- factor(input$ProteinName)
  
	return(input)
}

