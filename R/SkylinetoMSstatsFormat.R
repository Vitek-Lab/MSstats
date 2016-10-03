
## Pre-processing for Skyline output
## columns from Skyline : ProteinName, PeptideSequence, PeptideModifiedSequence,
##                        PrecursorCharge, PrecursorMz, FragmentIon, ProductCharge, ProductMz, IsotopeLabelType,   
##                        Condition, BioReplicate, FileName, Area, StandardType, Truncated, DetectionQValue 

#' @export
SkylinetoMSstatsFormat <- function(input, 
                                   annotation = NULL,
                                   removeiRT = TRUE, 
                                   removeProtein_with1Peptide = FALSE,
                                   filter_with_Qvalue = TRUE,
                                   qvalue_cutoff = 0.01){
  
  ##############################
  ### 1. Rename column names
  ##############################
  
  if( is.element(c('Protein.Name'), colnames(input)) ){
    colnames(input)[colnames(input) == 'Protein.Name'] <- 'ProteinName'
  }
  
  if( is.element(c('Peptide.Modified.Sequence'), colnames(input)) ){
    colnames(input)[colnames(input) == 'Peptide.Modified.Sequence'] <- 'PeptideModifiedSequence'
  }
  
  if( is.element(c('Precursor.Charge'), colnames(input)) ){
    colnames(input)[colnames(input) == 'Precursor.Charge'] <- 'PrecursorCharge'
  }
  
  if( is.element(c('Fragment.Ion'), colnames(input)) ){
    colnames(input)[colnames(input) == 'Fragment.Ion'] <- 'FragmentIon'
  }
  
  if( is.element(c('Product.Charge'), colnames(input)) ){
    colnames(input)[colnames(input) == 'Product.Charge'] <- 'ProductCharge'
  }
  
  if( is.element(c('Isotope.Label.Type'), colnames(input)) ){
    colnames(input)[colnames(input) == 'Isotope.Label.Type'] <- 'IsotopeLabelType'
  }
  
  if( is.element(c('File.Name'), colnames(input)) ){
    colnames(input)[colnames(input) == 'File.Name'] <- 'FileName'
  }
  
  if( is.element(c('Standard.Type'), colnames(input)) ){
    colnames(input)[colnames(input) == 'Standard.Type'] <- 'StandardType'
  }
  
  
  ## remove PeptideSequence,
  if( sum(is.element(c('PeptideSequence', 'PeptideModifiedSequence'), colnames(input))) == 2 ){
    input <- input[, -which(colnames(input) %in% 'PeptideSequence'), ]    
  }

  ## use PeptideModifiedSequence for PeptideSequence,
  colnames(input)[colnames(input) == 'PeptideModifiedSequence'] <- 'PeptideSequence'
  ## replace 'FileName' with Run
  colnames(input)[colnames(input) == 'FileName'] <- 'Run'
  ## replace 'Area' with Intensity
  colnames(input)[colnames(input) == 'Area'] <- 'Intensity'
  
  ##############################
  ### 2. Remove decoy protein name
  ##############################
  
  proname <- unique(input$ProteinName)
  decoy1 <- proname[grep('DECOY', proname)] 
  decoy2 <- proname[grep('Decoys', proname)] 
  decoy <- c(as.character(decoy1), decoy2)
  
  if(length(decoy) > 0) {
    
    input <- input[-which(input$ProteinName %in% decoy), ]
    
    message('** Proteins, which names include DECOY, are removed.')
  }

  ##############################
  ### 3. Remove iRT proteins
  ##############################
  
  if( removeiRT ){
    irt <- unique(input$StandardType)
  
    if( sum(is.element(irt, 'iRT')) > 0 ){
    
      input <- input[-which(input$StandardType %in% c("iRT")), ]
      message('** iRT proteins/peptides are removed.')
    }
  }

  ##############################
  ### 4. class of intensity is factor, change it as numeric
  ##############################

  input$Intensity <- as.numeric(as.character(input$Intensity))
  
  
  ##############################
  ###  5. remove truncated peaks with NA
  ##############################
  
  if ( sum(input$Truncated == 'True') > 0 ){
    
    input[input$Truncated == "True", "Intensity"] <- NA
    message('** Truncated peaks are replaced with NA.')
  }
  
  ##############################
  ###  6. Sum for isotopic peaks per peptide and precursor charge for DDA
  ##############################

  DDA <- FALSE
  
  ## check whether the dataset for DDA or not
  input$FragmentIon <- factor(input$FragmentIon)
  checkDDA <- setdiff(c('precursor', 'precursor [M+1]', 'precursor [M+2]'), levels(input$FragmentIon))
  
  if( length(checkDDA) < 3 ){
    
    DDA <- TRUE
    ## add the column for unique peptide and precursor
    input$pepprecursor <- paste(input$PeptideSequence, input$PrecursorCharge, sep="_")
    input <- input[!is.na(input$Intensity), ]

    ## sum of mooisotopic peaks
    ## zero is kept as zero, missing rows replace with NA (fill option)
    data_w <- dcast( pepprecursor ~ Run, data=input, value.var='Intensity', fun.aggregate=function(x) sum(x, na.rm=TRUE), fill=NA_real_) 
    
    ## make long format
    newdata <- melt(data_w, id.vars=c('pepprecursor'))
    colnames(newdata)[colnames(newdata) %in% c("variable","value")] <- c('Run','Intensity')
    
    ## assignn protein name
    uniinfo <- unique(input[, c("ProteinName", "PeptideSequence", "PrecursorCharge", "pepprecursor")])	
    
    ## get annotation
    if( is.null(annotation) ){
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
    
    input <- input[, c(4,5,6,9,10,11,7,8,1,3)]
    
    message('** For DDA datasets, three isotopic peaks per feature and run are summed.')
    
  }
  

  ##############################
  ###  7. remove proteins with only one peptide and charge per protein
  ##############################
	
	if(removeProtein_with1Peptide){
	  ######## remove protein which has only one peptide
	  input$feature <- paste(input$PeptideSequence, input$PrecursorCharge, input$FragmentIon, input$ProductCharge, sep="_")
	  
	  tmp <- unique(input[, c("ProteinName", 'feature')])
	  tmp$ProteinName <- factor(tmp$ProteinName)
	  count <- xtabs( ~ ProteinName, data=tmp)
    lengthtotalprotein <- length(count)
    
	  removepro <- names(count[count <= 1])
	  
	  if (length(removepro) > 0) {
	    
	    input <- input[-which(input$ProteinName %in% removepro), ]
	    message(paste("*** ", length(removepro), ' proteins, which have only one feature in a protein, are removed among ', lengthtotalprotein, ' proteins.', sep=""))
	  }
	  
	  input <- input[, -which(colnames(input) %in% c('feature'))]
	}
  
  ##############################
  ###  8. filter by Qvalue
  ##############################
  
  if( !DDA & filter_with_Qvalue ){
    
    if(! is.element(c('DetectionQValue'), colnames(input)) ){
      
      stop('** DetectionQValue column is needed in order to filter out by Qvalue. Please add DectionQValue column in the input.')
      
    } else {
    
      ## make Q value as numeric
      input$DetectionQValue <- as.numeric(as.character(input$DetectionQValue))
    
      ## when qvalue > qvalue_cutoff, replace with zero for intensity
      input[!is.na(input$DetectionQValue) & input$DetectionQValue > qvalue_cutoff, "Intensity"] <- 0
    
      message(paste('Intensities with great than ', qvalue_cutoff, ' in DetectionQValue are replaced with zero.', sep=''))
    }
  }
  
  input$ProteinName <- input$ProteinName
  
	return(input)
}




