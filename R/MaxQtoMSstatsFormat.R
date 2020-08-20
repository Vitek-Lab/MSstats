
## evidence : evidence.txt
## annotation : annotation.txt - Raw.file, Condition, BioReplicate, Run, (IsotopeLabelType)
## proteinGroups : proteinGroups.txt . if proteinGroups=NULL, use 'Proteins'. if not, use proteinGroups information for matching Protein group ID

## proteinID : Proteins or Leading.razor.protein
## useUniquePeptide : remove peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
## summaryforMultipleRows : max or sum - when there are multiple measurements for certain feature and certain fun, use highest or sum of all.
## fewMeasurements : if 1 or 2 measurements across runs per feature, 'remove' will remove those featuares. It can affected for unequal variance analysis.
## remove_m_sequencepeptides : remove the peptides including 'M' sequence
## experiment : "DDA" or "SILAC"

#' @export
MaxQtoMSstatsFormat <- function(evidence, 
                                annotation,
                                proteinGroups, 
                                proteinID="Proteins", 
                                useUniquePeptide=TRUE, 
                                summaryforMultipleRows=max, 
                                fewMeasurements="remove", 
                                removeMpeptides=FALSE,
                                removeOxidationMpeptides=FALSE,
                                removeProtein_with1Peptide=FALSE){
	
    if (is.null(fewMeasurements)) {
        stop('** Please select \'remove\' or \'keep\' for \'fewMeasurements\'.')
    }
    
    if (!is.element(fewMeasurements, c('remove', 'keep'))) {
        stop('** Please select \'remove\' or \'keep\' for \'fewMeasurements\'.')
    }
    
    if (!is.element(proteinID, c('Proteins', 'Leading.razor.protein'))) {
        stop('** Please select \'Proteins\' or \'Leading.razor.proteins\' for \'proteinID\'.')
    }
    
	experiment <- "DDA"
	
	## evidence.txt file
	infile <- evidence
	
	## annotation.txt : Raw.file, Condition, BioReplicate, (IsotopeLabelType)
	annot <- annotation
	
	## check annotation
	required.annotation <- c('Raw.file', 'Condition', 'BioReplicate', 'IsotopeLabelType')
	
	if (!all(required.annotation %in% colnames(annot))) {
	    
	    missedAnnotation <- which(!(required.annotation %in% colnames(annot)))
	    
	    stop(paste("**", toString(required.annotation[missedAnnotation]), 
	               "is not provided in Annotation. Please check the annotation file."))
	}
	
	## check annotation information
	## get annotation
    annotinfo <- unique(annot[, c("Raw.file", "Condition", 'BioReplicate')])	
	
	## Each Run should has unique information about condition and bioreplicate
	check.annot <- xtabs(~Raw.file, annotinfo)
	if ( any(check.annot > 1) ) {
	    stop('** Please check annotation. Each MS run (Raw.file) can\'t have multiple conditions or BioReplicates.' )
	}
	
	################################################
	## 1.1 remove contaminant, reverse proteinID 
	## Contaminant, Reverse column in evidence
	if (is.element("Contaminant", colnames(infile)) & 
	    is.element("+", unique(infile$Contaminant))) {
		infile <- infile[-which(infile$Contaminant %in% "+"), ]
	}
	
	if (is.element("Potential.contaminant", colnames(infile)) & 
	    is.element("+", unique(infile$Potential.contaminant))) {
	    infile <- infile[-which(infile$Potential.contaminant %in% "+"), ]
	}
	
	if (is.element("Reverse", colnames(infile)) & 
	    is.element("+", unique(infile$Reverse))) {
		infile <- infile[-which(infile$Reverse %in% "+"), ]
	}

	## ? Only.identified.by.site column in proteinGroupID? : sometimes, it is not in evidence.txt
	if (is.element("Only.identified.by.site", colnames(infile)) &
	    is.element("+", unique(infile$Only.identified.by.site))) {
		infile <- infile[-which(infile$Only.identified.by.site %in% "+"), ]
	}
	
	message('** + Contaminant, + Reverse, + Only.identified.by.site, proteins are removed.')
	
	
	################################################
	## 1.1.2 matching proteinGroupID protein list 

	## need to check proteinGroupID in evidence and proteinGroup.txt the same
	## 'id' in proteinGroups.txt vs 'Protein.group.IDs' in infile 
	## possible to have some combination in Protein.group.IDs in infile, such as 64;1274;1155;1273 instead of 64, 1274.. separately. 
	## combination of some ids seems not to be used for intensity
	## 2015/02/03
	
	## first, remove contaminants
	if (is.element("Contaminant", colnames(proteinGroups)) & 
	    is.element("+",unique(proteinGroups$Contaminant))) {
		proteinGroups <- proteinGroups[-which(proteinGroups$Contaminant %in% "+"), ]
	}
	
	if (is.element("Potential.contaminant", colnames(proteinGroups)) & 
	    is.element("+",unique(proteinGroups$Potential.contaminant))) {
		proteinGroups <- proteinGroups[-which(proteinGroups$Potential.contaminant %in% "+"), ]
	}

	if (is.element("Reverse", colnames(proteinGroups)) & 
	    is.element("+",unique(proteinGroups$Reverse))) {
		proteinGroups <- proteinGroups[-which(proteinGroups$Reverse %in% "+"), ]
	}

	## ? Only.identified.by.site column in proteinGroupID? : sometimes, it is not in evidence.txt
	if (is.element("Only.identified.by.site", colnames(proteinGroups)) & 
	    is.element("+", unique(proteinGroups$Only.identified.by.site))) {
		proteinGroups <- proteinGroups[-which(proteinGroups$Only.identified.by.site %in% "+"), ]
	}
	
	tempprotein <- proteinGroups
	
	## then take proteins which are included
	infile <- infile[which(infile$Protein.group.IDs %in% unique(tempprotein$id)), ]

	## then use 'protein.IDs' in proteinGroups.txt
	## because if two 'proteins' in evidence.txt are used in one protein ID, need to use certain protein name in infile.
	## for example, protein.IDs in proteinGroups.txt are P05204;O00479. but, two 'proteins in evidence.txt, such as P05204;O00479, and P05204.
	
	tempname <- unique(tempprotein[,c("Protein.IDs", "id")])
	colnames(tempname) <- c("uniqueProteins", "Protein.group.IDs")
	
	infile <- merge(infile, tempname, by="Protein.group.IDs")
	
	## get useful information
	## ? can remove Retention.time column later
	if (experiment == "SILAC") {
		
		infile <- infile[c("uniqueProteins", "Protein.group.IDs", "Sequence", 
		                   "Modified.sequence", "Charge", "Raw.file", 
		                   "Intensity.L", "Intensity.H", "Retention.time", "id")]
		infile.l <- infile[, !(colnames(infile) %in% "Intensity.H")]
		infile.h <- infile[, !(colnames(infile) %in% "Intensity.L")]

		colnames(infile.l)[colnames(infile.l) == "Intensity.L"] <- "Intensity"
		colnames(infile.h)[colnames(infile.h) == "Intensity.H"] <- "Intensity"
		
		## new IsotopeLabelType column
		infile.l$IsotopeLabelType <- "L"
		infile.h$IsotopeLabelType <- "H"

		infile <- rbind(infile.l, infile.h)
		
		rm(infile.l)
		rm(infile.h)
		
	} else {
		
	    get.column <- c("Protein.group.IDs", 
	                    "Sequence", "Modified.sequence", "Modifications", "Charge", 
	                    "Raw.file", "Intensity", "Retention.time", "id")
	    
	    if (proteinID == 'Proteins') {
	        get.column <- c(get.column, 'uniqueProteins')
	    } else {
	        get.column <- c(get.column, 'Leading.razor.protein')
	    }
	    
		infile <- infile[, get.column]
	}
	
	if (proteinID == 'Proteins') {
	    colnames(infile)[colnames(infile) == "uniqueProteins"] <- "Proteins"
	} else {
	    colnames(infile)[colnames(infile) == "Leading.razor.protein"] <- "Proteins"
	}
	
	## remove "_" at the beginning and end
	infile$Modified.sequence <- gsub("_", "", infile$Modified.sequence)
	
	################################################
	## 1.2 remove the peptides including M sequence
	if (removeMpeptides) {
		remove_m_sequence <- unique(infile[grep("M", infile$Modified.sequence), "Modified.sequence"])
 
 		if (length(remove_m_sequence) > 0) {
 			infile <- infile[-which(infile$Modified.sequence %in% remove_m_sequence), ]
 		}
		
		message('** Peptides including M in the sequence are removed.')
		
	}

    if (removeOxidationMpeptides) {
        remove_oxim_sequence <- unique(infile[grep("Oxidation", infile$Modifications), "Modifications"])
    
        if (length(remove_oxim_sequence) > 0) {
            infile <- infile[-which(infile$Modifications %in% remove_oxim_sequence), ]
        }
    
        message('** Peptides including oxidation (M) in the sequence are removed.')
    
    }

    ################################################
	## 2. remove peptides which are used in more than one protein
	## we assume to use unique peptide
	################################################
	if (useUniquePeptide) {
		
		pepcount <- unique(infile[, c("Proteins","Modified.sequence")]) ## Protein.group.IDs or Sequence
		pepcount$Modified.sequence <- factor(pepcount$Modified.sequence)
		
		## count how many proteins are assigned for each peptide
		structure <- aggregate(Proteins ~ ., data=pepcount, length)
		remove_peptide <- structure[structure$Proteins != 1, ]
		
		## remove the peptides which are used in more than one protein
		if (length(remove_peptide$Proteins != 1) != 0) {
			infile <- infile[-which(infile$Modified.sequence %in% remove_peptide$Modified.sequence), ]
			message('** Peptides, that are used in more than one proteins, are removed.')
		}
	}
	
	################################################
	## 3. duplicated(multiple) rows for certain feature and certain runs
	## 	3.1) take highest intensity
	##  3.2) take sum of intensities
	################################################

	## Let's find duplicates
	## first remove NA intensity
	infile <- infile[!is.na(infile$Intensity), ]

	#########################
	### 2.1) general Label-free : one measurement for a feature and a run
	if (experiment == "DDA") {
		## count the number of intensities for feature by runs
		##infile$Feature <- paste(infile$Modified.sequence, infile$Charge, sep="_")
		##structure <- dcast(Feature ~ Raw.file, data=infile, value.var='Intensity')
		##flagduplicate = sum(structure>1)>0	
	
		## take the highest intensity among duplicated or sum of intensities 
		## summaryforMultipleRows="max" or "sum
		infile_w <- .cast_maxquant_to_wide_glf(infile, aggregateFun=summaryforMultipleRows)
	
		## *** remove features which has less than 2 measurements across runs
		## !!! for MSstats v3, we don't need to remove them.
		## good to remove before reformatting to long-format
	
		if (fewMeasurements == "remove") {
			infile_w <- .remove_feature_with_few(infile_w)
			
			message('** Peptide and charge, that have 1 or 2 measurements across runs, are removed.')
		}
	
		## then, go back to long-format
		## good to fill rows with NAs, then now can have balanced data-structure.
		infile_l <- .melt_maxquant_to_long_glf(infile_w)
		
		## need to set 'IsotopeLabelType' because SILAC already has it.
		infile_l$IsotopeLabelType  <-  "L"
	}
	
	#########################
	## 2.2) label-free : however, several runs for a sample.
	
	
	#########################
	## 2.3) SILAC : two measurements for a feature and a run -> one measurements for a feature and a run and condition
	
	if (experiment == "SILAC") {
		## count the number of intensities for feature by runs
		##infile$Feature <- paste(infile$Modified.sequence, infile$Charge, sep="_")
		##structure <- dcast(Feature ~ Raw.file, data=infile, value.var='Intensity')
		##flagduplicate = sum(structure>1)>0	
	
		## take the highest intensity among duplicated or sum of intensities 
		## summaryforMultipleRows="max" or "sum
		infile_w <- .cast_maxquant_to_wide_silac(infile, aggregateFun=summaryforMultipleRows)
	
		## *** remove features which has less than 2 measurements across runs
		## good to remove before reformatting to long-format
	
		if (fewMeasurements=="remove") {
			
			## it is the same across experiments. # measurement per feature. 
			infile_w <- .remove_feature_with_few(infile_w)
		}
	
		## then, go back to long-format
		# good to fill rows with NAs, then now can have balanced data-structure.
		infile_l <- .melt_maxquant_to_long_silac(infile_w)
	}

	################################################
	## 4. remove proteins with only one peptide and charge per protein
	################################################
	
	if (removeProtein_with1Peptide) {
	    ## remove protein which has only one peptide
	    infile_l$feature <- paste(infile_l$Modified.sequence, infile_l$Charge, sep="_")
	  
	    tmp <- unique(infile_l[, c("Proteins", 'feature')])
	    tmp$Proteins <- factor(tmp$Proteins)
	    count <- xtabs( ~ Proteins, data=tmp)
        lengthtotalprotein <- length(count)
    
	    removepro <- names(count[count <= 1])
	  
	    infile_l <- infile_l[-which(infile_l$Proteins %in% removepro), ]
	  
	    message(paste0("** ", length(removepro), 
	                   ' proteins, which have only peptide and charge in a protein, are removed among ', 
	                   lengthtotalprotein, ' proteins.'))
	}
	
	################################################	
	## merge all information
	
	colnames(infile_l)[1] <- "ProteinName"
	colnames(infile_l)[2] <- "PeptideSequence"
	colnames(infile_l)[3] <- "PrecursorCharge"

	## Add in columns for FramentIon & ProductCharge (all values are NA)
	## Add column for IsotopeLabelType (all "L")
	infile_l$FragmentIon <- NA
	infile_l$ProductCharge <- NA

	## Create Condition & Bioreplicate columns; TODO: fill in with correct values
	infile_l <- merge(infile_l, annot, by=c("Raw.file", "IsotopeLabelType"))

	infile_l.final <- infile_l[, c(c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                             "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                             "Condition", "BioReplicate", "Raw.file", "Intensity"))]
	colnames(infile_l.final)[9] <- "Run"

	if (any(is.element(colnames(infile_l), 'Fraction'))) {
	    infile_l.final <- data.frame(infile_l.final,
	                                 "Fraction" = infile_l$Fraction)
	}
	
	
	if (any(is.element(colnames(infile_l), 'TechReplicate'))) {
	    infile_l.final <- data.frame(infile_l.final,
	                                 "TechReplicate" = infile_l$TechReplicate)
	}
	
	infile_l <- infile_l.final
	rm(infile_l.final)
	
	infile_l$PeptideSequence <- factor(infile_l$PeptideSequence)
	infile_l$ProteinName <- factor(infile_l$ProteinName)

	return(infile_l)
}


.cast_maxquant_to_wide_glf <- function(d_long, aggregateFun=aggregateFun){
	data_w <- dcast( Proteins + Modified.sequence + Charge ~ Raw.file, data=d_long, 
                   value.var='Intensity', 
                   fun.aggregate=aggregateFun, na.rm=T,
                   keep=TRUE) 
	## keep=TRUE : will keep the data.frame value as 1 even though there is no values for certain feature and certain run.
	## when there is completely missing in certain feature and certain run, '1' will be filled. Therefore put NA instead of 1.
	data_w[data_w == 1] <- NA
	return(data_w)
}

.cast_maxquant_to_wide_silac <- function(d_long, aggregateFun=aggregateFun){
	
	## check any cell has more than 1
	##data_w = dcast( Proteins + Modified.sequence + Charge + IsotopeLabelType ~ Raw.file, data=d_long, value.var='Intensity') 
	##temp <- data_w[,c(5:ncol(data_w))]
	##head(temp)
	##sum(temp>1)
	##which(temp>1, arr.ind=TRUE)
	##data_w[16300,]
	##d_long[d_long$Modified.sequence=="HIILVLSGK" & d_long$Charge=="2" & d_long$IsotopeLabelType=="L",]
	
	data_w <- dcast( Proteins + Modified.sequence + Charge + IsotopeLabelType ~ Raw.file, data=d_long, 
	                 value.var='Intensity', 
	                 fun.aggregate=aggregateFun, na.rm=T, 
	                 keep=TRUE) 
	## keep=TRUE : will keep the data.frame value as 1 even though there is no values for certain feature and certain run.
  
	## when there is completely missing in certain feature and certain run, '1' will be filled. Therefore put NA instead of 1.
	data_w[data_w == 1] <- NA
	return(data_w)
}

.melt_maxquant_to_long_glf <- function(d_wide){
	data_l <- melt(d_wide, id.vars=c('Proteins', 'Modified.sequence', 'Charge'))
	colnames(data_l)[colnames(data_l) %in% c("variable", "value")] <- c('Raw.file', 'Intensity')
	return(data_l)
}


.melt_maxquant_to_long_silac <- function(d_wide){
	data_l <- melt(d_wide, id.vars=c('Proteins', 'Modified.sequence', 'Charge', "IsotopeLabelType"))
	colnames(data_l)[colnames(data_l) %in% c("variable", "value")] <- c('Raw.file', 'Intensity')
	return(data_l)
}


.remove_feature_with_few <- function(x){
	count_measure <- apply (x[, !(colnames(x) %in% c("Proteins", "Modified.sequence", "Charge"))], 1, 
	                        function ( x ) length ( x[!is.na(x)] ) ) 
	remove_feature_name <- x[count_measure < 3, c("Proteins", "Modified.sequence", "Charge")]

	x$Feature <- paste(x$Proteins, x$Modified.sequence, x$Charge, sep="_")
	remove_feature_name$Feature <- paste(remove_feature_name$Proteins, 
	                                     remove_feature_name$Modified.sequence, 
	                                     remove_feature_name$Charge, sep="_")

	x <- x[-which(x$Feature %in% remove_feature_name$Feature), ]
	x <- x[, -ncol(x)]
	
	return(x)
}

