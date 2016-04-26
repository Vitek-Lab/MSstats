
#############################################
## dataProcess
#############################################

dataProcess  <-  function(raw,
						logTrans=2, 
						normalization="equalizeMedians",
						nameStandards=NULL, 
						betweenRunInterferenceScore=FALSE,
						address="", 
						fillIncompleteRows=TRUE, 
						featureSubset="all",
						remove_proteins_with_interference=FALSE,
            n_top_feature=3,
						summaryMethod="TMP",
						equalFeatureVar=TRUE, 
						filterLogOfSum=TRUE,
						censoredInt="NA",
						cutoffCensored="minFeature",
						MBimpute=TRUE,
						original_scale=FALSE,
						logsum=FALSE,
						remove50missing=FALSE,
						skylineReport=FALSE) {
  
	## save process output in each step
    allfiles <- list.files()
  
	num <- 0
	filenaming <- "msstats"
	finalfile <- "msstats.log"
  
	while(is.element(finalfile,allfiles)) {
		num <- num+1
		finalfile <- paste(paste(filenaming,num,sep="-"),".log",sep="")
	}
  
	session <- sessionInfo()
	sink("sessionInfo.txt")
	print(session)
	sink()
  
	processout <- as.matrix(read.table("sessionInfo.txt", header=T, sep="\t"))
	write.table(processout, file=finalfile, row.names=FALSE)
  
	processout <- rbind(processout, as.matrix(c(" "," ","MSstats - dataProcess function"," "),ncol=1))
  
	## make case-insensitive for function options
  	## ------------------------------------------
  
	normalization <- toupper(normalization)

	## Check correct option or input
	## check right column in input
  
	requiredinput <- c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                     "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                     "Condition", "BioReplicate", "Run", "Intensity")

	## [THT: disambiguation for PeptideSequence & PeptideModifiedSequence - begin]
	## PeptideModifiedSequence is also allowed.
	requiredInputUpper <- toupper(requiredinput)
	providedInputUpper <- toupper(colnames(raw))
	if (all(requiredInputUpper %in% providedInputUpper)) {
    	processout <- rbind(processout, c("The required input : provided - okay"))
    	write.table(processout, file = finalfile, row.names = FALSE)
	} else if (all(setdiff(requiredInputUpper, "PEPTIDESEQUENCE") %in% providedInputUpper) 
             && "PEPTIDEMODIFIEDSEQUENCE" %in% providedInputUpper) {
		processout <- rbind(processout, c("The required input : provided - okay"))
    	write.table(processout, file = finalfile, row.names = FALSE)
    	# if PeptideModifiedSequence is provided instead of PeptideSequence, 
    	# change the column name as PeptideSequence
    	colnames(raw)[which(providedInputUpper == "PEPTIDEMODIFIEDSEQUENCE")]  <-  "PeptideSequence"
  	} else {
		missedInput <- which(!(requiredInputUpper %in% providedInputUpper))
		processout <- rbind(processout, c(paste("ERROR : The required input : ", 
                            paste(requiredinput[missedInput], collapse = ", "), 
                            " are not provided in input - stop")))
		write.table(processout, file = finalfile, row.names = FALSE)
    	stop("Please check the required input. The required input needs (ProteinName, PeptideSequence (or PeptideModifiedSequence), PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity)")
	}
	## [THT: disambiguation for PeptideSequence & PeptideModifiedSequence - end]
  
	## check whether class of intensity is factor or chaterer, if yes, neec to chage as numeric
	if (is.factor(raw$Intensity) | is.character(raw$Intensity)) {	
		suppressWarnings(raw$Intensity <- as.numeric(as.character(raw$Intensity)))
	}
  
	## check whether the intensity has 0 value or negative value
  	if (length(which(raw$Intensity<=0))>0 & !skylineReport) {
  	
		if (is.null(censoredInt)) {
			processout <- rbind(processout,c("ERROR : There are some intensities which are zero or negative values. need to change them. - stop"))
			write.table(processout, file=finalfile,row.names=FALSE)
    
			stop("Intensity has 0 or negative values. Please check these intensities and change them. \n")
  		
		} else if (censoredInt=="NA") {
			
			processout <- rbind(processout,c("ERROR : There are some intensities which are zero or negative values. need to change them. - stop"))
			write.table(processout, file=finalfile,row.names=FALSE)
    
			stop("Intensity has 0 or negative values. Please check these intensities and change them. \n")
  			
  		}
	}
  
	## however skyline report, keep zero and replace with 1,then relace with NA for truncated
	if (skylineReport) {
  	
  		## if there is 'Truncated' column,
  		if (is.element(toupper("Truncated"), toupper(colnames(raw)))) {
  			## remove truncated peaks
  			raw[!is.na(raw$Truncated) & raw$Truncated==TRUE,"Intensity"] <- NA
  	
    		processout <- rbind(processout,c("There are some truncated peaks. They replaced with NA."))
    		write.table(processout, file=finalfile,row.names=FALSE)
    	}
	}
  
	## check unexpected token(":") : retired 2015.04.21, from v3.3.11
	#if (length(grep(":",raw$PeptideSequence))!=0) {
  #  	processout <- rbind(processout,c("ERROR : colon(:) is invalid in PeptideSequence column. Please replace with other entry. - stop"))
  #  	write.table(processout, file=finalfile,row.names=FALSE)
    
  #  	stop("Colon(:) is invalid in PeptideSequence column. Please replace with other entry. \n")
	#}
  
	#if (length(grep(":",raw$FragmentIon))!=0) {
  #  	processout <- rbind(processout,c("ERROR : colon(:) is invalid in FragmentIon column. Please replace with other entry. - stop"))
  #  	write.table(processout, file=finalfile,row.names=FALSE)
    
  #  	stop("Colon(:) is invalid in FragmentIon column. Please replace with other entry. \n")
	#}
  
	## check logTrans is 2,10 or not
	if (logTrans!=2 & logTrans!=10) {
    	processout <- rbind(processout,c("ERROR : Logarithm transformation : log2 or log10 only - stop"))
    	write.table(processout, file=finalfile,row.names=FALSE)
    
    	stop("Only log2 or log10 are posssible.\n")
	}
  
	## check no row for some feature : balanced structure or not  
	if (!(fillIncompleteRows==TRUE | fillIncompleteRows==FALSE) | !is.logical(fillIncompleteRows)) {
		processout <- rbind(processout,c(paste("The required input - fillIncompleteRows : 'fillIncompleteRows' value is wrong. It should be either TRUE or FALSE. - stop")))
		write.table(processout, file=finalfile, row.names=FALSE)
    
		stop("'fillIncompleteRows' must be one of TRUE or FALSE as a logical value.")
	}
  
	## check input for summaryMethod
  
	if (sum(summaryMethod == c("linear", "TMP", "logOfSum")) == 0) {
    	processout <- rbind(processout,c("The required input - summaryMethod : 'summaryMethod' value is wrong. It should be one of 'linear','TMP','logOfSum'. - stop"))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop("'summaryMethod' value is wrong. It should be one of 'linear','TMP','logOfSum'.")
	} else {
		processout <- rbind(processout, c(paste("summaryMethod : ", as.character(summaryMethod), sep="")))
		write.table(processout, file=finalfile, row.names=FALSE)
	}
  
	## check input for cutoffCensored
	if (sum(cutoffCensored==c("minFeature","minRun","minFeatureNRun"))==0) {
		processout <- rbind(processout,c("The required input - cutoffCensored : 'cutoffCensored' value is wrong. It should be one of 'minFeature','minRun','minFeatureNRun'. - stop"))
		write.table(processout, file=finalfile, row.names=FALSE)
    
		stop("'cutoffCensored' value is wrong. It should be one of 'minFeature','minRun','minFeatureNRun'.")
	} else {
		processout <- rbind(processout,c(paste("cutoffCensored : ",as.character(cutoffCensored), sep="")))
		write.table(processout, file=finalfile, row.names=FALSE)
	}
  
	## check input for censoredInt
	if (sum(censoredInt==c("0","NA"))==0 & !is.null(censoredInt)) {
    	processout <- rbind(processout,c("The required input - censoredInt : 'censoredInt' value is wrong. It should be one of '0','NA', NULL. - stop"))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop("'censoredInt' value is wrong. It should be one of '0','NA',NULL.")
  	} else {
  		processout <- rbind(processout,c(paste("censoredInt : ",as.character(censoredInt), sep="")))
    	write.table(processout, file=finalfile, row.names=FALSE)
	}
	
	## [THT: if (!all(normalization %in% c("NONE", "FALSE", "EQUALIZEMEDIANS", "QUANTILE", "GLOBALSTANDARDS")))]
	## [THT: send a warning message if the user mixes "NONE" with any of the last three choices]
	if (!(normalization=="NONE" | normalization=="FALSE" | normalization=="EQUALIZEMEDIANS" | normalization=="QUANTILE" | normalization=="GLOBALSTANDARDS")) {
		processout <- rbind(processout,c(paste("The required input - normalization : 'normalization' value is wrong. - stop")))
		write.table(processout, file=finalfile, row.names=FALSE)
    
		stop("'normalization' must be one of \"None\", \"FALSE\", \"equalizeMedians\", \"quantile\", or \"globalStandards\". Please assign 'normalization' again.")
	} 

	## need the names of global standards
	if (!is.element("NONE",normalization) & !is.element("FALSE",normalization) & is.element("GLOBALSTANDARDS",normalization) & is.null(nameStandards)) {
    
		processout <- rbind(processout,c("ERROR : For normalization with global standards, the names of global standards are needed. Please add 'nameStandards' input."))
		write.table(processout, file=finalfile,row.names=FALSE)
    
		stop ("For normalization with global standards, the names of global standards are needed. Please add 'nameStandards' input." )
	}
  

	## here, need to get standard protein name
	## column name : standardtype..
	## what value it has, normzalition, unique(proteinname)
	## if normalition== "standard" & no normalizaion selection, error message

	## For Skyline
	## required cols : ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, and Condition, BioReplicate, Run, Intensity
  
	## make letters case-insensitive
	colnames(raw) <- toupper(colnames(raw))
	raw.temp <- raw[,c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", "FRAGMENTION", "PRODUCTCHARGE", "ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN", "INTENSITY")]
  
	## before remove, get PeptideSequence and combination of PeptideSequence and precursorcharge for global standard normalization
	tempPeptide <- unique(raw[, c("PEPTIDESEQUENCE", "PRECURSORCHARGE")])
	tempPeptide$PEPTIDE <- paste(tempPeptide$PEPTIDESEQUENCE, tempPeptide$PRECURSORCHARGE, sep="_")
  
	rm(raw)
  
	## assign peptide, transition
	raw.temp <- data.frame(raw.temp,PEPTIDE=paste(raw.temp$PEPTIDESEQUENCE,raw.temp$PRECURSORCHARGE,sep="_"), TRANSITION=paste(raw.temp$FRAGMENTION, raw.temp$PRODUCTCHARGE,sep="_"))
  
	if (length(unique(raw.temp$ISOTOPELABELTYPE))>2) {
    	processout <- rbind(processout,c("ERROR : There are more than two levels of labeling. So far, only label-free or reference-labeled experiment are supported. - stop"))
    	write.table(processout, file=finalfile,row.names=FALSE)
    
    	stop("Statistical tools in MSstats are only proper for label-free or with reference peptide experiments.")
	}
  
	## change light, heavy -> L,H
	## [THT: should check if users really provide light/heavy, L/H, l/h, or something else ]
	## [THT: should also check if users provide only H (instead of L)]
	raw.temp$ISOTOPELABELTYPE <- factor(raw.temp$ISOTOPELABELTYPE)
	if (nlevels(raw.temp$ISOTOPELABELTYPE)==2) {
    	levels(raw.temp$ISOTOPELABELTYPE) <- c("H","L")
  	}
  	if (nlevels(raw.temp$ISOTOPELABELTYPE)==1) {
		levels(raw.temp$ISOTOPELABELTYPE) <- c("L")
	}
  
	raw.temp <- raw.temp[, c("PROTEINNAME", "PEPTIDE", "TRANSITION", "ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE","RUN", "INTENSITY")]
  
	colnames(raw.temp) <- c("Protein", "Peptide", "Transition", "Label", "Condition", "Sample", "Run", "Intensity")
  
  
	## create work data for quant analysis
	## -----------------------------------
	
	work <- data.frame(PROTEIN=raw.temp[,"Protein"], PEPTIDE=raw.temp[,"Peptide"], TRANSITION=raw.temp[,"Transition"], FEATURE=paste(raw.temp[,"Peptide"], raw.temp[,"Transition"],sep="_"), LABEL=raw.temp[,"Label"], GROUP_ORIGINAL=raw.temp[,"Condition"], SUBJECT_ORIGINAL=raw.temp[,"Sample"], RUN=raw.temp[,"Run"], GROUP=0,SUBJECT=0)
  
	work$GROUP_ORIGINAL <- factor(work$GROUP_ORIGINAL)
	work$SUBJECT_ORIGINAL <- factor(work$SUBJECT_ORIGINAL, levels=unique(work$SUBJECT_ORIGINAL))
	work$LABEL <- factor(work$LABEL, levels=levels(work$LABEL))

	work[work$LABEL=="L", "GROUP"] <- work[work$LABEL=="L", "GROUP_ORIGINAL"]
	work[work$LABEL=="L", "SUBJECT"] <- work[work$LABEL=="L", "SUBJECT_ORIGINAL"]
  
	work <- data.frame(work, SUBJECT_NESTED=paste(work[, "GROUP"], work[, "SUBJECT"], sep="."))
  
	processout <- rbind(processout, c("New input format : made new columns for analysis - okay"))
	write.table(processout, file=finalfile, row.names=FALSE)
  
	work <- data.frame(work, INTENSITY=raw.temp$Intensity)
  
  
	## log transformation
	work$ABUNDANCE <- work$INTENSITY
  
	## now, INTENSITY keeps original values.
    
	### [THT: should check intensity between 0-1 as well...]
	## change zero with 1 for log transformation
	## NA means no observation. assume that spectral tools are not report if no observation. zero means detected but zero. 
	work[!is.na(work$ABUNDANCE) & work$ABUNDANCE==0,"ABUNDANCE"] <- 1
  	
	processout <- rbind(processout,c("There are some intensities which are zero. Intensities with zero are replaced with 1 in order to do log transformation."))
	write.table(processout, file=finalfile,row.names=FALSE)
    
	## based on logTrans option, assign log transformation
	## remove log2 or log10 intensity
	### [THT: add one more conidtion to have the program complain if a user 
	### provide unexpected value for logTrans]
	if (logTrans == 2) {
    	work$ABUNDANCE <- log2(work$ABUNDANCE)
	} else if (logTrans == 10) {
    	work$ABUNDANCE <- log10(work$ABUNDANCE)
	} 	
  
	processout <- rbind(processout,c(paste("Logarithm transformation: log",logTrans," transformation is done - okay", sep="")))
	write.table(processout, file=finalfile,row.names=FALSE)
   
	## Check multi-method or not : multiple run for a replicate
  
	##standardFeature <- unique(work[work$RUN=="1","FEATURE"]) ## if some feature are missing for this spedific run, it could be error. that is why we need balanced design.
  
	##countdiff = tapply (work$FEATURE, work$RUN, function ( x ) length(setdiff(unique(x),standardFeature)) ) 
  
	## whether multirun or not : we assume that different method has completely different feature
	work$RUN <- factor(work$RUN)
	multirun <- .countMultiRun(work)
	checkMultirun <- any(multirun==0)
  

	## if multirun, make new column 'method'
    ## -------------------------------------
	if (checkMultirun) { ## when checkMultirun is TRUE, it means there are more than 1 method.
    
		## make column 'method'
    	work$METHOD <- 1
    	numincreasing=1
    
    	## get run which has different feature names from run1
    	while(length(multirun[multirun==0])!=0) { ## until there is no more unique feature per run
			nextmethod <- names(multirun[multirun == 0])
			numincreasing <- numincreasing+1
			work[which(work$RUN %in% nextmethod), "METHOD"] <- numincreasing
      
			worktemp <- work[which(work$RUN %in% nextmethod),]
			worktemp$RUN <- factor(worktemp$RUN)
      
			multirun <- .countMultiRun(worktemp)	
    	}
    
    	processout <- rbind(processout, c(paste("Multiple methods are existed : ", length(unique(work$METHOD)), "methods per MS replicate.")))
    	write.table(processout, file=finalfile,row.names=FALSE)
    
  	}else{
    	work$METHOD <- 1
  	}
  
	## check messingness for multirun 
  
	## check no value for some feature : balanced structure or not
	## need to separate label-free or label-based
  
	processout <- rbind(processout,c(paste("fillIncompleteRows = ",fillIncompleteRows,sep="")))
	write.table(processout, file=finalfile, row.names=FALSE)
  
	## [THT: better to write a function for single method, and call that function
	## here and for the case with multuple methods]
	## only 1 method
  
	if (!checkMultirun) {
    
		## label-free experiments
    	if (nlevels(work$LABEL) == 1) {
      
      		## get feature by Run count of data
      		structure = tapply ( work$ABUNDANCE, list ( work$FEATURE, work$RUN ) , function ( x ) length ( x ) ) 
      
      		## structure value should be 1 for label-free, if not there are missingness. if more there are duplicates.
      
      		flagmissing = sum(is.na(structure))>0
      		flagduplicate = sum(structure[!is.na(structure)]>1)>0
      
      		### if there is missing rows
      		if ( flagmissing ) {
        		processout <- rbind(processout,c("CAUTION: the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with 'NA'. The incomplete rows are listed below."))
        		write.table(processout, file=finalfile,row.names=FALSE)
        
        		message("CAUTION : the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with 'NA'. The incomplete rows are listed below.")
        
        		## first, which run has missing	
        		runstructure <- apply ( structure, 2, function ( x ) sum ( is.na ( x ) ) ) > 0
        
        		## get the name of Run
        		runID <- names(runstructure[runstructure==TRUE])
        
        		## for missign row, need to assign before looping
       			missingwork <- NULL
        
        		## then for each run, which features are missing,
				for(j in 1:length(runID)) {
          
					## get subject, group information for this run
					nameID <- unique(work[work$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
          
					## get feature ID
					featureID <- structure[,colnames(structure)==runID[j]]
          
					## get feature ID which has no measuremnt.
					finalfeatureID <- featureID[is.na(featureID)]
          
					## print features ID	 	
					message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
          
					## save in process file.
					processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
					write.table(processout, file=finalfile,row.names=FALSE)
          
					## add missing rows if option is TRUE
					if (fillIncompleteRows) {
            
            			tempTogetfeature <- work[which(work$FEATURE %in% names(finalfeatureID)), ]
            
            			## get PROTEIN and FEATURE infomation
            			tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
            
            			## merge feature info and run info as 'work' format
            			tempmissingwork <- data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
            
            			## merge with tempary space, missingwork
            			missingwork <- rbind(missingwork,tempmissingwork)
          			} # end fillIncompleteRows options
        		} # end loop for run ID
        
        		## [THT: this part can probably be merged into the above. 
        		## Also, it might be better to check fillIncompleteRows earlier
        		## and terminate the process when it's FALSE]
        		if (fillIncompleteRows) {
          
          			## merge with work
          			## in future, use rbindlist?? rbindlist(list(work, missingwork))
          			work <- rbind(work, missingwork)
          
          			## print message
          			message("\n DONE : Incomplete rows for missing peaks are added with intensity values=NA. \n")
          
          			## save in process file.
          			processout <- rbind(processout, "Incomplete rows for missing peaks are added with intensity values=NA. - done, Okay")
          			write.table(processout, file=finalfile, row.names=FALSE)
          
        		} else {
          
          			## save in process file.
          			processout <- rbind(processout,"Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
          			write.table(processout, file=finalfile,row.names=FALSE)
          
          			stop("Please check whether features in the list are generated from spectral processing tool or not. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
          
        		}
      		} # end for flag missing
      		
			## if there are duplicates measurements
			if (flagduplicate) {
        
        		## first, which run has duplicates
        		runstructure <- apply ( structure, 2, function ( x ) sum (x[!is.na(x)] > 1 ) > 0 )
        
        		runID <- names(runstructure[runstructure==TRUE])
        
        		## then for each run, which features have duplicates,
        		for(j in 1:length(runID)) {
          
          			nameID <- unique(work[work$RUN == runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP","SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
          
          			featureID <- structure[, colnames(structure)==runID[j]]
          			finalfeatureID <- featureID[!is.na(featureID) & featureID > 1]
          
          			message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
          
          			## save in process file.
          			processout <- rbind(processout, c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) , ", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some features (", paste(names(featureID[is.na(featureID)]), collapse=", "), ")", sep="" )))
          			write.table(processout, file=finalfile, row.names=FALSE)
        		}
        
        		## save in process file.
        		processout <- rbind(processout,"Please remove duplicate rows in the list above. ")
        		write.table(processout, file=finalfile,row.names=FALSE)
        
        		stop("Please remove duplicate rows in the list above.\n")		
			} # end flag duplicate
      
			## no missing and no duplicates
      		if (!flagmissing & !flagduplicate) {
        		processout <- rbind(processout, c("Balanced data format with NA for missing feature intensities - okay"))
        		write.table(processout, file=finalfile, row.names=FALSE)
      		} 
      
      		## end label-free
    	}else{ 
      
			## label-based experiment
      
      		## count the reference and endobenous separately
      		work.l <- work[work$LABEL=="L", ]
      		work.h <- work[work$LABEL=="H", ]
       
      		## get feature by Run count of data
      		structure.l = tapply(work.l$ABUNDANCE, list(work.l$FEATURE, work.l$RUN), function (x) length (x) ) 
      		structure.h = tapply(work.h$ABUNDANCE, list(work.h$FEATURE, work.h$RUN), function (x) length (x) ) 
           
			## first, check some features which completely missing across run
      		missingcomplete.l <- NULL	
      		missingcomplete.h <- NULL	
      
      		## 1. reference peptides
      		featurestructure.h <- apply(structure.h, 1, function (x) sum(is.na(x)))
      
      		## get feature ID of reference which are completely missing across run
      		featureID.h <- names(featurestructure.h[featurestructure.h == ncol(structure.h)])
      
      		if (length(featureID.h) > 0) {
        		## print message
        		message(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "), ". Please check whether features in the list are correctly generated from spectral processing tool. \n", sep=""))
        
				## save in process file.
				processout <- rbind(processout,c(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "), ". Please check whether features in the list are correctly generated from spectral processing tool.", sep="")))
				write.table(processout, file=finalfile, row.names=FALSE)
        
        		## add missing rows if option is TRUE
        		if (fillIncompleteRows) {
          
          			## get unique Run information
          			nameID <- unique(work.h[, c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
          
          			## get PROTEIN and FEATURE information
           			## here use whole work dataset
          			tempTogetfeature <- work[which(work$FEATURE %in% featureID.h), ]
          			tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
          
          			## then generate data.frame for missingness,
          			#for(j in 1:nrow(nameID)) {
            
           			#	## merge feature info and run info as 'work' format
            		#	tempmissingwork <- data.frame(tempfeatureID, LABEL="H",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
            
            		#	## merge with tempary space, missingwork
            		#	missingcomplete.h <- rbind(missingcomplete.h, tempmissingwork)
          			#}
                
          			# MC : 2016.04.21 : use merge for simplicity
          			tmp <- merge(nameID, tempfeatureID, by=NULL)
          			missingcomplete.h <- data.frame(PROTEIN=tmp$PROTEIN, PEPTIDE=tmp$PEPTIDE, TRANSITION=tmp$TRANSITION, FEATURE=tmp$FEATURE, LABEL="H", GROUP_ORIGINAL=tmp$GROUP_ORIGINAL, SUBJECT_ORIGINAL=tmp$SUBJECT_ORIGINAL, RUN=tmp$RUN, GROUP=tmp$GROUP, SUBJECT=tmp$SUBJECT, SUBJECT_NESTED=tmp$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=tmp$METHOD)
          			rm(tmp)
                
        		}	# end fillIncompleteRows option     
      		} # end for reference peptides
      
      		## 2. endogenous peptides
      		featurestructure.l <- apply(structure.l, 1, function (x) sum(is.na(x)))
      
      		## get feature ID of reference which are completely missing across run
      		featureID.l <- names(featurestructure.l[featurestructure.l==ncol(structure.l)])
      
      		if (length(featureID.l) > 0) {
        		## print message
        		message(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENDOGENOUS features are ", paste(featureID.l, collapse=", "), ". Please check whether features in the list are correctly generated from spectral processing tool. \n", sep=""))
        
        		## save in process file.
        		processout <- rbind(processout,c(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENDOGENOUS features are ", paste(featureID.l, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool. \n", sep="")))
        		write.table(processout, file=finalfile,row.names=FALSE)
        
        		## add missing rows if option is TRUE
        		if (fillIncompleteRows) {
          
          			## get unique Run information
          			nameID <- unique(work.l[, c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
          
          			## get PROTEIN and FEATURE information
          			## here use whole work dataset
          			tempTogetfeature <- work[which(work$FEATURE %in% featureID.l), ]
          			tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
          
          			## then generate data.frame for missingness,
          			#for (j in 1:nrow(nameID)) {
            
            		#	## merge feature info and run info as 'work' format
            		#	tempmissingwork <- data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
            
            		#	## merge with tempary space, missingwork
            		#	missingcomplete.l <- rbind(missingcomplete.l, tempmissingwork)
          			#}
                
                # MC : 2016.04.21 : use merge for simplicity
                tmp <- merge(nameID, tempfeatureID, by=NULL)
          			missingcomplete.l <- data.frame(PROTEIN=tmp$PROTEIN, PEPTIDE=tmp$PEPTIDE, TRANSITION=tmp$TRANSITION, FEATURE=tmp$FEATURE, LABEL="L", GROUP_ORIGINAL=tmp$GROUP_ORIGINAL, SUBJECT_ORIGINAL=tmp$SUBJECT_ORIGINAL, RUN=tmp$RUN, GROUP=tmp$GROUP, SUBJECT=tmp$SUBJECT, SUBJECT_NESTED=tmp$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=tmp$METHOD)
        		    rm(tmp)
            } # end fillIncompleteRows option
      		} # end endogenous peptides
          
			## second, check other some missingness
      
      		## for missign row, need to assign before looping. need to assign at the beginning because it need either cases, with missingness or not
      		missingwork.l <- NULL
      		missingwork.h <- NULL
      
     		## structure value should be 1 for reference and endogenous separately, if not there are missingness. if more there are duplicates.
      
      		## if count of NA is not zero and not number of run (excluding complete missingness across runs)
      
      		missing.l <- names(featurestructure.l[featurestructure.l != ncol(structure.l) & featurestructure.l != 0])
      		missing.h <- names(featurestructure.h[featurestructure.h != ncol(structure.h) & featurestructure.h != 0])
      
			    flagmissing.l = length(missing.l) > 0
			    flagmissing.h = length(missing.h) > 0
      
			    ## structure value is greater than 1, there are duplicates
			    flagduplicate.l = sum(structure.l[!is.na(structure.l)] > 1) > 0
			    flagduplicate.h = sum(structure.h[!is.na(structure.h)] > 1) > 0
      
      		## if there is missing rows for endogenous
      		if ( flagmissing.l | flagmissing.h ) {
        		processout <- rbind(processout,c("CAUTION: the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with 'NA'. The incomplete rows are listed below."))
        		write.table(processout, file=finalfile, row.names=FALSE)
        
       			message("CAUTION : the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with 'NA'. The incomplete rows are listed below.")
        
        		## endogenous intensities
        		if (flagmissing.l) {
          
                if (length(missing.l) > 1){
          			  runstructure <- apply ( structure.l[which(rownames(structure.l) %in% missing.l), ], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
          			} else if (length(missing.l) == 1) {
          			  runstructure <- is.na ( structure.l[which(rownames(structure.l) %in% missing.l), ]) > 0
          			}
          			
          			## get the name of Run
          			runID <- names(runstructure[runstructure==TRUE])
          
          			## then for each run, which features are missing,
          			for(j in 1:length(runID)) {
            
            			## get subject, group information for this run
            			nameID <- unique(work.l[work.l$RUN==runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
            
            			# MC : 2016/04/21. if there is one row, can't catch up data.frame
            			## get feature ID
            			if (length(missing.l) > 1){
            			  featureID <- structure.l[which(rownames(structure.l) %in% missing.l), colnames(structure.l) == runID[j]]
            			  
            			  ## get feature ID which has no measuremnt.
            			  finalfeatureID <- names(featureID[is.na(featureID)])
            			  
            			  ## print features ID	 	
            			  message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(finalfeatureID, collapse=", "),")", sep="" ))
            			  
            			  ## save in process file.
            			  processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(finalfeatureID, collapse=", "),")", sep="" )))
            			  write.table(processout, file=finalfile,row.names=FALSE)
            			  
            			} else if (length(missing.l) == 1) {
            			  
            			  finalfeatureID <- missing.l
                    
            			  ## print features ID   	
            			  message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", finalfeatureID,")", sep="" ))
            			  
            			  ## save in process file.
            			  processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", finalfeatureID,")", sep="" )))
            			  write.table(processout, file=finalfile,row.names=FALSE)
            			  
            			}	
            
            			## add missing rows if option is TRUE
            			if (fillIncompleteRows) {
              
              				tempTogetfeature <- work.l[which(work.l$FEATURE %in% finalfeatureID), ]
              
              				## get PROTEIN and FEATURE infomation
              				tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
              
              				## merge feature info and run info as 'work' format
              				tempmissingwork <- data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
              
              				## merge with tempary space, missingwork
              				missingwork.l <- rbind(missingwork.l,tempmissingwork)
            			} # end fillIncompleteRows options
          			} # end loop for run ID
        		} # end for endogenous
        
        		## reference intensities
        		if (flagmissing.h) {
          
          			## first, which run has missing	
                if (length(missing.h) > 1){
          			    runstructure <- apply ( structure.h[which(rownames(structure.h) %in% missing.h), ], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
                } else if (length(missing.h) == 1) {
                    runstructure <- is.na ( structure.h[which(rownames(structure.h) %in% missing.h), ]) > 0
                }
                
          			## get the name of Run
          			runID <- names(runstructure[runstructure==TRUE])
          
          			## then for each run, which features are missing,
          			for(j in 1:length(runID)) {
            
            			## get subject, group information for this run
            			nameID <- unique(work.h[work.h$RUN==runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
            
            			# MC : 2016/04/21. if there is one row, can't catch up data.frame
            			## get feature ID
            			if (length(missing.h) > 1){
                    featureID <- structure.h[which(rownames(structure.h) %in% missing.h), colnames(structure.h) == runID[j] ]
                    
                    ## get feature ID which has no measuremnt.
                    finalfeatureID <- names(featureID[is.na(featureID)])
                    
                    ## print features ID	 	
                    message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(finalfeatureID, collapse=", "),")", sep="" ))
                    
                    ## save in process file.
                    processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(finalfeatureID, collapse=", "),")", sep="" )))
                    write.table(processout, file=finalfile,row.names=FALSE)
                    
            			} else if (length(missing.h) == 1) {
            			 
            			  finalfeatureID <- missing.h
                    
            			  ## print features ID   	
            			  message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", finalfeatureID,")", sep="" ))
            			  
            			  ## save in process file.
            			  processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", finalfeatureID,")", sep="" )))
            			  write.table(processout, file=finalfile,row.names=FALSE)
                    
            			}
            			           			                         	            
            			## add missing rows if option is TRUE
            			if (fillIncompleteRows) {
              
              				tempTogetfeature <- work.h[which(work.h$FEATURE %in% finalfeatureID), ]
              
              				## get PROTEIN and FEATURE infomation
              				tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
              
              				## merge feature info and run info as 'work' format
              				tempmissingwork <- data.frame(tempfeatureID, LABEL="H",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
              
              				## merge with tempary space, missingwork
              				missingwork.h <- rbind(missingwork.h, tempmissingwork)
                                    				
            			} # end fillIncompleteRows options
          			} # end loop for run ID
        		} # end for endogenous  
      		} # end for flag missing
      
      		## merge missing rows if fillIncompleteRows=TRUE or message.
      		if (fillIncompleteRows) {
        
        		## merge with work
        		## in future, use rbindlist?? rbindlist(list(work, missingwork))
        		work <- rbind(work,missingcomplete.l, missingcomplete.h, missingwork.l, missingwork.h)
        
        		## print message
        		message("\n DONE : Incomplete rows for missing peaks are added with intensity values=NA. \n")
        
        		## save in process file.
        		processout <- rbind(processout,"Incomplete rows for missing peaks are added with intensity values=NA. - done, Okay")
        		write.table(processout, file=finalfile, row.names=FALSE)
        
      		} else if (!is.null(missingcomplete.l) | !is.null(missingcomplete.h) | !is.null(missingwork.l) | !is.null(missingwork.l) ) {
        
        		## save in process file.
       	 		processout <- rbind(processout,"Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        		write.table(processout, file=finalfile,row.names=FALSE)
        
        		stop("Please check whether features in the list are generated from spectral processing tool or not. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        
			}
      
			## if there are duplicates measurements
      		if (flagduplicate.h) {
        
        		## first, which run has duplicates
        		runstructure <- apply ( structure.h, 2, function ( x ) sum ( x[!is.na(x)] > 1 )>0 )
        
        		runID <- names(runstructure[runstructure==TRUE])
        
        		## then for each run, which features have duplicates,
        		for(j in 1:length(runID)) {
          
          			nameID <- unique(work[work$RUN==runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
          
          			featureID <- structure.h[,colnames(structure.h)==runID[j]]
          			finalfeatureID <- featureID[!is.na(featureID) & featureID > 1]
          
          			message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some REFERENCE features (", paste(names(finalfeatureID), collapse=", "), ")", sep="" ))
          
          			## save in process file.
         	 		processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some REFERENCE features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
         		 	write.table(processout, file=finalfile,row.names=FALSE)
        		}
        
        		## save in process file.
        		processout <- rbind(processout,"Please remove duplicate rows in the list above. ")
        		write.table(processout, file=finalfile, row.names=FALSE)
        
        		stop("Please remove duplicate rows in the list above.\n")		
			} # end flag duplicate for reference
      
			if (flagduplicate.l) {
        
        		## first, which run has duplicates
        		runstructure <- apply ( structure.l, 2, function ( x ) sum ( x[!is.na(x)] > 1 )>0 )
        
        		runID <- names(runstructure[runstructure == TRUE])
        
        		## then for each run, which features have duplicates,
        		for (j in 1:length(runID)) {
          
          			nameID <- unique(work[work$RUN==runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
          
          			featureID <- structure.l[,colnames(structure.l)==runID[j]]
          			finalfeatureID <- featureID[!is.na(featureID) & featureID > 1]
          
          			message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some ENDOGENOUS features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
          
          			## save in process file.
          			processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some ENDOGENOUS features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
          			write.table(processout, file=finalfile,row.names=FALSE)
        		}
        
				## save in process file.
				processout <- rbind(processout,"ERROR : Please remove duplicate rows in the list above. ")
				write.table(processout, file=finalfile,row.names=FALSE)
        
				stop("ERROR : Please remove duplicate rows in the list above.\n")		
			} # end flag duplicate for endogenous
      
			## no missing and no duplicates
			if (!flagmissing.h & !flagmissing.l & !flagduplicate.h & !flagduplicate.l) {
        		processout <- rbind(processout, c("Balanced data format with NA for missing feature intensities - okay"))
        		write.table(processout, file=finalfile, row.names=FALSE)
			} 		 
		} # end 1 method
    
	} else { # multiple methods
    
    	allflagmissing <- NULL
		allflagduplicate <- NULL
    
    	## check each method
    	for (k in 1:length(unique(work$METHOD))) {
      
			worktemp <- work[work$METHOD==k,]
			worktemp$RUN <- factor(worktemp$RUN)
			worktemp$FEATURE <- factor(worktemp$FEATURE)
      
			structure = tapply ( worktemp$ABUNDANCE, list ( worktemp$FEATURE, worktemp$RUN ) , function ( x ) length ( x ) ) 
      
			## structure value should be 2 for labeled, 1 for label-free, if not there are missingness
			if (nlevels(worktemp$LABEL)==2) { ## label-based
				flag = sum(is.na(structure)) > 0 | sum(structure[!is.na(structure)] < 2) > 0
			} else {  ## label-free
				flag = sum(is.na(structure)) > 0
			}
      
			allflagmissing <- c(allflagmissing,flag)
      
      		## for duplicate
			if (nlevels(worktemp$LABEL)==2) { # label-based
        		worktemp.h <- worktemp[worktemp$LABEL=="H", ]
        		worktemp.l <- worktemp[worktemp$LABEL=="L", ]
        
        		structure.h = tapply ( worktemp.h$ABUNDANCE, list ( worktemp.h$FEATURE, worktemp.h$RUN ) , function ( x ) length ( x ) ) 
        		structure.l = tapply ( worktemp.l$ABUNDANCE, list ( worktemp.l$FEATURE, worktemp.l$RUN ) , function ( x ) length ( x ) ) 
        
        		flagduplicate <- sum(structure.h[!is.na(structure.h)] > 1) > 0 | sum(structure.l[!is.na(structure.l)] > 1) > 0
        
      		} else {  # label-free
        		flagduplicate <- sum(structure[!is.na(structure)]>1)>0
      		}
      
      		allflagduplicate <- c(allflagduplicate,flag)
      
    	} # end to check any flag among methods
    
    	if ( sum(allflagmissing) != 0 ) {
			processout <- rbind(processout, c("CAUTION: the input dataset has incomplete rows. Missing feature intensities should be present in the dataset, and their intensities should be indicated with 'NA'. The incomplete rows are listed below."))
      		write.table(processout, file=finalfile, row.names=FALSE)
      
			message("CAUTION : the input dataset has incomplete rows. Missing feature intensities should be present in the dataset, and their intensities should be indicated with 'NA'. The incomplete rows are listed below.")
      
			## for missign row, need to assign before looping
			missingwork <- NULL
      
			missingcomplete.h <- NULL
			missingcomplete.l <- NULL
			missingwork.h <- NULL
			missingwork.l <- NULL
      
			for (k in 1:length(unique(work$METHOD))) {
        
        		## see which method has missing rows
        		if (allflagmissing[k]) {
          			worktemp <- work[work$METHOD==k, ]
          			worktemp$RUN <- factor(worktemp$RUN)
          			worktemp$FEATURE <- factor(worktemp$FEATURE)
          
          			if (nlevels(worktemp$LABEL) == 1) { ## label-free
            
            			structure = tapply ( worktemp$ABUNDANCE, list ( worktemp$FEATURE, worktemp$RUN ) , function ( x ) length ( x ) ) 
            
            			## first, which run has missing	
            			runstructure <- apply ( structure, 2, function ( x ) sum ( is.na ( x ) ) ) > 0
            
            			## get the name of Run
            			runID <- names(runstructure[runstructure==TRUE])
            
            			## then for each run, which features are missing,
            			for (j in 1:length(runID)) {
              
              				nameID <- unique(worktemp[worktemp$RUN==runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
              
              				## get feature ID
              				featureID <- structure[, colnames(structure)==runID[j]]
              
              				## get feature ID which has no measuremnt.
              				finalfeatureID <- featureID[is.na(featureID)]
              
              				## print features ID	 	
              				message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
              
              				## save in process file.
              				processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
              				write.table(processout, file=finalfile, row.names=FALSE)
              
              				## add missing rows if option is TRUE
              				if (fillIncompleteRows) {
                
                				tempTogetfeature <- work[which(work$FEATURE %in% names(finalfeatureID)), ]
                
                				## get PROTEIN and FEATURE infomation
                				tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
                
                				## merge feature info and run info as 'work' format
                				tempmissingwork <- data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
                
               			 		## merge with tempary space, missingwork
                				missingwork <- rbind(missingwork, tempmissingwork)
              				} # end fillIncompleteRows options
            			} # end loop for run
            
          			} else { # end label-free
            
            			## label-based
            			## count the reference and endobenous separately
            			work.l <- worktemp[worktemp$LABEL=="L", ]
            			work.h <- worktemp[worktemp$LABEL=="H", ]
            
            			## get feature by Run count of data
            			structure.l <- tapply ( work.l$ABUNDANCE, list(work.l$FEATURE, work.l$RUN), function (x) length (x) ) 
            			structure.h <- tapply ( work.h$ABUNDANCE, list(work.h$FEATURE, work.h$RUN), function (x) length (x) ) 
            
            			## 1. reference peptides
            			featurestructure.h <- apply(structure.h, 1, function (x) sum(is.na(x)))
            
            			## get feature ID of reference which are completely missing across run
            			featureID.h <- names(featurestructure.h[featurestructure.h==ncol(structure.h)])
            
            			if (length(featureID.h) > 0) {
              				## print message
              				message(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool. \n", sep=""))
              
              				## save in process file.
              				processout <- rbind(processout,c(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool.", sep="")))
              				write.table(processout, file=finalfile, row.names=FALSE)
              
              				## add missing rows if option is TRUE
              				if (fillIncompleteRows) {
                
                				## get unique Run information
                				nameID <- unique(work.h[, c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
                
                   				## get PROTEIN and FEATURE information
                				## here use whole worktemp dataset
                				tempTogetfeature <- worktemp[which(worktemp$FEATURE %in% featureID.h), ]
                				tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
                
                				## then generate data.frame for missingness,
                				for (j in 1:nrow(nameID)) {
                  
                  					## merge feature info and run info as 'work' format
                  					tempmissingwork <- data.frame(tempfeatureID, LABEL="H",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
                  
                  					## merge with tempary space, missingwork
                  					missingcomplete.h <- rbind(missingcomplete.h, tempmissingwork)
                				}
              				} # end fillIncompleteRows option
						} # end for reference peptides
            
            			## 2. endogenous peptides
            			featurestructure.l <- apply(structure.l, 1, function (x) sum(is.na(x)))
            
            			## get feature ID of reference which are completely missing across run
            			featureID.l <- names(featurestructure.l[featurestructure.l==ncol(structure.l)])
            
            			if (length(featureID.l) > 0) {
              				## print message
              				message(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENDOGENOUS features are ", paste(featureID.l, collapse=", "), ". Please check whether features in the list are correctly generated from spectral processing tool. \n", sep=""))
              
              				## save in process file.
              				processout <- rbind(processout,c(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENCOGENOUS features are ", paste(featureID, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool. \n", sep="")))
              				write.table(processout, file=finalfile, row.names=FALSE)
              
              				## add missing rows if option is TRUE
              				if (fillIncompleteRows) {
                
                				## get unique Run information
                				nameID <- unique(work.l[, c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
                
                				## get PROTEIN and FEATURE information
                				## here use whole worktemp dataset
                				tempTogetfeature <- worktemp[which(worktemp$FEATURE %in% featureID.l), ]
                				tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
                
                				## then generate data.frame for missingness,
                				for(j in 1:nrow(nameID)) {
                  
                  					## merge feature info and run info as 'work' format
                  					tempmissingwork <- data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
                  
                  					## merge with tempary space, missingwork
                  					missingcomplete.l <- rbind(missingcomplete.l, tempmissingwork)
                				}
              				} # end fillIncompleteRows option
            			} # end endogenous peptides
            
						## second, check other some missingness
            
            			## structure value should be 1 for reference and endogenous separately, if not there are missingness. if more there are duplicates.
            
            			## if count of NA is not zero and not number of run (excluding complete missingness across runs)
            			missing.l <- names(featurestructure.l[featurestructure.l!=ncol(structure.l) & featurestructure.l != 0])
            			missing.h <- names(featurestructure.h[featurestructure.h!=ncol(structure.h) & featurestructure.h != 0])
            
            			flagmissing.l <- length(missing.l) > 0
            			flagmissing.h <- length(missing.h) > 0
            
            			## structure value is greater than 1, there are duplicates
           	 			flagduplicate.l <- sum(structure.l[!is.na(structure.l)] > 1) > 0
            			flagduplicate.h <- sum(structure.h[!is.na(structure.h)] > 1) > 0
            
            			## if there is missing rows for endogenous
            			if (flagmissing.l | flagmissing.h) {
              				processout <- rbind(processout,c("CAUTION: the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with 'NA'. The incomplete rows are listed below."))
              				write.table(processout, file=finalfile, row.names=FALSE)
              
              				message("CAUTION : the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with 'NA'. The incomplete rows are listed below.")
              
             				## endogenous intensities
              				if (flagmissing.l) {
                
                				## first, which run has missing	
                				runstructure <- apply ( structure.l[-which(rownames(structure.l) %in% featureID.l),], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
                
                				## get the name of Run
                				runID <- names(runstructure[runstructure==TRUE])
                
                				## then for each run, which features are missing,
                				for (j in 1:length(runID)) {
                  
                  					## get subject, group information for this run
                  					nameID <- unique(work.l[work.l$RUN==runID[j], c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT", "SUBJECT_NESTED", "RUN", "METHOD")])
                  
                  					## get feature ID
                  					featureID <- structure.l[-which(rownames(structure.l) %in% featureID.l), colnames(structure.l)==runID[j]]
                  
                  					## get feature ID which has no measuremnt.
                  					finalfeatureID <- featureID[is.na(featureID)]
                  
                  					## print features ID	 	
                  					message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[, "GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
                  
                  					## save in process file.
                  					processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
                  					write.table(processout, file=finalfile, row.names=FALSE)
                  
                  					## add missing rows if option is TRUE
                  					if (fillIncompleteRows) {
                    
                    					tempTogetfeature <- work.l[which(work.l$FEATURE %in% names(finalfeatureID)), ]
                    
                    					## get PROTEIN and FEATURE infomation
                    					tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
                    
                    					## merge feature info and run info as 'work' format
                    					tempmissingwork <- data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
                    
                    					## merge with tempary space, missingwork
                    					missingwork.l <- rbind(missingwork.l, tempmissingwork)
                  					} # end fillIncompleteRows options
                				} # end loop for run ID
              				} # end for endogenous
              
              				## reference intensities
              				if (flagmissing.h) {
                
                				## first, which run has missing	
                				runstructure <- apply ( structure.h[-which(rownames(structure.h) %in% featureID.h),], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
                
                				## get the name of Run
                				runID <- names(runstructure[runstructure==TRUE])
                
                				## then for each run, which features are missing,
                				for (j in 1:length(runID)) {
                  
                  					## get subject, group information for this run
                  					nameID <- unique(work.h[work.h$RUN==runID[j], c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
                  
                  					## get feature ID
                  					featureID <- structure.h[-which(rownames(structure.h) %in% featureID.h), colnames(structure.h)==runID[j]]
                  
                  					## get feature ID which has no measuremnt.
                  					finalfeatureID <- featureID[is.na(featureID)]
                  
                  					## print features ID	 	
                  					message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
                  
                 					## save in process file.
                  					processout <- rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
                  					write.table(processout, file=finalfile,row.names=FALSE)
                  
                  					## add missing rows if option is TRUE
                  					if (fillIncompleteRows) {
                    
                    					tempTogetfeature <- work.h[which(work.h$FEATURE %in% names(finalfeatureID)), ]
                    
                    					## get PROTEIN and FEATURE infomation
                    					tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
                    
                    					## merge feature info and run info as 'work' format
                    					tempmissingwork <- data.frame(tempfeatureID, LABEL="H",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
                    
                    					## merge with tempary space, missingwork
                    					missingwork.h <- rbind(missingwork.h, tempmissingwork)
                  					} # end fillIncompleteRows options
                				} # end loop for run ID
              				} # end for endogenous
            			} # end any missingness
					} # end label-based
        		} # if only any flag for method
      		} # end loop for methods
      
			if (fillIncompleteRows) {
        
        		## merge with work
        		## in future, use rbindlist?? rbindlist(list(work, missingwork))
        		if (nlevels(worktemp$LABEL) == 1) {
         			work <- rbind(work, missingwork)
        		} else {
          			work <- rbind(work, missingcomplete.l, missingcomplete.h, missingwork.l, missingwork.h)
        		}	
        
        		## print message
        		message("\n DONE : Incomplete rows for missing peaks are added with intensity values=NA. \n")
        
        		## save in process file.
        		processout <- rbind(processout, "Incomplete rows for missing peaks are added with intensity values=NA. - done, Okay")
        		write.table(processout, file=finalfile,row.names=FALSE)
        
      		} else if (!is.null(missingcomplete.l) | !is.null(missingcomplete.h) | !is.null(missingwork.l) | !is.null(missingwork.l) | !is.null(missingwork)) {
        
        		## save in process file.
        		processout <- rbind(processout,"Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        		write.table(processout, file=finalfile,row.names=FALSE)
        
        		stop("Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        
			}  
      
    	} else {
      		processout <- rbind(processout,c("Balanced data format with NA for missing feature intensities - okay"))
      		write.table(processout, file=finalfile,row.names=FALSE)
    	}
    
    ## for duplicate, in future
    
	} # end multiple method
  

	## factorize GROUP, SUBJECT, GROUP_ORIGINAL, SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN
  	## -------------------------------------------------------------------------------------------------

	work$PROTEIN <- factor(work$PROTEIN)
	work$PEPTIDE <- factor(work$PEPTIDE)
	work$TRANSITION <- factor(work$TRANSITION)
  
	work <- work[with(work, order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, PROTEIN, PEPTIDE, TRANSITION)),]
  
	work$GROUP <- factor(work$GROUP)
	work$SUBJECT <- factor(work$SUBJECT)
	## SUBJECT_ORIGINAL_NESTED will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL
  
	work$SUBJECT_NESTED <- factor(work$SUBJECT_NESTED, levels=unique(work$SUBJECT_NESTED))
  
	## FEATURE will sorted as PROTEIN, PEPTIDE, TRANSITION
	work$FEATURE <- factor(work$FEATURE, levels=unique(work$FEATURE))

	## RUN will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN
	work$RUN <- factor(work$RUN, levels=unique(work$RUN), labels=seq(1, length(unique(work$RUN))))
  
	processout <- rbind(processout, c("Factorize in columns(GROUP, SUBJECT, GROUP_ORIGINAL, SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN) - okay"))
	write.table(processout, file=finalfile, row.names=FALSE)
  
  
	## Normalization ##
	## ------------- ##
	   
	## Normalization : option 0. none
	if (is.element("NONE",normalization) | is.element("FALSE",normalization)) { # after 'toupper', FALSE becomes character.
		processout <- rbind(processout, c("Normalization : no normalization - okay"))
		write.table(processout, file=finalfile, row.names=FALSE)
	}
  

	## Normalization : option 1. constant normalization , equalize medians ##
	## -------------------------------------------------------------------
	if (!is.element("NONE", normalization) & !is.element("FALSE", normalization) & is.element("EQUALIZEMEDIANS", normalization)) {
    
    #	if (!checkMultirun) {
    #		if (nlevels(work$LABEL)==1) {
    #			### Constant normalization by endogenous
    #			median.run <- tapply(work$ABUNDANCE,work[,"RUN"], function(x) median(x,na.rm=TRUE))
    #			median.all <- median(work$ABUNDANCE, na.rm=TRUE)
    
    #			for (i in 1:length(unique(work$RUN))) {
    #				## ABUNDANCE is normalized
    #				work[work$RUN==i,"ABUNDANCE"] <- work[work$RUN==i,"ABUNDANCE"]-median.run[i]+median.all
    #			}
    #		}
    
    #		if (nlevels(work$LABEL)==2) {
    ### Constant normalization by heavy standard
    #			h <- work[work$LABEL=="H",]
    #			median.run <- tapply(h$ABUNDANCE,h[,"RUN"], function(x) median(x,na.rm=TRUE))
    #			median.all <- median(h$ABUNDANCE, na.rm=TRUE)
    
    #			for (i in 1:length(unique(work$RUN))) {
    #				## ABUNDANCE is normalized
    #				work[work$RUN==i,"ABUNDANCE"] <- work[work$RUN==i,"ABUNDANCE"]-median.run[i]+median.all
    #			}		
    #		}
    #	}else{ ## for multi-method case
    	
		if (nlevels(work$LABEL) == 1) {
			## Constant normalization by endogenous per method
      
			## [MC : use median of medians]
			median.run.method  <-  aggregate(ABUNDANCE ~ RUN + METHOD, data = work, median, na.rm = TRUE)
			median.method  <-  tapply(median.run.method$ABUNDANCE, median.run.method$METHOD, median, na.rm = TRUE)
      
			nmethod <- unique(work$METHOD)
      
			for(j in 1:length(nmethod)) {
				namerun <- unique(work[work$METHOD == nmethod[j], "RUN"])
        
				for (i in 1:length(namerun)) {
					## ABUNDANCE is normalized
					work[work$RUN == namerun[i], "ABUNDANCE"] <- work[work$RUN == namerun[i], "ABUNDANCE"] - median.run.method[median.run.method$RUN == namerun[i], "ABUNDANCE"] + median.method[j]
				}
			}
		}
    
    
		if (nlevels(work$LABEL) ==2 ) {
      
			## Constant normalization by heavy standard per method
			h <- work[work$LABEL == "H", ]
      		      		
      		## [MC : use median of medians]
			median.run.method  <-  aggregate(ABUNDANCE ~ RUN + METHOD, data = h, median, na.rm = TRUE)
			median.method  <-  tapply(median.run.method$ABUNDANCE, median.run.method$METHOD, median, na.rm = TRUE)
      
      		nmethod <- unique(work$METHOD)
      
      		for(j in 1:length(nmethod)) {
        		namerun <- unique(work[work$METHOD==nmethod[j],"RUN"])
        
        		for (i in 1:length(namerun)) {
          			## ABUNDANCE is normalized
          			work[work$RUN == namerun[i], "ABUNDANCE"] <- work[work$RUN == namerun[i], "ABUNDANCE"] - median.run.method[median.run.method$RUN == namerun[i], "ABUNDANCE"] + median.method[j]
        		}
      		} # end loop method		
    	} # for labe-based 
    
    
		processout <- rbind(processout,c("Normalization : Constant normalization (equalize medians) - okay"))
    	write.table(processout, file=finalfile,row.names=FALSE)
  	} ## end equaliemedian normalization	
  
	## Normalization : option 2. quantile normalization ##
	## ------------------------------------------------ ##
	if (!is.element("NONE", normalization) & !is.element("FALSE", normalization) & is.element("QUANTILE", normalization)) {
    
		if (nlevels(work$LABEL) == 1) {
      
			## for label-free, just use endogenous
      
      		nmethod <- unique(work$METHOD)
      		quantileall <- NULL
      
      		for (j in 1:length(nmethod)) {
         		namerun <- unique(work[work$METHOD == nmethod[j],"RUN"])
        
        		worktemp <- work[which(work$RUN %in% namerun & !is.na(work$INTENSITY)),]
       	 		worktemp$RUN <- factor(worktemp$RUN)
        		worktemp$FEATURE <- factor(worktemp$FEATURE)
        
        		quantiletemp <- as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp))
        
        		## need to put NA for missing value in endogenous
        		quantiletemp[quantiletemp == 0] <- NA
        
        		## using preprocessCore library
        		quantiledone <- normalize.quantiles(quantiletemp)
        		rownames(quantiledone) <- rownames(quantiletemp)
        		colnames(quantiledone) <- colnames(quantiletemp)
        
				## get quantiled to long format for apply difference endogenous
        		quantilelong <- melt(quantiledone, id=rownames(quantiledone))
        		colnames(quantilelong) <- c("FEATURE", "RUN", "ABUNDANCE_quantile")
        		rm(quantiledone)
        
        		## quantileall <- rbindlist(list(quantileall,quantilelong))
        		quantileall <- rbind(quantileall, quantilelong)
        
        		rm(quantilelong)
      		}
      
      		work <- merge(work, quantileall, by=c("FEATURE", "RUN"))
      		rm(quantileall)
      
      		## reorder
      		work <- data.frame("PROTEIN"=work$PROTEIN, "PEPTIDE"=work$PEPTIDE, "TRANSITION"=work$TRANSITION, "FEATURE"=work$FEATURE, "LABEL"=work$LABEL, "GROUP_ORIGINAL"=work$GROUP_ORIGINAL, "SUBJECT_ORIGINAL"=work$SUBJECT_ORIGINAL, "RUN"=work$RUN, "GROUP"=work$GROUP, "SUBJECT"=work$SUBJECT, "SUBJECT_NESTED"=work$SUBJECT_NESTED, "INTENSITY"=work$INTENSITY, "ABUNDANCE"=work$ABUNDANCE_quantile, "METHOD"=work$METHOD)
      
      		work <- work[with(work, order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, PROTEIN, PEPTIDE, TRANSITION)), ]
      
    	}
    
    	if (nlevels(work$LABEL) == 2) {
      
      		nmethod <- unique(work$METHOD)
      		quantileall <- NULL
      
      		for (j in 1:length(nmethod)) {
        		namerun <- unique(work[work$METHOD == nmethod[j], "RUN"])
        
        		## for label-based, make quantile normalization for reference
        		##worktemp <- work[which(work$RUN %in% namerun & work$LABEL=="H" & !is.na(work$INTENSITY)),] ## because for sparse of reference
        		worktemp <- work[which(work$RUN %in% namerun & work$LABEL == "H"),]
        		worktemp$RUN <- factor(worktemp$RUN)
        		worktemp$FEATURE <- factor(worktemp$FEATURE)
        
        		quantiletemp <- as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp))
        		rm(worktemp)
        
        		## need to put NA for missing value in endogenous
        		quantiletemp[quantiletemp==0] <- NA
        
       			## using preprocessCore library
        		quantiledone <- normalize.quantiles(quantiletemp)
        		rownames(quantiledone) <- rownames(quantiletemp)
        		colnames(quantiledone) <- colnames(quantiletemp)
        
        		## get quantiled to long format for apply difference endogenous
       			quantilelong.h <- melt(quantiledone, id=rownames(quantiledone))
        		colnames(quantilelong.h) <- c("FEATURE","RUN","ABUNDANCE_quantile")
        		quantilelong.h <- data.frame(quantilelong.h, LABEL="H")
        
        		## endogenous, in order to applying
        		##worktemp.l <- work[which(work$RUN %in% namerun & work$LABEL=="L" & !is.na(work$INTENSITY)),] ## because for sparse of reference
        		worktemp.l <- work[which(work$RUN %in% namerun & work$LABEL=="L"),]
        		worktemp.l$RUN <- factor(worktemp.l$RUN)
        		worktemp.l$FEATURE <- factor(worktemp.l$FEATURE)
        
        		quantiletemp.l <- as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp.l))
        		rm(worktemp.l)
        
        		## need to put NA for missing value in endogenous
        		quantiletemp.l[quantiletemp.l==0] <- NA
        
        		## apply the difference from reference
        		quantiledone.l <- quantiletemp.l-(quantiletemp-quantiledone)
        
        		## get quantiled to long format for apply difference endogenous
        		quantilelong.l <- melt(quantiledone.l, id=rownames(quantiledone.l))
        		colnames(quantilelong.l) <- c("FEATURE", "RUN", "ABUNDANCE_quantile")
        		quantilelong.l <- data.frame(quantilelong.l, LABEL="L")
        
        		rm(quantiletemp)
        		rm(quantiledone)
        		rm(quantiletemp.l)
        		rm(quantiledone.l)
        
        		# quantileall <- rbindlist(list(quantileall,quantilelong.h, quantilelong.l))
        		quantileall <- rbind(quantileall,quantilelong.h, quantilelong.l)
        
      		}
      
			## merge with original data
     		 work <- merge(work,quantileall, by=c("FEATURE","RUN","LABEL"))
      
      		## reorder
      		work <- data.frame("PROTEIN"=work$PROTEIN, "PEPTIDE"=work$PEPTIDE, "TRANSITION"=work$TRANSITION, "FEATURE"=work$FEATURE, "LABEL"=work$LABEL, "GROUP_ORIGINAL"=work$GROUP_ORIGINAL, "SUBJECT_ORIGINAL"=work$SUBJECT_ORIGINAL, "RUN"=work$RUN, "GROUP"=work$GROUP, "SUBJECT"=work$SUBJECT, "SUBJECT_NESTED"=work$SUBJECT_NESTED, "INTENSITY"=work$INTENSITY, "ABUNDANCE"=work$ABUNDANCE_quantile,"METHOD"=work$METHOD)
      
      		work <- work[with(work,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL,RUN,PROTEIN,PEPTIDE,TRANSITION)),]
      
    	}
    
    	processout <- rbind(processout, c("Normalization : Quantile normalization - okay"))
    	write.table(processout, file=finalfile, row.names=FALSE)
  	}
  
  
	## Normalization : option 3. global standards - for endogenous ##
	## ----------------------------------------------------------- ##
	if (!is.element("NONE", normalization) & !is.element("FALSE", normalization) & is.element("GLOBALSTANDARDS", normalization)) {
    
		work$RUN <- factor(work$RUN)
		combine <- data.frame(RUN=levels(work$RUN))
		allPeptide <- unique(work$PEPTIDE)
		allProtein <- unique(work$PROTEIN)
    
		for (i in 1:length(nameStandards)) {
      
			## if Peptides
			## namePeptide <- allPeptide[grep(nameStandards[i],allPeptide)] ## cannot grep for modified peptide sequence, [,],+ sign
      		namePeptide <- tempPeptide[tempPeptide$PEPTIDESEQUENCE == nameStandards[i], "PEPTIDE"]
      
			if (length(namePeptide)!=0) {
				tempStandard <- work[work$PEPTIDE == namePeptide,]
			} else {
        
        		## if Proteins
        		nameProtein <- allProtein[allProtein == nameStandards[i]] # if we use 'grep', can' find the proteins name with some symbol, such as 'sp|P30153|2AAA_HUMAN'
        
        		if (length(nameProtein)!=0) {
          			tempStandard <- work[work$PROTEIN==nameProtein,]
        		} else {
             		processout <- rbind(processout,c(paste("global standard peptides or proteins, ",nameStandards[i] ,", is not in dataset. Please check whether 'nameStandards' input is correct or not.")))
          			write.table(processout, file=finalfile,row.names=FALSE)
          
          			stop(paste("global standard peptides or proteins, ",nameStandards[i] ,", is not in dataset. Please check whether 'nameStandards' input is correct or not."))
        		}	
      		}
      
      		## here, by RUN, but need to check !!!
      		tempStandard <- tempStandard[tempStandard$GROUP!="0",]
      		tempStandard$RUN <- factor(tempStandard$RUN)
      
      		tempStandard <- tempStandard[!is.na(tempStandard$ABUNDANCE),]
      		meanStandard <- tapply(tempStandard$ABUNDANCE, tempStandard$RUN, function(x) mean(x, na.rm=TRUE))
      
      		meanStandard <- data.frame(RUN=names(meanStandard),meanStandard)
      		combine <- merge(combine, meanStandard, by="RUN", all=TRUE)
      		colnames(combine)[i+1] <- paste("meanStandard",i,sep="")
    	}
    
    	rownames(combine) <- combine$RUN
    	combine <- subset(combine, select=-c(RUN))
    
    	## get mean among global standards
    	allmean <- apply(combine,1, function(x) mean(x, na.rm=TRUE))
    	## allmean[is.na(allmean)] <- 0
    
    	allmeantemp <- data.frame(RUN=names(allmean),allmean)
    	allrun <- unique(work[,c("RUN","METHOD")])
    
    	allmeantemp <- merge(allmeantemp, allrun,by="RUN")
    	median.all <- tapply(allmeantemp$allmean, allmeantemp$METHOD, function(x) median(x,na.rm=TRUE))
    
    	## adjust
    	nmethod <- unique(work$METHOD)
    
    	for(j in 1:length(nmethod)) {
      		namerun <- unique(work[work$METHOD==nmethod[j], "RUN"])
      
      		for (i in 1:length(namerun)) {
        		## ABUNDANCE is normalized			
        		if (!is.na(allmean[names(allmean)==namerun[i]])) work[work$RUN==namerun[i] & work$LABEL=="L","ABUNDANCE"] <- work[work$RUN==namerun[i] & work$LABEL=="L","ABUNDANCE"]-allmean[names(allmean)==namerun[i]]+median.all[j]
      		}
   		} # end loop method
    
    
   	 	processout <- rbind(processout, c("Normalization : normalization with global standards protein - okay"))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
  	}
  
	## BetweenRunInterferenceScore
	## need to make new function
  
  	if (betweenRunInterferenceScore) {
    
    	## only output light
    
   		l <- subset(work,LABEL=="L")
    
    	## add ProtFeature and ProtPeptide, because the shared peptides appear in multiple proteins
    	l$ProtFeature <- paste(l$PROTEIN,l$FEATURE,sep="/")	
    	l$ProtPeptide <- paste(l$PROTEIN,l$PEPTIDE,sep="/")	
    
    	temp <- tapply(l$ABUNDANCE,l[,c("RUN","ProtPeptide")],function(x) mean(x,na.rm=TRUE))
    	temp1 <- data.frame(ProtPeptide=rep(colnames(temp),each=dim(temp)[1]),RUN=rep(rownames(temp),dim(temp)[2]),meanPEPTIDE=as.numeric(unlist(temp)))
    
    	temp2 <- merge(l[, c("PROTEIN", "PEPTIDE", "FEATURE", "ProtPeptide", "ProtFeature", "RUN", "ABUNDANCE")], temp1, by=c("ProtPeptide","RUN"))
    
    	temp3 <- temp2[!is.na(temp2$ABUNDANCE),]
    
    	temp4 <- tapply(rownames(temp3),temp3[,c("ProtFeature")], function(x) cor(temp3[x,"ABUNDANCE"],temp3[x,"meanPEPTIDE"]))
    
    	names <- unique(temp2[,c("PROTEIN","PEPTIDE","FEATURE","ProtFeature")])
    	names <- names[with(names,order(ProtFeature)),]
    	BetweenRunInterferenceFile <- data.frame(names[,c("PROTEIN","PEPTIDE","FEATURE")],BetweenRunInterferenceScore=temp4)
    
    	BetweenRunInterferenceFile <- BetweenRunInterferenceFile[with(BetweenRunInterferenceFile,order(PROTEIN,PEPTIDE,FEATURE)),]
    
   		write.table(BetweenRunInterferenceFile,file=paste(address,"BetweenRunInterferenceFile.txt",sep=""))
    
    	processout <- rbind(processout,c("Between Run Interference Score is calculated and saved in .csv file - okay"))
    	write.table(processout, file=finalfile,row.names=FALSE)
    
  	} else {
    	processout <- rbind(processout, c("Between Run Interference Score is not calculated."))
    	write.table(processout, file=finalfile, row.names=FALSE)
  	}
  
  	#Below two lines were merely for in-house testing and comparisons when needed
	#work.NoImpute <- work
	#AbundanceAfterImpute <- .Imputation(work, cutoffCensored, censoredInt, remove50missing, MBimpute, original_scale)
	
	
	## featureSubset ##
	## ------------- ##
  	##  !! need to decide how to present : keep original all data and make new column to mark, or just present selected subset    
  
	if (featureSubset == "all") {
 	 	message("* Use all features that the dataset origianally has.")
 	 
 	 	processout <- rbind(processout,c("* Use all features that the dataset origianally has."))
     	write.table(processout, file=finalfile, row.names=FALSE)
  	} 

	if (featureSubset == "highQuality") {
	  message("* Selecting high quality features temporarily defaults to featureSubset = top3. Updates for this option will be available in the next release.")
    
    featureSubset <- 'top3'
	  
	  processout <- rbind(processout, c("* Selecting high quality features temporarily defaults to featureSubset = top3. Updates for this option will be available in the next release."))

	  write.table(processout, file=finalfile, row.names=FALSE)
	  
	  #message("* Use feature selection algorithm in order to remove features with interference.")
	  
	  #processout <- rbind(processout,c("* Use feature selection algorithm in order to get high quality features."))
	  #write.table(processout, file=finalfile, row.names=FALSE)
	  
	  ### 2016.04.25. MC
	  ### there is the possibility to remain features which have completely missing in the certain condition after imputation
	  ### Therefore, remove the features which are completely missing in the certain condition before imputation
	  #removeissuefeature <- .getfeatureID(work)
	  
	  #if (length(removeissuefeature) != 0) {
	  #  work <- work[-which(work$FEATURE %in% removeissuefeature), ]
	  #}
	  
	  ###Impute the missing valuess before feature selection
	  #AbundanceAfterImpute <- .Imputation(work, cutoffCensored, censoredInt, MBimpute, remove50missing)
	  
	  #work <- AbundanceAfterImpute
	  
	  ### 
	  #work <- .feature_selection(work, remove_proteins_with_interference)
	  ##SelectionAfterImpute <- work
	  
	  ### 20160425-MC : after selecting feature, original ABUNCANCE should be used.
	  #work$ABUNDANCE <- work$ABUNDANCE.O
	  #work <- work[, -which(colnames(work) %in% c("ABUNDANCE.O", "feature.label", "run.label", "cen", "pred", "ref", "Protein_Peptide"))]
	}
  
	if (featureSubset == "top3") {
  		message("* Use top3 features that have highest average of log2(intensity) across runs.")
  		
  		processout <- rbind(processout, c("* Use top3 features that have highest average of log2(intensity) across runs."))
        write.table(processout, file=finalfile, row.names=FALSE)
	 
  	 	## INTENSITY vs ABUNDANCE? [THT: make more sense to use ABUNDANCE]
  		## how to decide top3 for DIA?
	
		temp1 <- aggregate(INTENSITY~PROTEIN+FEATURE,data=work, function(x) mean(x, na.rm=TRUE))

		temp2 <- split(temp1, temp1$PROTEIN)

		temp3 <- lapply(temp2, function(x) { 
			x <- x[order(x$INTENSITY, decreasing=T),]
			x <- x$FEATURE[1:3]
			})
	
		selectfeature <- unlist(temp3, use.names=FALSE)
		selectfeature <- selectfeature[!is.na(selectfeature)]
		
		## get subset
		work <- work[which(work$FEATURE %in% selectfeature),]	

  	}
  
	if (featureSubset == "topN") {
    
    ## check whether there is the input for 'N'
    
	  message(paste("* Use top", n_top_feature, " features that have highest average of log2(intensity) across runs.", sep=""))
	  
	  processout <- rbind(processout, c(paste("* Use top", n_top_feature, " features that have highest average of log2(intensity) across runs.", sep="")))
	  write.table(processout, file=finalfile, row.names=FALSE)
	  
	  ## INTENSITY vs ABUNDANCE? [THT: make more sense to use ABUNDANCE]
	  ## how to decide top3 for DIA?
	  
    worktemp <- work[!is.na(work$ABUNDANCE) & work$ABUNDANCE != 0, ]
	  temp1 <- aggregate(INTENSITY ~ PROTEIN+FEATURE, data=worktemp, function(x) mean(x, na.rm=TRUE))
	  
	  temp2 <- split(temp1, temp1$PROTEIN)
	  
	  temp3 <- lapply(temp2, function(x) { 
	    x <- x[order(x$INTENSITY, decreasing=T), ]
	    x <- x$FEATURE[1:n_top_feature]
	  })
	  
	  selectfeature <- unlist(temp3, use.names=FALSE)
	  selectfeature <- selectfeature[!is.na(selectfeature)]
	  
	  ## get subset
	  work <- work[which(work$FEATURE %in% selectfeature),]	
	  
	}
  
	## check missingness 
	## transitions are completely missing in one condition : missingness ##
	if (nlevels(work$LABEL) == 1) {
    #Use the data frame before imputation to summarize the missingness
    all.work <- work	
    test <- tapply(is.na(work[, "ABUNDANCE"]), work[, c("GROUP_ORIGINAL","FEATURE")], function(x) sum(x, na.rm=TRUE))
    numObs <- tapply(work[, "ABUNDANCE"], work[, c("GROUP_ORIGINAL","FEATURE")], function(x) length(x))
    test1 <- test == numObs
    test2 <- apply(test1, 2, function(x) sum(x, na.rm=TRUE))
    filterList <- names(test2)[test2 > 0]
    final.decision <- ifelse(test2>0, 1, 0)
	}	
  
	if (nlevels(work$LABEL) == 2) {
		
		#Use the data frame before imputation to summarize the missingness
    ## first, remove NA
    all.work <- work   # with all NA observations
    work.miss <- na.omit(work)
    
		## draw table
		light <- subset(work.miss,LABEL=="L")
		heavy <- subset(work.miss,LABEL=="H")
    
		## use FEATURE because the name of transition can be used in other peptide
		count.light <- xtabs(~FEATURE+GROUP_ORIGINAL, light)
		count.heavy <- xtabs(~FEATURE+GROUP_ORIGINAL, heavy)
    
		count.light <- count.light==0
		count.heavy <- count.heavy==0
    
		count.light <- as.data.frame(count.light)
		count.heavy <- as.data.frame(count.heavy)
    
		## summary of missingness
		decision <- count.light
		decision[] <- 0
    
		for (i in 1:ncol(decision)) {
			for (j in 1:nrow(decision)) {
        
				## either light or heavy has no obs -> subject to filter
        		if (count.light[j,i]==TRUE || count.heavy[j,i]==TRUE) { 
        			decision[j,i] <- 1 
        		}
      		}
    	}
    
    	final.decision <- apply(decision,1,sum)
    
    	## assign "subject to filter" column
    	work <- data.frame(work, "SuggestToFilter"=0)
    
    	for(i in 1:length(final.decision)) {
      		## assign subject_to_filter=1 for entire transition
      		if (final.decision[i]!=0) {
      			work[work$FEATURE==names(final.decision[i]), "SuggestToFilter"] <- 1
      		} 
		}	
	}
  
	## output : summary ##
	## ---------------- ##

	## output for label
  	processout <- rbind(processout, c(paste(length(unique(work$LABEL)), " level of Isotope type labeling in this experiment", sep="")))
	write.table(processout, file=finalfile, row.names=FALSE)
  
	temp <- data.frame("Summary of Features :")
	colnames(temp) <- " "
	rownames(temp) <- " "
	print(temp)
  
	summary.f <- matrix(NA,nrow=3)
	summary.f[1] <- nlevels(work$PROTEIN)
  
	temp <- unique(work[, c("PROTEIN", "PEPTIDE")])
	temp1 <- xtabs(~PROTEIN, data=temp)
	temp2 <- summary(as.numeric(temp1))
	summary.f[2] <- paste(temp2["Min."], temp2["Max."], sep="-")
  
	temp <- unique(work[, c("PEPTIDE", "FEATURE")])
	temp1 <- xtabs(~PEPTIDE, data=temp)
	temp2 <- summary(as.numeric(temp1))
	summary.f[3] <- paste(temp2["Min."], temp2["Max."], sep="-")
  
	colnames(summary.f) <- "count"
	rownames(summary.f) <- c("# of Protein", "# of Peptides/Protein", "# of Transitions/Peptide")
  
	print(as.data.frame(summary.f))
  
	## output for process
	processout <- rbind(processout, c("Summary of Features :"))
	processout <- rbind(processout, c(paste(rownames(summary.f)[1]," : ", summary.f[1], sep="")))
	processout <- rbind(processout, c(paste(rownames(summary.f)[2]," : ", summary.f[2], sep="")))
	processout <- rbind(processout, c(paste(rownames(summary.f)[3]," : ", summary.f[3], sep="")))
  
	write.table(processout, file=finalfile, row.names=FALSE)
  
	## protein list with 1 feature
	temp <- unique(work[, c("PROTEIN", "FEATURE")])
	temp1 <- xtabs(~PROTEIN, data=temp)
	temp2 <- as.data.frame(temp1[temp1 == 1])
	if (nrow(temp2) > 0) {
		message("\n","** Protein (",paste(rownames(temp2),collapse = ", "),") has only single transition : Consider excluding this protein from the dataset.", "\n")
	}
  
	temp <- data.frame("Summary of Samples :")
	colnames(temp) <- " "
	rownames(temp) <- " "
	print(temp)
  
	summary.s <- matrix(NA,ncol=nlevels(work$GROUP_ORIGINAL),nrow=3)
  
	## # of MS runs
	temp <- unique(work[, c("GROUP_ORIGINAL", "RUN")])
	temp1 <- xtabs(~GROUP_ORIGINAL, data=temp)
	summary.s[1,] <- temp1
  
	## # of biological replicates
	temp <- unique(work[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
	temp1 <- xtabs(~GROUP_ORIGINAL, data=temp)
	summary.s[2,] <- temp1
  
	## # of technical replicates
	c.tech <- round(summary.s[1,] / (summary.s[2,] * length(unique(work$METHOD))))
	##summary.s[3,] <- ifelse(c.tech==1,0,c.tech)
	summary.s[3,] <- c.tech
  
	colnames(summary.s) <- unique(work$GROUP_ORIGINAL)
	rownames(summary.s) <- c("# of MS runs","# of Biological Replicates", "# of Technical Replicates")
  
	print(summary.s)
  
	message("\n Summary of Missingness :\n" )
	message("  # transitions are completely missing in one condition: ", sum(final.decision!=0), "\n")
	if (sum(final.decision!=0)!=0) {
		message("    -> ", paste(names(final.decision[final.decision!=0]),collapse = ", "))
	}
  
	without <- xtabs(~RUN, work)
	withall <- xtabs(~RUN, all.work)
	run.missing <- without / withall
	message("\n  # run with 75% missing observations: ", sum(run.missing<0.25), "\n")
	if (sum(run.missing<0.25)!=0) {
		message("    -> ", paste("RUN",names(without[run.missing<0.25]),sep=" "))
	}
  
	## output process
	processout <- rbind(processout,  c("Summary of Missingness :"))
	processout <- rbind(processout,c(paste("  # transitions are completely missing in one condition: ", sum(final.decision!=0), sep="")))
	if (sum(final.decision!=0)!=0){
		processout <- rbind(processout,"    -> ", paste(names(final.decision[final.decision != 0]), collapse = ", "))
	} 
  
	processout <- rbind(processout, c(paste("  # run with 75% missing observations: ", sum(run.missing < 0.25), sep="")))
	if (sum(run.missing<0.25)!=0) {
		processout <- rbind(processout, "    -> ", paste("RUN", names(without[run.missing < 0.25]), sep=" "))
  	}
  	write.table(processout, file=finalfile, row.names=FALSE)
  
	## check any protein has only light for labeled-experiment
	if (nlevels(work$LABEL) == 2) {
    	temp <- unique(work[, c("PROTEIN", "LABEL")])
    	temp1 <- xtabs(~PROTEIN, data=temp)
    
    	if (any(temp1 != 2)) {
      		## check that is L or H
      		namepro <- names(temp1[temp1!=2])
      		
      		for(j in 1:length(namepro)) {
      			
        		if (unique(work[work$PROTEIN == namepro[j], "LABEL"]) == "L") {
          			message("\n *** ", namepro[j], " has only endogeneous intensities in label-based experiment. Please check this protein or remove it.")
          		}
        
        		if (unique(work[work$PROTEIN == namepro[j], "LABEL"]) == "H") {
          			message("\n *** ", namepro[j], " has only reference intensities in label-based experiment. Please check this protein or remove it.")
         		}
      		}
    	}	
	}
  
	processout <- rbind(processout, c("Processing data for analysis is done. - okay"))
	write.table(processout, file=finalfile, row.names=FALSE)
  
	## after normalization, zero intensity could be negative
  work[!is.na(work$ABUNDANCE) & work$ABUNDANCE < 0, "ABUNDANCE"] <- 0
    
  work[!is.na(work$INTENSITY) & work$INTENSITY == 0, "ABUNDANCE"] <- 0


 	## get the summarization per subplot (per RUN)   
	## -------------------------------------------
	
	message("\n == Start the summarization per subplot...")

	rqresult <- try(.runQuantification(work, summaryMethod, equalFeatureVar, filterLogOfSum, cutoffCensored, censoredInt, remove50missing, MBimpute, original_scale, logsum, featureSubset), silent=TRUE)

	if (class(rqresult) == "try-error") {
		message("*** error : can't summarize per subplot with ", summaryMethod, ".")
     
      	processout <- rbind(processout,c(paste("error : can't summarize per subplot with ", summaryMethod, ".", sep = "")))
      	write.table(processout, file=finalfile, row.names=FALSE)
	
	  	rqall <- NULL
	  	rqmodelqc <- NULL
	  	workpred <- NULL
	  
	 } else {
      
		label <- nlevels(work$LABEL) == 2
	
		if (sum(is.element(colnames(rqresult$rqdata), "RUN")) == 0) {
			## logsum is summarization per subject
			lab <- unique(work[, c("GROUP", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT")])
	
			if (label) lab <- lab[lab$GROUP != 0,]

			rqall <- merge(rqresult$rqdata, lab, by="SUBJECT_ORIGINAL")
			
		} else {
			lab <- unique(work[, c("RUN", "GROUP", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT")])
	
			if (label) lab <- lab[lab$GROUP != 0,]

			rqall <- merge(rqresult$rqdata, lab, by="RUN")
		}
		
		rqall$GROUP <- factor(rqall$GROUP)
		rqall$Protein <- factor(rqall$Protein)
		
		rqmodelqc <- rqresult$ModelQC
		
		workpred <- rqresult$PredictedBySurvival
		
		message("\n == the summarization per subplot is done.")
		
		processout <- rbind(processout, c(paste("the summarization per subplot is done.- okay : ", summaryMethod, sep="")))
		write.table(processout, file=finalfile, row.names=FALSE)

	 }
	 
  	## return work data.frame	and run quantification
	
	#Align the run quantification data
	rqall <- rqall[order(rqall$Protein, as.numeric(as.character(rqall$RUN))),]
    rownames(rqall) <- NULL
	
	#Mike: Below is for in-house verification occasionally
	#processedquant <- list(ProcessedData=work.NoImpute, RunlevelData=rqall, SummaryMethod=summaryMethod, ModelQC=rqmodelqc, PredictBySurvival=workpred, ImputedData=AbundanceAfterImpute)
	processedquant <- list(ProcessedData=work, RunlevelData=rqall, SummaryMethod=summaryMethod, ModelQC=rqmodelqc, PredictBySurvival=workpred)
	
    return(processedquant)
  
}




########################################################
.runQuantification <- function(data, summaryMethod, equalFeatureVar, filterLogOfSum, cutoffCensored, censoredInt, remove50missing, MBimpute, original_scale, logsum, featureSubset) {
	
    ##Since the imputation has been done before feature selection, delete the columns of censoring indicator to avoid imputing the same intensity again	
    #if(featureSubset == "highQuality") {
    #	data$cen <- NULL; data$pred <- NULL; data$INTENSITY <- 2^data$ABUNDANCE
    #} 
    
    ##If we want to impute again after the feature selection
    #if(featureSubset == "highQuality" & ImputeAgain==TRUE) {
    #	data$ABUNDANCE <- data$ABUNDANCE.O	
    #}

    data$LABEL <- factor(data$LABEL)
    label <- nlevels(data$LABEL)==2
    
   	# set ref which is distinguish reference and endogenous. any reference=0. endogenous is the same as RUN
	if (label) {
		data$ref <- 0
		data$ref[data$LABEL!="H"] <- data$RUN[data$LABEL!="H"]
		data$ref <- factor(data$ref)
#		unique(data[,c("RUN","LABEL","GROUP","ref")])
	}
	      
#    finalresult <- data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$RUN)),RUN=rep(c(levels(data$RUN)),nlevels(data$PROTEIN)),Condition=NA, BioReplicate=NA,LogIntensities=NA,NumFeature=NA,NumPeaks=NA)

	# for saving predicting value for impute option
	predAbundance <- NULL
	
	###################################
	## method 1 : model based summarization
	if (summaryMethod=="linear"  & is.null(censoredInt)) {
		
		data <- data[!is.na(data$ABUNDANCE),]
    	data$PROTEIN <- factor(data$PROTEIN)
    	data$RUN <- factor(data$RUN)
    
		result <- NULL
	    dataafterfit <- NULL
	    
		for(i in 1: nlevels(data$PROTEIN)) {
        
     		sub <- data[data$PROTEIN==levels(data$PROTEIN)[i],]

#      sub$GROUP <- factor(sub$GROUP)
#      sub$SUBJECT <- factor(sub$SUBJECT)
#      sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
#      sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)
     		sub$SUBJECT_NESTED <- factor(sub$SUBJECT_NESTED)
      		sub$FEATURE <- factor(sub$FEATURE)	
      		sub$RUN <- factor(sub$RUN)	        
      
      		if (!label) {
      			temp <- data.frame(xtabs(~RUN, data=sub))
	
      			sub.result <- data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$RUN)),RUN=rep(c(levels(sub$RUN)),1),LogIntensities=NA, NumFeature=length(unique(sub$FEATURE)),NumPeaks=temp$Freq)
      	
      		}else{
      	
      			sub$ref <- factor(sub$ref)			

      			temp <- data.frame(xtabs(~ref, data=sub))

      			sub.result <- data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$ref)),RUN=rep(c(levels(sub$ref)[-1],"Ref"),1),LogIntensities=NA, NumFeature=length(unique(sub$FEATURE)),NumPeaks=c(temp[-1,"Freq"],temp[1,"Freq"]))
      	
      		}
      
      		singleFeature <- .checkSingleFeature(sub)
      		singleSubject <- .checkSingleSubject(sub)
      		TechReplicate <- .checkTechReplicate(sub) ## use for label-free model
        
      		##### fit the model
      		message(paste("Getting the summarization per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
        
      		fit <- try(.fit.quantification.run(sub,singleFeature,singleSubject, TechReplicate,labeled=label, equalFeatureVar), silent=TRUE)
             
      		if (class(fit)=="try-error") {
        		message("*** error : can't fit the model for ", levels(data$PROTEIN)[i])
          
        		result <- rbind(result, sub.result)
          
         		if (nrow(sub)!=0) {
        			sub$residuals <- NA
        			sub$fitted <- NA
      			}
          
      		}else{
          
       			if (class(fit)=="lm") {
          			cf  <-  summary(fit)$coefficients
       	 		}else{
         			cf  <-  fixef(fit)
         		}
          
        		# calculate sample quantification for all levels of sample
        		a=1	
          
        		for(j in 1:nlevels(sub$RUN)) {
         			contrast.matrix <- rep(0,nlevels(sub$RUN))
          			contrast.matrix[j] <- 1
            
          			contrast <- .make.contrast.run.quantification(fit,contrast.matrix,sub, labeled=label)
            
         			if (class(fit)=="lm") {
           				sub.result[a,3] <- .estimableFixedQuantification(cf,contrast)
          			}else{
            			sub.result[a,3] <- .estimableRandomQuantification(cf,contrast)
          			}
          
          			a=a+1
        		}
          
          		## for label-based case, need reference quantification
        		if (label) {
        			contrast <- .make.contrast.run.quantification.reference(fit,contrast.matrix,sub)
        			if (class(fit)=="lm") {
              			sub.result[a,3] <- .estimableFixedQuantification(cf,contrast)
          			}else{
              			sub.result[a,3] <- .estimableRandomQuantification(cf,contrast)
          			}
        		}
          
        		result <- rbind(result, sub.result)

 				if (class(fit)=="lm") {  ### lm model
 					sub$residuals <- fit$residuals
      				sub$fitted <- fit$fitted.values
  				}else{   ### lmer model
    				sub$residuals <- resid(fit)
    				sub$fitted <- fitted(fit)
  				}
  			
     			dataafterfit <- rbind(dataafterfit,sub)

          	}
        
      	} ## end-loop for each protein	
	} ## for linear model summary
	
	###################################
	## Method 2 : Tukey Median Polish	
	if (summaryMethod=="TMP") {
		
		#data <- data[!is.na(data$ABUNDANCE),]
   	 	data$PROTEIN <- factor(data$PROTEIN)
    	data$RUN <- factor(data$RUN)
    
		result <- NULL
	  
		for(i in 1: nlevels(data$PROTEIN)) {
              		
     		sub <- data[data$PROTEIN==levels(data$PROTEIN)[i],]
     		
     		message(paste("Getting the summarization by Tukey's median polish per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

      		sub$FEATURE <- factor(sub$FEATURE)	
      		sub$feature.label <- paste(sub$FEATURE, sub$LABEL, sep="_")
      		sub$run.label <- paste(sub$RUN, sub$LABEL, sep="_")
      		
      		## if all measurements are NA,
      		if (nrow(sub)==sum(is.na(sub$ABUNDANCE))) {
       			message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),") because all measurements are NAs."))
        		next()
      		}
      		
      		## remove features which are completely NAs
			subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
			countfeature <- xtabs(~FEATURE, subtemp)
			namefeature <- names(countfeature)[countfeature==0]
				
			if (length(namefeature)!=0) {
				sub <- sub[-which(sub$FEATURE %in% namefeature), ]
				
				if (nrow(sub) == 0) {
					message(paste("Can't summarize for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(data$PROTEIN)), ") because all measurements are NAs."))
        			next()
        		
				} else {
					sub$FEATURE <- factor(sub$FEATURE)
				}
			}
			
			## remove features which have only 1 measurement.
			namefeature1 <- names(countfeature)[countfeature == 1]
				
			if (length(namefeature1)!=0) {
				sub <- sub[-which(sub$FEATURE %in% namefeature1), ]
				 
				if (nrow(sub) == 0) {
					message(paste("Can't summarize for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(data$PROTEIN)), ") because features have only one measurement across MS runs."))
        			next()
        		
				} else {
					sub$FEATURE <- factor(sub$FEATURE)
				}
			}
      		
      		## remove run which has no measurement at all 
			subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
			count <- aggregate(ABUNDANCE~RUN,data=subtemp, length)
			norun <- setdiff(unique(data$RUN),count$RUN)
				
			if (length(norun)!=0 & length(intersect(norun, as.character(unique(sub$RUN))))) { # removed NA rows already, if there is no overlapped run, error
				sub <- sub[-which(sub$RUN %in% norun),]
				sub$RUN <- factor(sub$RUN)
			}

			
			if (remove50missing) {
				# count # feature per run
        			if (!is.null(censoredInt)) {
						if (censoredInt=="NA") {
							subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY),]
						}
					
						if (censoredInt=="0") {
							subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						}
					}
					
					numFea <- xtabs(~RUN, subtemp) ## RUN or run.label?
					numFea <- numFea/length(unique(subtemp$FEATURE))
					numFea <- numFea<=0.5
					removerunid <- names(numFea)[numFea]
					
					## if all measurements are NA,
      				if (length(removerunid)==length(numFea)) {
       					message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),") because all runs have more than 50% NAs and are removed with the option, remove50missing=TRUE."))
        				next()
      				}

			}
			
			##### how to decide censored or not
			if (!is.null(censoredInt)) {
				## 1. censored 
				if (censoredInt=="0") {
					sub$cen <- ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY==0,0,1)
				}
				
				### 2. all censored missing
				if (censoredInt=="NA") {
					sub$cen <- ifelse(is.na(sub$INTENSITY),0,1)
				}
				
				### check whether we need to impute or not.
				if (sum(sub$cen==0)>0) {
			
					##### cutoffCensored
					## 1. put 0 to censored
					#if (cutoffCensored=="0") {
					#	if (censoredInt=="NA") {
					#		sub[is.na(sub$INTENSITY),"ABUNDANCE"] <- 0
					#	}
				
					#	if (censoredInt=="0") {
					#		sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- 0
					#	}
					#}
			
					## 2. put minimum in feature level to NA
					if (cutoffCensored=="minFeature") {
						if (censoredInt=="NA") {
							cut <- aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
							## cutoff for each feature is little less than minimum abundance in a run.
							cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
						
							## remove runs which has more than 50% missing values
							if (remove50missing) {
								if (length(removerunid)!=0) {
									sub <- sub[-which(sub$RUN %in% removerunid),]
									sub$RUN <- factor(sub$RUN)
								}
							}

							for(j in 1:length(unique(cut$feature.label))) {
								sub[is.na(sub$INTENSITY) & sub$feature.label==cut$feature.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
							}
						}
					
						if (censoredInt=="0") {
							subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
							cut <- aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
							## cutoff for each feature is little less than minimum abundance in a run.
							cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
						
						
							## remove runs which has more than 50% missing values
							if (remove50missing) {
								if (length(removerunid)!=0) {
									sub <- sub[-which(sub$RUN %in% removerunid),]
									sub$RUN <- factor(sub$RUN)
								}
							}
						
							for(j in 1:length(unique(cut$feature.label))) {
								sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0  & sub$feature.label==cut$feature.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
							}
						}
					}
				
					## 3. put minimum in RUN to NA
					if (cutoffCensored=="minRun") {
				
						## remove runs which has more than 50% missing values
						if (remove50missing) {
							if (length(removerunid)!=0) {
								sub <- sub[-which(sub$RUN %in% removerunid),]
								sub$RUN <- factor(sub$RUN)
							}
						}
						
						if (censoredInt=="NA") {
							cut <- aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
							## cutoff for each Run is little less than minimum abundance in a run.
							cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

							for(j in 1:length(unique(cut$run.label))) {
								sub[is.na(sub$INTENSITY) & sub$run.label==cut$run.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
							}
						}
					
						if (censoredInt=="0") {
							subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
							cut <- aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
							cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

							for(j in 1:length(unique(cut$run.label))) {
								sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$run.label==cut$run.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
							}
						}
					}
				
					## 20150829 : 4. put minimum RUN and FEATURE
					if (cutoffCensored=="minFeatureNRun") {
						if (censoredInt=="NA") {
						
							## cutoff for each feature is little less than minimum abundance in a run.

							cut.fea <- aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
							cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
						
							## remove runs which has more than 50% missing values
							## before removing, need to contribute min feature calculation
							if (remove50missing) {
								if (length(removerunid)!=0) {
									sub <- sub[-which(sub$RUN %in% removerunid),]
									sub$RUN <- factor(sub$RUN)
								}
							}
						
							## cutoff for each Run is little less than minimum abundance in a run.

							cut.run <- aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
							cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
						
						
							if (length(unique(cut.fea$feature.label))>1) {
								for(j in 1:length(unique(cut.fea$feature.label))) {
									for(k in 1:length(unique(cut.run$run.label))) {
										# get smaller value for min Run and min Feature
										finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
										sub[is.na(sub$INTENSITY) & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"] <- finalcut
									}
								}
							}
							# if single feature, not impute
						}
					
						if (censoredInt=="0") {
							subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]

							cut.fea <- aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
							cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
												
							## remove runs which has more than 50% missing values
							## before removing, need to contribute min feature calculation
							if (remove50missing) {
								if (length(removerunid)!=0) {
									sub <- sub[-which(sub$RUN %in% removerunid),]
									sub$RUN <- factor(sub$RUN)
								}
							}

							cut.run <- aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
							cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
	
							if (length(unique(cut.fea$feature.label))>1) {
								for(j in 1:length(unique(cut.fea$feature.label))) {
									for(k in 1:length(unique(cut.run$run.label))) {
										# get smaller value for min Run and min Feature
										finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
										sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"] <- finalcut
									}
								}
							}else{ # single feature
					
								sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut.fea$ABUNDANCE
							
							}
						}
					}
				
					if (MBimpute) {
					
						if (nrow(sub[sub$cen==0,])>0) {
							## impute by survival model
							subtemp <- sub[!is.na(sub$ABUNDANCE),]
							countdf <- nrow(subtemp)<(length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
				
							### fit the model
							if (length(unique(sub$FEATURE))==1) {
								fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
							}else{
								if (countdf) {
									fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
								}else{
									fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN,data=sub, dist='gaussian')
								}
							}
					
							# get predicted value from survival
							sub <- data.frame(sub, pred=predict(fittest, newdata=sub, type="response"))
					
							# the replace censored value with predicted value
							sub[sub$cen==0,"ABUNDANCE"] <- sub[sub$cen==0,"pred"]	
					
							# save predicted value
							predAbundance <- c(predAbundance,predict(fittest, newdata=sub, type="response"))				
						}
					}	
				}
			}
			
			## then, finally remove NA in abundance
			sub <- sub[!is.na(sub$ABUNDANCE),]
					       
					        
      		if (nlevels(sub$FEATURE)>1) { ## for more than 1 features
      			if (!label) { ## label-free
      			
      				data_w = dcast(RUN ~ FEATURE, data=sub, value.var='ABUNDANCE', keep=TRUE)
  					rownames(data_w) <- data_w$RUN
  					data_w <- data_w[,-1]
  					data_w[data_w==1] <- NA
  					
  					if (!original_scale) {
  						
  						meddata  <-  medpolish(data_w,na.rm=TRUE,trace.iter = FALSE)
						tmpresult <- meddata$overall + meddata$row
						
  					} else { # original_scale
  						data_w <- 2^(data_w)
  						
  						if (logsum) {
  							meddata  <-  medpolish(data_w,na.rm=TRUE,trace.iter = FALSE)
							tmpresult <- data_w - meddata$residuals
							tmpresult <- apply(tmpresult,1, function(x) sum(x, na.rm=TRUE))
							
  						} else {
  							meddata  <-  medpolish(data_w,na.rm=TRUE,trace.iter = FALSE)
							tmpresult <- meddata$overall + meddata$row
						}
						
  						tmpresult <- log2(tmpresult)
  					}
  					
					# count # feature per run
					if (!is.null(censoredInt)) {
						if (censoredInt=="NA") {
							subtemp <- sub[!is.na(sub$INTENSITY),]
							subtempimpute <- sub[is.na(sub$INTENSITY),]
							subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE), ]
						}
					
						if (censoredInt=="0") {
							subtemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
							subtempimpute <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,]
							subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE) & subtempimpute$ABUNDANCE!=0, ]
						}
						
						numFea <- xtabs(~RUN, subtemp)
						numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
						numFeaTF <- numFeaPercentage >= 0.5
					
						numimpute <- xtabs(~RUN, subtempimpute)
					
						sub.result <- data.frame(Protein=unique(sub$PROTEIN),LogIntensities=tmpresult, RUN=names(tmpresult), NumMeasuredFeature = as.vector(numFea), MissingPercentage=as.vector(numFeaPercentage), more50missing=numFeaTF, NumImputedFeature = as.vector(numimpute))
					
					}else{
						subtemp <- sub[!is.na(sub$INTENSITY),]
						
						numFea <- xtabs(~RUN, subtemp)
						numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
						numFeaTF <- numFeaPercentage >= 0.5
										
						sub.result <- data.frame(Protein=unique(sub$PROTEIN),LogIntensities=tmpresult, RUN=names(tmpresult), NumMeasuredFeature = as.vector(numFea), MissingPercentage=as.vector(numFeaPercentage), more50missing=numFeaTF)
						
					}
								
      				result <- rbind(result, sub.result)
      			
      			}else{ ## labeled
					
					data_w = dcast(run.label ~ FEATURE, data=sub, value.var='ABUNDANCE', keep=TRUE)
  					rownames(data_w) <- data_w$run.label
  					data_w <- data_w[,-1]
  					#data_w[data_w==1] <- NA
  					
  					meddata  <-  medpolish(data_w,na.rm=TRUE,trace.iter = FALSE)
					tmpresult <- meddata$overall + meddata$row
					
					reformresult <- data.frame(tmpresult)
					end <- nchar(rownames(reformresult))
					reformresult$LABEL <- substr(rownames(reformresult),end,end)
					reformresult$RUN <- substr(rownames(reformresult),1,end-2)
					colnames(reformresult)[1] <- "ABUNDANCE"
					
					## now single feature, adjust reference feature difference
      				h <- reformresult[reformresult$LABEL=="H",]
 					allmed <- median(h$ABUNDANCE, na.rm=TRUE)

        			for (i in 1:length(unique(reformresult$RUN))) {
          				## ABUNDANCE is normalized
          				reformresult[reformresult$RUN==unique(reformresult$RUN)[i],"ABUNDANCE"] <- reformresult[reformresult$RUN==unique(reformresult$RUN)[i],"ABUNDANCE"]-reformresult[reformresult$RUN==unique(reformresult$RUN)[i] & reformresult$LABEL=="H","ABUNDANCE"]+allmed
        			}
        			
        			reformresult <- reformresult[reformresult$LABEL=="L",]
        			
        			subtemp <- reformresult[!is.na(reformresult$ABUNDANCE),]
        			
        			# count # feature per run
        			if (!is.null(censoredInt)) {
						if (censoredInt=="NA") {
							subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY), ]
							subtempimpute <- sub[sub$LABEL=="L" & is.na(sub$INTENSITY),]
							subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE), ]
						}
					
						if (censoredInt=="0") {
							subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0, ]
							subtempimpute <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY==0, ]
							subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE) & subtempimpute$ABUNDANCE!=0, ]
						}
						
      					numFea <- xtabs(~RUN, subtemp)
						numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
						numFeaTF <- numFeaPercentage >= 0.5
					
						numimpute <- xtabs(~RUN, subtempimpute)
					
						sub.result <- data.frame(Protein=unique(sub$PROTEIN),LogIntensities=reformresult$ABUNDANCE, RUN=reformresult$RUN, NumMeasuredFeature = as.vector(numFea), MissingPercentage=as.vector(numFeaPercentage), more50missing=numFeaTF, NumImputedFeature = as.vector(numimpute))

					}else{
						subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY), ]
						
						numFea <- xtabs(~RUN, subtemp)
						numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
						numFeaTF <- numFeaPercentage >= 0.5
					
						sub.result <- data.frame(Protein=unique(sub$PROTEIN),LogIntensities=reformresult$ABUNDANCE, RUN=reformresult$RUN, NumMeasuredFeature = as.vector(numFea), MissingPercentage=as.vector(numFeaPercentage), more50missing=numFeaTF)

					}
      				
      				result <- rbind(result, sub.result)
      			}
      		
      		}else{ ## single feature
      			if (label) { ## label-based
      				## single feature, adjust reference feature difference
      				h <- sub[sub$LABEL=="H",]
 					allmed <- median(h$ABUNDANCE, na.rm=TRUE)

        			for (i in 1:length(unique(sub$RUN))) {
          				## ABUNDANCE is normalized
          				sub[sub$RUN==unique(sub$RUN)[i],"ABUNDANCE"] <- sub[sub$RUN==unique(sub$RUN)[i],"ABUNDANCE"]-sub[sub$RUN==unique(sub$RUN)[i] & sub$LABEL=="H","ABUNDANCE"]+allmed
        			}
        			
        			sub <- sub[sub$LABEL=="L",]
      			}
      			
				## single feature, use original values
				
      			subtemp <- sub[!is.na(sub$ABUNDANCE),]
      			
      			if (!is.null(censoredInt)) {
      				if (censoredInt=="NA") {
						subtempcount <- sub[!is.na(sub$INTENSITY),]
						subtempimpute <- sub[is.na(sub$INTENSITY),]
						subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE), ]
					}
					
					if (censoredInt=="0") {
						subtempcount <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						subtempimpute <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,]
						subtempimpute <- subtempimpute[!is.na(subtempimpute$ABUNDANCE) & subtempimpute$ABUNDANCE!=0, ]
					}
						
					numFea <- xtabs(~RUN, subtempcount)
					numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
					numFeaTF <- numFeaPercentage >= 0.5
					
					numimpute <- xtabs(~RUN, subtempimpute)
			
					sub.result <- data.frame(Protein=subtemp$PROTEIN,LogIntensities=subtemp$ABUNDANCE, RUN=subtemp$RUN, NumMeasuredFeature = as.vector(numFea), MissingPercentage=as.vector(numFeaPercentage), more50missing=numFeaTF, NumImputedFeature = as.vector(numimpute))
      			
				}else{
					subtempcount <- subtemp
					
					numFea <- xtabs(~RUN, subtempcount)
					numFeaPercentage <- 1 - numFea / length(unique(subtemp$FEATURE))
					numFeaTF <- numFeaPercentage >= 0.5
					
      				sub.result <- data.frame(Protein=subtemp$PROTEIN,LogIntensities=subtemp$ABUNDANCE, RUN=subtemp$RUN, NumMeasuredFeature = as.vector(numFea), MissingPercentage=as.vector(numFeaPercentage), more50missing=numFeaTF)
      				
				}

      			result <- rbind(result, sub.result)
      		}
      	}  ## loop for proteins
      	
      	dataafterfit <- NULL
	}
		
	###################################
	## Method 3 : log sum	
	if (summaryMethod=="logOfSum") {
		
		if (label) {
			message("* For log2(sum of intensities) summarization with label-based experiment, ratio between endogenous intensity and reference intensity is used.")
		}
		
		data <- data[!is.na(data$ABUNDANCE),]
    	data$PROTEIN <- factor(data$PROTEIN)
    	data$RUN <- factor(data$RUN)
    
		result <- NULL
		
		single <- unique(data[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
	  	singlesubject <- any(xtabs(~GROUP_ORIGINAL,single)==1)
	  	
		for(i in 1: nlevels(data$PROTEIN)) {
              		
     		sub <- data[data$PROTEIN==levels(data$PROTEIN)[i],]
     		
     		message(paste("Getting the summarization by log2 (sum of intensities) per subject for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

 			## 1.remove runs which have any missing values : same as Skyline Filtering
 			if (filterLogOfSum) {
 				sub$FEATURE <- factor(sub$FEATURE)
				countmissing <- xtabs(~RUN, data=sub) ## already removed NAs
			
				if (label) {
					nomissingrun <- names(countmissing[countmissing==(2*length(levels(sub$FEATURE)))])

				}else{
					nomissingrun <- names(countmissing[countmissing==length(levels(sub$FEATURE))])
				}
			
			###### if there is any run which have any missing values, remove those runs

				sub <- sub[which(sub$RUN %in% nomissingrun),]
			}
			
			if (nrow(sub)==0) {
				
				message(paste("* All replicates in ",levels(data$PROTEIN)[i], " have at least one missing value. Can't summarize by log2 (sum of intensities)."))
				
				next
			}
			
		
			## 2. sum of normalized intensities per run
		
			#### back to original scale from normalized abundance
			sub$INTENSITY <- 2^(sub$ABUNDANCE)
		
			#### 2.1 label based
			if (label) {
				sub.l <- sub[sub$LABEL=="L",]
				sub.h <- sub[sub$LABEL=="H",]
				
				data_w.l = dcast(SUBJECT_ORIGINAL+RUN ~ FEATURE, data=sub.l, value.var='ABUNDANCE', keep=TRUE)

				data_w.h = dcast(SUBJECT_ORIGINAL+RUN ~ FEATURE, data=sub.h, value.var='ABUNDANCE', keep=TRUE)
				
				if (length(unique(sub$FEATURE))==1) {
					
					data_w.l$sumInt <- data_w.l[,3]/data_w.h[,3]
					
					subsum <- data_w.l[,c("SUBJECT_ORIGINAL","RUN","sumInt")]
				
				}else{
					
					data_w <- data_w.l[,c(3:ncol(data_w.l))]/data_w.h[,c(3:ncol(data_w.h))]
					
					data_w <- cbind(SUBJECT_ORIGINAL=data_w.l$SUBJECT_ORIGINAL,RUN=data_w.l$RUN,data_w)
				
					data_w$sumInt <- apply(data_w[,-c(1:2)],1,sum)
				
					subsum <- data_w[,c("SUBJECT_ORIGINAL","RUN","sumInt")]

				}
				
			}else{
			#### 2.2 label-free

				subsum <- aggregate(sub$INTENSITY,list(SUBJECT_ORIGINAL=sub$SUBJECT_ORIGINAL,RUN=sub$RUN),sum, na.rm=TRUE)
				colnames(subsum)[3] <- "sumInt"
			}
			
			## 3. log2
			subsum$logsum <- log2(subsum$sumInt)
		
			if (!singlesubject) {
				
				## this is what Skyline does
				## 4. average per subject
               	sublogsum <- aggregate(subsum$logsum,list(SUBJECT_ORIGINAL=subsum$SUBJECT_ORIGINAL), mean, na.rm=TRUE)
				colnames(sublogsum)[2] <- "LogIntensities"
		
				## 5. make output
				sub.result <- data.frame(Protein=unique(sub$PROTEIN),LogIntensities=sublogsum$LogIntensities, SUBJECT_ORIGINAL=sublogsum$SUBJECT_ORIGINAL)

			}else{
				
				## make output per run
				sub.result <- data.frame(Protein=unique(sub$PROTEIN),LogIntensities=subsum$logsum, RUN=subsum$RUN)

			}
			
      		result <- rbind(result, sub.result)
		} ## end loop for each protein
		
		dataafterfit <- NULL
	}
	
	###################################
	## method 4 : survival model for censored missing values
	if (summaryMethod=="linear" & !is.null(censoredInt)) {
			
		#data <- data[!is.na(data$ABUNDANCE),]
    	data$PROTEIN <- factor(data$PROTEIN)
    	data$RUN <- factor(data$RUN)
    	
		if (label) {
			
			result <- NULL

			for(i in 1:length(unique(data$PROTEIN))) {
	
				sub <- data[data$PROTEIN==unique(data$PROTEIN)[i],]
				
				message(paste("Getting the summarization for censored missing values per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
				
				sub$FEATURE <- factor(sub$FEATURE)
				sub$feature.label <- paste(sub$FEATURE, sub$LABEL, sep="_")
				sub$run.label <- paste(sub$RUN, sub$LABEL, sep="_")

				## if all measurements are NA,
      			if (nrow(sub)==sum(is.na(sub$ABUNDANCE))) {
       				message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),") because all measurements are NAs."))
        			next()
      			}
      			
				## remove run which has no measurement at all
				subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY),]
				count <- aggregate(ABUNDANCE~RUN,data=subtemp, length)
				norun <- setdiff(unique(data$RUN),count$RUN)
				
				if (length(norun)!=0 & length(intersect(norun, as.character(unique(sub$RUN))))) { # removed NA rows already, if there is no overlapped run, error
					sub <- sub[-which(sub$RUN %in% norun),]
					sub$RUN <- factor(sub$RUN)
				}
				
				if (length(unique(sub$RUN))==1) {
				
					message(paste("* Only 1 MS run in ",levels(data$PROTEIN)[i], " has measurement. Can't summarize with censored intensities."))
				
					next
				}	
				
				## remove features which are completely NAs or zero
				subtemp <- sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
				countfeature <- xtabs(~FEATURE, subtemp)
				namefeature <- names(countfeature)[countfeature==0]
				
				if (length(namefeature)!=0) {
					sub <- sub[-which(sub$FEATURE %in% namefeature),]
					sub$FEATURE <- factor(sub$FEATURE)
				}		
				
				##### how to decide censored or not
				## 1. censored 
				if (censoredInt=="0") {
					sub$cen <- ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY==0,0,1)
				}
				
				### 2. all censored missing
				if (censoredInt=="NA") {
					sub$cen <- ifelse(is.na(sub$INTENSITY),0,1)
				}
				
				##### cutoffCensored
				## 1. put minimum in protein level to NA
				#if (cutoffCensored=="minEachProtein") {
				#	if (censoredInt=="NA") {
				#		cut <- min(sub$ABUNDANCE, na.rm=TRUE) 
				#		sub[is.na(sub$INTENSITY),"ABUNDANCE"] <- cut
				#	} 
					
				#	if (censoredInt=="0") {
				#		cut <- min(sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,"ABUNDANCE"]) 
				#		sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut
				#	}
				#}
				
				## 2. put minimum in feature level to NA
				if (cutoffCensored=="minFeature") {
					if (censoredInt=="NA") {
						cut <- aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$feature.label))) {
							sub[is.na(sub$INTENSITY) & sub$feature.label==cut$feature.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
					
					if (censoredInt=="0") {
						subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						cut <- aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$feature.label))) {
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0  & sub$feature.label==cut$feature.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
				}
				
				## 3. put minimum in RUN to NA
				if (cutoffCensored=="minRun") {
					if (censoredInt=="NA") {
						cut <- aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$run.label))) {
							sub[is.na(sub$INTENSITY) & sub$run.label==cut$run.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
					
					if (censoredInt=="0") {
						subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						cut <- aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$run.label))) {
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$run.label==cut$run.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
				}
				
				## 20150829 : 4. put minimum RUN and FEATURE
				if (cutoffCensored=="minFeatureNRun") {
					if (censoredInt=="NA") {
						
						## cutoff for each feature is little less than minimum abundance in a run.

						cut.fea <- aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
						cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
						
						## cutoff for each Run is little less than minimum abundance in a run.

						cut.run <- aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
						cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
						
						
						if (length(unique(sub$feature.label))>1) {
							for(j in 1:length(unique(sub$feature.label))) {
								for(k in 1:length(unique(sub$run.label))) {
									# get smaller value for min Run and min Feature
									finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[is.na(sub$INTENSITY) & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"] <- finalcut
								}
							}
						}
							# if single feature, not impute
					}
					
					if (censoredInt=="0") {
						subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]

						cut.fea <- aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
						cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
												
						## remove runs which has more than 50% missing values
						## before removing, need to contribute min feature calculation
						if (remove50missing) {
							if (length(removerunid)!=0) {
								sub <- sub[-which(sub$RUN %in% removerunid),]
								sub$RUN <- factor(sub$RUN)
							}
						}

						cut.run <- aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
						cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
	
						if (length(unique(sub$feature.label))>1) {
							for(j in 1:length(unique(sub$feature.label))) {
								for(k in 1:length(unique(sub$run.label))) {
									# get smaller value for min Run and min Feature
									finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"] <- finalcut
								}
							}
						}else{ # single feature
					
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut.fea$ABUNDANCE
							
						}
					}
				}	
	
				## when number of measurement is less than df, error for fitting
				subtemp <- sub[!is.na(sub$ABUNDANCE),]
				countdf <- nrow(subtemp)<(length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
				
				
				### fit the model
				if (length(unique(sub$FEATURE))==1) {
					
					# with single feature, not converge, wrong intercept
					# need to check
					fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN+ref,data=sub, dist='gaussian')
				}else{
					if (countdf) {
						fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN+ref,data=sub, dist='gaussian')
					}else{
						fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN+ref,data=sub, dist='gaussian')
					}
				}
				
		
				sub.result <- data.frame(Protein=unique(sub$PROTEIN),RUN=rep(c(levels(sub$RUN)),1),LogIntensities=NA)

				# get the parameters
				cf  <-  summary(fittest)$coefficients

				# calculate sample quantification for all levels of sample
      			a=1	
          
       			for(j in 1:nlevels(sub$RUN)) {
       				contrast.matrix <- rep(0,nlevels(sub$RUN))
        			contrast.matrix[j] <- 1
          					contrast <- .make.contrast.run.quantification.Survival(fittest,contrast.matrix,sub, labeled=TRUE)
	
         			 sub.result[a,3] <- .estimableFixedQuantificationSurvival(cf,contrast)
         			 a=a+1
        		}

				result <- rbind(result, sub.result)
			}

			datamat = dcast( Protein ~ RUN, data=result, value.var='LogIntensities', keep=TRUE)
			datamat = melt(datamat, id.vars=c('Protein'))
			colnames(datamat) <- c('Protein','RUN','LogIntensities')
			result <- datamat
			
		}else{
			
			result <- NULL

			for(i in 1:length(unique(data$PROTEIN))) {
	
				sub <- data[data$PROTEIN==unique(data$PROTEIN)[i],]
				
				message(paste("Getting the summarization for censored missing values per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
				
				sub$FEATURE <- factor(sub$FEATURE)
	
				## if all measurements are NA,
      			if (nrow(sub)==sum(is.na(sub$ABUNDANCE))) {
       				message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),") because all measurements are NAs."))
        			next
      			}
      			
				## remove run which has no measurement at all
				subtemp <- sub[!is.na(sub$INTENSITY),]
				count <- aggregate(ABUNDANCE~RUN,data=subtemp, length)
				norun <- setdiff(unique(data$RUN),count$RUN)
				
				if (length(norun)!=0 & length(intersect(norun, as.character(unique(sub$RUN))))!=0) { # removed NA rows already, if there is no overlapped run, error
					sub <- sub[-which(sub$RUN %in% norun),]
					sub$RUN <- factor(sub$RUN)
				}
				
				if (length(unique(sub$RUN))==1) {
				
					message(paste("* Only 1 MS run in ",levels(data$PROTEIN)[i], " has measurement. Can't summarize with censored intensities."))
				
					next
				}	
						
				
				## remove features which are (completely NAs or zero) 
				subtemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0 ,]
				countfeature <- xtabs(~FEATURE, subtemp)
				namefeature <- names(countfeature)[countfeature==0]
				
				if (length(namefeature)!=0) {
					sub <- sub[-which(sub$FEATURE %in% namefeature),]
					sub$FEATURE <- factor(sub$FEATURE)
				}
				
				if (nrow(sub)==0) {
				
					message(paste("* All measurements are NAs or only one measurement per feature in ",levels(data$PROTEIN)[i], ". Can't summarize with censored intensities."))
				
					next
				}	

				##### how to decide censored or not
				## 1. censored 
				if (censoredInt=="0") {
					sub$cen <- ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY==0,0,1)
				}
				
				### 2. all censored missing
				if (censoredInt=="NA") {
					sub$cen <- ifelse(is.na(sub$INTENSITY),0,1)
				}

				## 3. decide random above some point
	
				## get runs which has all features
				#subtemp <- sub[!is.na(sub$INTENSITY),]
				#count <- aggregate(ABUNDANCE~RUN,data=subtemp, length)
	
				#completerun <- count[count$ABUNDANCE==length(unique(sub$FEATURE)),"RUN"]
	
				#if (length(completerun)!=0) {
				#	subtemp <- sub[which(sub$RUN %in% completerun),]
				#}else{
				#	subtemp <- sub[which(sub$RUN %in% count[which.max(count$ABUNDANCE),"RUN"]),]
				#}
				
				# get feature mean and make order of feature
				# mean or median?
				#featureorder <- aggregate(ABUNDANCE~FEATURE,data=subtemp, mean)
				#featureorder <- featureorder[with(featureorder, order(ABUNDANCE, decreasing=T)),]
	
				# runs which has any missing
				#if (length(completerun)!=0) {
				#	incompleterun <- count[count$ABUNDANCE!=length(unique(sub$FEATURE)),"RUN"]
				#}else{
				#	incompleterun <- count[-which.max(count$ABUNDANCE),"RUN"]
				#}
	
				#if (length(incompleterun)!=0) {
				#	for(j in 1:length(incompleterun)) {
				#		temp <- sub[sub$RUN==incompleterun[j],]
				#		temptemp <- temp[!is.na(temp$INTENSITY),]
		
				#		minfeature <- temptemp[which.min(temptemp$ABUNDANCE),"FEATURE"]
				#		abovefeature <- featureorder[1:which(featureorder$FEATURE==minfeature),"FEATURE"]
	
				#		sub[which(sub$RUN==incompleterun[j] & sub$FEATURE %in% abovefeature & is.na(sub$INTENSITY)),"ABUNDANCE"] <- NA
				#		sub[which(sub$RUN==incompleterun[j] & sub$FEATURE %in% abovefeature & is.na(sub$INTENSITY)),"cen"] <- 1
				#	}
				#}
				
				##### cutoffCensored
				## 1. put minimum in protein level to NA
				#if (cutoffCensored=="minEachProtein") {
				#	if (censoredInt=="NA") {
				#		cut <- min(sub$ABUNDANCE, na.rm=TRUE) 
				#		sub[is.na(sub$INTENSITY),"ABUNDANCE"] <- cut
				#	} 
					
				#	if (censoredInt=="0") {
				#		cut <- min(sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,"ABUNDANCE"]) 
				#		sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut
				#	}
				#}
				
				## 2. put minimum in feature level to NA
				if (cutoffCensored=="minFeature") {
					if (censoredInt=="NA") {
						cut <- aggregate(ABUNDANCE~FEATURE,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$FEATURE))) {
							sub[is.na(sub$INTENSITY) & sub$FEATURE==cut$FEATURE[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
					
					if (censoredInt=="0") {
						subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0 ,]
						cut <- aggregate(ABUNDANCE~FEATURE,data=subtemptemp, FUN=min)
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$FEATURE))) {
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$FEATURE==cut$FEATURE[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
				}
				
				## 3. put minimum in RUN to NA
				if (cutoffCensored=="minRun") {
					if (censoredInt=="NA") {
						cut <- aggregate(ABUNDANCE~RUN,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$RUN))) {
							sub[is.na(sub$INTENSITY) & sub$RUN==cut$RUN[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
					
					if (censoredInt=="0") {
						subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						cut <- aggregate(ABUNDANCE~RUN,data=subtemptemp, FUN=min)
						
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE <- 0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$RUN))) {
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$RUN==cut$RUN[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
						}
					}
				}	
				
				## 20150829 : 4. put minimum RUN and FEATURE
				if (cutoffCensored=="minFeatureNRun") {
					if (censoredInt=="NA") {
						
						## cutoff for each feature is little less than minimum abundance in a run.
						cut.fea <- aggregate(ABUNDANCE~FEATURE,data=sub, function(x) min(x, na.rm=TRUE))
						cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
												
						## cutoff for each Run is little less than minimum abundance in a run.

						cut.run <- aggregate(ABUNDANCE~RUN,data=sub, function(x) min(x, na.rm=TRUE))
						cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
						
						if (length(unique(sub$FEATURE))>1) {
							for(j in 1:length(unique(sub$FEATURE))) {
								for(k in 1:length(unique(sub$RUN))) {
									# get smaller value for min Run and min Feature
									finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[is.na(sub$INTENSITY) & sub$FEATURE==cut.fea$FEATURE[j] & sub$RUN==cut.run$RUN[k],"ABUNDANCE"] <- finalcut
								}
							}
						}
							# if single feature, not impute
					}
					
					if (censoredInt=="0") {
						subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]

						cut.fea <- aggregate(ABUNDANCE~FEATURE,data=subtemptemp, FUN=min)
						cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
						
						cut.run <- aggregate(ABUNDANCE~RUN,data=subtemptemp, FUN=min)
						cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
	
						if (length(unique(sub$FEATURE))>1) {
							for(j in 1:length(unique(sub$FEATURE))) {
								for(k in 1:length(unique(sub$RUN))) {
									# get smaller value for min Run and min Feature
									finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$FEATURE==cut.fea$FEATURE[j] & sub$RUN==cut.run$RUN[k],"ABUNDANCE"] <- finalcut
								}
							}
						}else{ # single feature
					
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut.fea$ABUNDANCE
							
						}
					}
				}
					
				
				## when number of measurement is less than df, error for fitting
				subtemp <- sub[!is.na(sub$ABUNDANCE),]
				countdf <- nrow(subtemp)<(length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
				
				### fit the model
				if (length(unique(sub$FEATURE))==1) {
					fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
				}else{
					if (countdf) {
						fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
					}else{
						fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN,data=sub, dist='gaussian')
					}
				}
				
		
				sub.result <- data.frame(Protein=unique(sub$PROTEIN),RUN=rep(c(levels(sub$RUN)),1),LogIntensities=NA)

				# get the parameters
				cf  <-  summary(fittest)$coefficients

				# calculate sample quantification for all levels of sample
      			a=1	
          
       			for(j in 1:nlevels(sub$RUN)) {
       				contrast.matrix <- rep(0,nlevels(sub$RUN))
        			contrast.matrix[j] <- 1
          					contrast <- .make.contrast.run.quantification.Survival(fittest,contrast.matrix,sub, labeled=FALSE)
	
         			 sub.result[a,3] <- .estimableFixedQuantificationSurvival(cf,contrast)
         			 a=a+1
        		}

				result <- rbind(result, sub.result)
			}

			datamat = dcast( Protein ~ RUN, data=result, value.var='LogIntensities', keep=TRUE)
			datamat = melt(datamat, id.vars=c('Protein'))
			colnames(datamat) <- c('Protein','RUN','LogIntensities')
			result <- datamat
		}
		dataafterfit <- NULL	
	}
	
	###################################
	## final result
	finalout <- list(rqdata=result,ModelQC=dataafterfit, PredictedBySurvival=predAbundance)
	return(finalout)
}



##########################################################################################
## updated v3
.fit.quantification.run <- function(sub,singleFeature,singleSubject, TechReplicate,labeled,equalFeatureVar) {
	
	if (!labeled) { ### label-free case
		## for single Feature, original value is the run quantification
		if (singleFeature) {
			fit.full <- lm(ABUNDANCE ~ RUN , data = sub)
		}else{
			fit.full <- lm(ABUNDANCE ~ FEATURE + RUN , data = sub)
		}
		
	}else{ ### labeled-based case
		### update v3
		if (singleFeature) {
			fit.full <- lm(ABUNDANCE ~ RUN+ref , data = sub)
		}else{ ## several subjects
			fit.full <- lm(ABUNDANCE ~ FEATURE+RUN+ref , data = sub)
		}
	}
	
	## make equal variance for feature : need to be updated
    if (!equalFeatureVar) {
       fit.full <- .iter.wls.fit.model(data=sub,fit=fit.full,nrepeats=1)
    }
	
	return(fit.full)
}




#############################################
# check whether there are multiple runs for a replicate
# if yes, normalization should be different way.
#############################################

.countMultiRun <- function(data) {
  
  standardFeature <- unique(data[data$RUN == unique(data$RUN[1]), "FEATURE"]) ## if some feature are missing for this spedific run, it could be error. that is why we need balanced design.
  
  ## get overlapped feature ID
  countdiff = tapply (data$FEATURE, data$RUN, function ( x ) length(intersect(unique(x), standardFeature)) ) 
  
  return(countdiff)
}



.Imputation <- function(data, cutoffCensored, censoredInt, MBimpute, remove50missing) {
	
  #ABUNDANCE.O preserves the intensity before imputations so that this can be used for drawing the profile plot later
  data$ABUNDANCE.O <- data$ABUNDANCE
  data$LABEL <- factor(data$LABEL)
  label <- nlevels(data$LABEL) == 2
    
  # set ref which is distinguish reference and endogenous. any reference=0. endogenous is the same as RUN
	if (label) {
    data$ref <- 0
		data$ref[data$LABEL!="H"] <- data$RUN[data$LABEL!="H"]
		data$ref <- factor(data$ref)
#		unique(data[,c("RUN","LABEL","GROUP","ref")])
	}
	      
#    finalresult <- data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$RUN)),RUN=rep(c(levels(data$RUN)),nlevels(data$PROTEIN)),Condition=NA, BioReplicate=NA,LogIntensities=NA,NumFeature=NA,NumPeaks=NA)

	# for saving predicting value for impute option
  predAbundance <- NULL
	AbundanceAfterImpute <- NULL
		
	###################################		
		#data <- data[!is.na(data$ABUNDANCE),]
  data$PROTEIN <- factor(data$PROTEIN)
  data$RUN <- factor(data$RUN)
    
  result <- NULL
	  
  for(i in 1: nlevels(data$PROTEIN)) {
              		
    sub <- data[data$PROTEIN==levels(data$PROTEIN)[i],]
     		
    message(paste("Imputing the censoring intensities for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

    sub$FEATURE <- factor(sub$FEATURE)	
    sub$feature.label <- paste(sub$FEATURE, sub$LABEL, sep="_")
    sub$run.label <- paste(sub$RUN, sub$LABEL, sep="_")
      		
    ## if all measurements are NA,
    if (nrow(sub)==sum(is.na(sub$ABUNDANCE))) {
      message(paste("Can't impute for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),") because all measurements are NAs."))
      next()
    }
      		
    ## remove features which are completely NAs
    subtemp <- sub[sub$LABEL == "L" & !is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
    countfeature <- xtabs(~FEATURE, subtemp)
    namefeature <- names(countfeature)[countfeature == 0]
				
    if (length(namefeature) != 0) {
      sub <- sub[-which(sub$FEATURE %in% namefeature), ]
				
      if (nrow(sub) == 0) {
        message(paste("Can't impute for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(data$PROTEIN)), ") because all measurements are NAs."))
        next()
        		
      } else {
        sub$FEATURE <- factor(sub$FEATURE)
      }
    }
			
    ## remove features which have only 1 measurement.
    namefeature1 <- names(countfeature)[countfeature == 1]
				
    if (length(namefeature1)!=0) {
      sub <- sub[-which(sub$FEATURE %in% namefeature1), ]
				 
      if (nrow(sub) == 0) {
        message(paste("Can't impute for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(data$PROTEIN)), ") because features have only one measurement across MS runs."))
        next()
        		
      } else {
        sub$FEATURE <- factor(sub$FEATURE)
			}
		}
    
    ## 20160425-MC : If a single feature, pass the imputation.
    ## If impute for this NA, extreme values can be imputed, zero or inf
    subtemp <- sub[!is.na(sub$ABUNDANCE), ]
		sigFeature <- .checkSingleFeature(subtemp)
    rm(subtemp)
    
    if (sigFeature) {
      ## no imputation
      sub$cen <- NA
      sub$pred <- NA
      AbundanceAfterImpute <- rbind(AbundanceAfterImpute, sub)
      
    } else {
      ## how to decide censored or not
      if (!is.null(censoredInt)) {
        ## 1. censored 
        if (censoredInt=="0") {
          sub$cen <- ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY==0,0,1)
        }
        
        ### 2. all censored missing
        if (censoredInt=="NA") {
          sub$cen <- ifelse(is.na(sub$INTENSITY),0,1)
        }
        
        ### check whether we need to impute or not.
        if (sum(sub$cen==0)>0) {
          
          ## 2. put minimum in feature level to NA
          if (cutoffCensored=="minFeature") {
            if (censoredInt=="NA") {
              cut <- aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
              ## cutoff for each feature is little less than minimum abundance in a run.
              cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
              
              ## remove runs which has more than 50% missing values
              if (remove50missing) {
                if (length(removerunid)!=0) {
                  sub <- sub[-which(sub$RUN %in% removerunid),]
                  sub$RUN <- factor(sub$RUN)
                }
              }
              
              for(j in 1:length(unique(cut$feature.label))) {
                sub[is.na(sub$INTENSITY) & sub$feature.label==cut$feature.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
              }
            }
            
            if (censoredInt=="0") {
              subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
              cut <- aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
              ## cutoff for each feature is little less than minimum abundance in a run.
              cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
              
              ## remove runs which has more than 50% missing values
              if (remove50missing) {
                if (length(removerunid)!=0) {
                  sub <- sub[-which(sub$RUN %in% removerunid),]
                  sub$RUN <- factor(sub$RUN)
                }
              }
              
              for(j in 1:length(unique(cut$feature.label))) {
                sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0  & sub$feature.label==cut$feature.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
              }
            }
          }
          
          ## 3. put minimum in RUN to NA
          if (cutoffCensored=="minRun") {
            
            ## remove runs which has more than 50% missing values
            if (remove50missing) {
              if (length(removerunid)!=0) {
                sub <- sub[-which(sub$RUN %in% removerunid),]
                sub$RUN <- factor(sub$RUN)
              }
            }
            
            if (censoredInt=="NA") {
              cut <- aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
              ## cutoff for each Run is little less than minimum abundance in a run.
              cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
              
              for(j in 1:length(unique(cut$run.label))) {
                sub[is.na(sub$INTENSITY) & sub$run.label==cut$run.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
              }
            }
            
            if (censoredInt=="0") {
              subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
              cut <- aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
              cut$ABUNDANCE <- 0.99*cut$ABUNDANCE
              
              for(j in 1:length(unique(cut$run.label))) {
                sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$run.label==cut$run.label[j],"ABUNDANCE"] <- cut$ABUNDANCE[j]
              }
            }
          }
          
          ## 20150829 : 4. put minimum RUN and FEATURE
          if (cutoffCensored=="minFeatureNRun") {
            if (censoredInt=="NA") {
              
              ## cutoff for each feature is little less than minimum abundance in a run.
              
              cut.fea <- aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
              cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
              
              ## remove runs which has more than 50% missing values
              ## before removing, need to contribute min feature calculation
              if (remove50missing) {
                if (length(removerunid)!=0) {
                  sub <- sub[-which(sub$RUN %in% removerunid),]
                  sub$RUN <- factor(sub$RUN)
                }
              }
              
              ## cutoff for each Run is little less than minimum abundance in a run.
              
              cut.run <- aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
              cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
              
              if (length(unique(cut.fea$feature.label))>1) {
                for(j in 1:length(unique(cut.fea$feature.label))) {
                  for(k in 1:length(unique(cut.run$run.label))) {
                    # get smaller value for min Run and min Feature
                    finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
                    
                    sub[is.na(sub$INTENSITY) & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"] <- finalcut
                  }
                }
              }
              # if single feature, not impute
            }
            
            if (censoredInt=="0") {
              subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
              
              cut.fea <- aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
              cut.fea$ABUNDANCE <- 0.99*cut.fea$ABUNDANCE
              
              ## remove runs which has more than 50% missing values
              ## before removing, need to contribute min feature calculation
              if (remove50missing) {
                if (length(removerunid)!=0) {
                  sub <- sub[-which(sub$RUN %in% removerunid),]
                  sub$RUN <- factor(sub$RUN)
                }
              }
              
              cut.run <- aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
              cut.run$ABUNDANCE <- 0.99*cut.run$ABUNDANCE
              
              if (length(unique(cut.fea$feature.label))>1) {
                for(j in 1:length(unique(cut.fea$feature.label))) {
                  for(k in 1:length(unique(cut.run$run.label))) {
                    # get smaller value for min Run and min Feature
                    finalcut <- min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
                    
                    sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"] <- finalcut
                  }
                }
              }else{ # single feature
                
                sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut.fea$ABUNDANCE
                
              }
            }
          }
          
          if (MBimpute) {
            
            if (nrow(sub[sub$cen==0,])>0) {
              ## impute by survival model
              subtemp <- sub[!is.na(sub$ABUNDANCE),]
              countdf <- nrow(subtemp)<(length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
              
              ## fit the model
              if (length(unique(sub$FEATURE))==1) {
                fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN, data=sub, dist='gaussian')
              } else {
                if (countdf) {
                  fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN, data=sub, dist='gaussian')
                } else {
                  fittest <- survreg(Surv(ABUNDANCE, cen, type='left') ~ FEATURE + RUN, data=sub, dist='gaussian')
                }
              }
              
              # get predicted value from survival
              sub <- data.frame(sub, pred=predict(fittest, newdata=sub, type="response"))
              
              # the replace censored value with predicted value
              sub[sub$cen==0,"ABUNDANCE"] <- sub[sub$cen==0,"pred"]	
              
              # save predicted value
              #predAbundance <- c(predAbundance,predict(fittest, newdata=sub, type="response"))
              AbundanceAfterImpute <- rbind(AbundanceAfterImpute, sub)				
            } 
          }	
        } else {  
          sub$pred <- NA; 
          AbundanceAfterImpute <- rbind(AbundanceAfterImpute, sub)
        }
      } ## end for censoredInt is not null
    } ## end for multiple feature, imputation  			
  }  ## loop for proteins
      		
	###################################
	## Output: data frame after imputation
	return(AbundanceAfterImpute)
} #End of function .Imputation


## to get feature ID which has completely missing in a certain condition, and not from the protein which has only one feature
.getfeatureID <- function(data){
  countfeature <- xtabs( ~ FEATURE + GROUP_ORIGINAL, data[!is.na(data$ABUNDANCE), ])

  getfeature.missingcondition <- apply(countfeature, 1, function(x) any(x==0))
  ## feature ID, which have any completely missing in a condition
  getissuefeature <- names(getfeature.missingcondition[getfeature.missingcondition])

  getproteinfeature <- unique(data[, c("PROTEIN", "FEATURE")])
  countfeature.perpro <- xtabs( ~ PROTEIN, getproteinfeature)
  ## feature ID, which comes from the proteins which have only one feature
  getsinglefeature <- getproteinfeature[which(getproteinfeature$PROTEIN %in% names(countfeature.perpro[countfeature.perpro == 1])), "FEATURE"]
  getissuefeature <- getissuefeature[-which(getissuefeature %in% getsinglefeature)]
  
  return(getissuefeature)
}

