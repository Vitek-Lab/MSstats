
#############################################
##  quantification
#############################################


quantification <- function(data,
							type="Sample",
							format="matrix") {
  
  	## save process output in each step
  	allfiles <- list.files()
  	filenaming <- "msstats"
  
  	if (length(grep(filenaming,allfiles)) == 0) {
    
    	finalfile <- "msstats.log"
    	processout <- NULL
    
  	} else {
    
    	num <- 0
    	finalfile <- "msstats.log"
    
    	while(is.element(finalfile, allfiles)) {
      		num <- num + 1
      		lastfilename <- finalfile ## in order to rea
      		finalfile <- paste(paste(filenaming, num, sep="-"), ".log", sep="")
    	}
    
    	finalfile <- lastfilename
    	processout <- as.matrix(read.table(finalfile, header=T,sep="\t"))
  	}
  
  	processout <- rbind(processout, as.matrix(c(" ", " ", "MSstats - quantification function", " "), ncol=1))
  
  
  	## data format
  	#rawinput <- c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
  
  	#if (length(setdiff(toupper(rawinput),toupper(colnames(data))))==0) {
    	#processout <- rbind(processout, c(paste("The required input - data : did not process from dataProcess function. - stop")))
    
    	#write.table(processout, file=finalfile, row.names=FALSE)
    
    	#stop("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison.")
  	#}
  
  	## check other options
  	if (!(toupper(type) == "SAMPLE" | toupper(type) == "GROUP")) {
    	processout <- rbind(processout, c(paste("The required input - type : 'type' value is wrong. - stop")))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop("'type' must be one of \"Sample\" or \"Group\". ")
  	} 
  
  	if (!(toupper(format) == "MATRIX" | toupper(format) == "LONG")) {
    	processout <- rbind(processout, c(paste("The required input - format : 'format' value is wrong. - stop")))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop("'format' must be one of \"matrix\" or \"long\". ")
  	} 
  
  	## all input
  	processout <- rbind(processout, c(paste("type of quantification = ", type, sep="")))
  	processout <- rbind(processout, c(paste("format of output = ", format, sep="")))
  
  	write.table(processout, file=finalfile, row.names=FALSE)
  
  
  	#################################	
  	## sample quantification	
  	#################################	
	
  	if (toupper(type) == "SAMPLE") {
    	datarun <- data$RunlevelData
    	datarun$Protein <- factor(datarun$Protein)
    	
   		datarun <- datarun[!is.na(datarun$LogIntensities), ]
       
   	 	datam <- dcast(Protein ~ GROUP_ORIGINAL + SUBJECT_ORIGINAL, data=datarun, value.var='LogIntensities', fun.aggregate=median)
   	 	
   	 	if (format=="long") {
   	 		data_l <- melt(datam, id.vars=c('Protein'))
			colnames(data_l)[colnames(data_l) %in% c("variable", "value")] <- c('Group_Subject', 'LogIntensity')

   	 	}
  
    	processout <- rbind(processout,c("Finish sample quantificiation - okay."))
    	write.table(processout, file=finalfile, row.names=FALSE)    
    
    	if (format == "long") {
    		return(data_l)
    	}
    	
    	if (format == "matrix") {
    		return(datam)
    	}
    	
  	} ## end sample quantification	
  	
  	#################################
  	## Group quantification
    #################################

  	if (toupper(type)=="GROUP") {
    
    	datarun <- data$RunlevelData
    	datarun$Protein <- factor(datarun$Protein)
    	
   		datarun <- datarun[!is.na(datarun$LogIntensities), ]
       
   	 	datam <- dcast(Protein + GROUP_ORIGINAL ~ SUBJECT_ORIGINAL, data=datarun, value.var='LogIntensities', fun.aggregate=median)

		datam2  <- melt(datam, id.vars=c('Protein', "GROUP_ORIGINAL"))
		colnames(datam2)[colnames(datam2) %in% c("variable", "value")] <- c('Subject', 'LogIntensity')

		datam3 <- dcast(Protein ~ GROUP_ORIGINAL , data=datam2, value.var='LogIntensity', fun.aggregate=median)
		
		rm(datam)
		rm(datam2)
   	 	
   	 	if (format == "long") {
   	 		data_l <- melt(datam3, id.vars=c('Protein'))
			colnames(data_l)[colnames(data_l) %in% c("variable", "value")] <- c('Group', 'LogIntensity')

   	 	}
  
    	processout <- rbind(processout,c("Finish group quantificiation - okay."))
    	write.table(processout, file=finalfile, row.names=FALSE)    
    
    	if (format == "long") {
    		return(data_l)
    	}
    	
    	if (format == "matrix") {
    		return(datam3)
    	}
  	} ## end group quantification	
	
}
