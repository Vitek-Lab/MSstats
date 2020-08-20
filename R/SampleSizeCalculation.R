
#############################################
## designSampleSize
#############################################

#' @export
designSampleSize <- function(data=data,
                             desiredFC=desiredFC,
                             FDR=0.05,
                             numSample=TRUE,
                             power=0.9) {
  
  	labeled <- FALSE
  	scopeOfBioReplication <- "expanded"
  	interference <- TRUE
  	equalFeatureVar <- TRUE
  
  	## save process output in each step
  	allfiles <- list.files()
  	filenaming <- "msstats"
  
  	if (length(grep(filenaming, allfiles)) == 0) {
    
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
    	processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
  	}
  
  	processout <- rbind(processout, as.matrix(c(" ", " ", "MSstats - designSampleSize function", " "), ncol=1))
  
  	processout <- rbind(processout, c(paste0("Desired fold change = ", paste(desiredFC, collapse=" - "))))
  	processout <- rbind(processout, c(paste0("FDR = ", FDR)))
  	processout <- rbind(processout, c(paste0("Power = ", power)))
  
  	write.table(processout, file=finalfile, row.names=FALSE)
  
  	## for label-free experiment
  	#if (!labeled) {
    
    sigma.error <- NULL	
    VarComponent <- data.frame("Protein"=seq(1, length(data)), 
                               "Error"=NA, 
                               "Subject"=NA, 
                               "GroupBySubject"=NA)
    
    for (i in 1:length(data)) {
      
      	# note: when run is fixed, we can obtain the same variance of error for both case-control and time course studies.
      
      	fit.full <- data[[i]]
      
      	## if fit.full==NA (class(fit.full)=="try-error)
      	if (is.null(fit.full)) {
        	## !!!!! but if we have NULL for last protein?
        	next
        
      	} else {
        
        	## get variance component
        
        	if (class(fit.full) != "lmerMod") {
          		VarComponent[i, "Error"] <- summary(fit.full)$sigma^2
        	} else {
          		stddev  <-  c(sapply(VarCorr(fit.full), function(el) attr(el, "stddev")), attr(VarCorr(fit.full), "sc"))
         	 	VarComponent[i, "Error"] <- stddev[names(stddev) == ""]^2
         	 	
          		if (sum(names(stddev) %in% "SUBJECT_NESTED.(Intercept)") > 0) {
            		VarComponent[i, "Subject"] <- stddev[names(stddev) == "SUBJECT_NESTED.(Intercept)"]^2
          		}
          		if (sum(names(stddev) %in% "SUBJECT.(Intercept)") > 0) {
            		VarComponent[i, "Subject"] <- stddev[names(stddev) == "SUBJECT.(Intercept)"]^2
          		}
          		if (sum(names(stddev) %in% "SUBJECT:GROUP.(Intercept)") > 0) {
            		VarComponent[i, "GroupBySubject"] <- stddev[names(stddev) == "SUBJECT:GROUP.(Intercept)"]^2
          		}
        	}
      	}  
    } ## end-loop
    
    ## VarComponent[is.na(VarComponent)] <- 0	
    ## for label-free DDA, there are lots of missingness and lots of zero SE. So, remove NA SE.
    median.sigma.error <- median(VarComponent[, "Error"], na.rm=TRUE)
    
   	if (sum(!is.na(VarComponent[, "GroupBySubject"])) > 0) {
    	
    	median.sigma.subject <- median(VarComponent[, "GroupBySubject"], na.rm=TRUE)
      
    } else {
      	if (sum(!is.na(VarComponent[, "Subject"])) > 0) {
        	median.sigma.subject <- median(VarComponent[, "Subject"], na.rm=TRUE)
     	} else {
        	median.sigma.subject <- 0
      	}
    }
    
    ##
    processout <- rbind(processout, "Calculated variance component. - okay")
    write.table(processout, file=finalfile, row.names=FALSE)
        
    ## power calculation
    if (isTRUE(power)) {
      	delta <- log2(seq(desiredFC[1], desiredFC[2], 0.025))
      	desiredFC <- 2^delta
      	m0_m1 <- 99
     	t <- delta / sqrt(2 * (median.sigma.error / numSample + median.sigma.subject / numSample))
      	#t <- delta/sqrt(2*median.sigma.error/numPep/numTran/numSample+median.sigma.subject/numSample)
     	powerTemp <- seq(0, 1, 0.01)
      
      	power <- NULL
      	for (i in 1:length(t)) {
       	 	diff <- qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
        	min(abs(diff), na.rm=TRUE)
        	power[i] <- powerTemp[order(abs(diff))][1]
      	}
      
      	CV <- round( (2 * (median.sigma.error / numSample + median.sigma.subject / numSample)) / desiredFC, 3)
     
      	###
     	processout <- rbind(processout, 
     	                    "Power is calculated. - okay")
      	write.table(processout, file=finalfile, row.names=FALSE)
      
      	out <- data.frame(desiredFC, numSample, FDR, power=power, CV)

      	return(out)
    }	
    
	if (is.numeric(power)) {
      
      	## Large portion of proteins are not changing
      	m0_m1 <- 99 ## it means m0/m1=99, m0/(m0+m1)=0.99
      	alpha <- power * FDR / (1 + (1 - FDR) * m0_m1)
      
      	## Num Sample calculation
      	if (isTRUE(numSample)) {
        	delta <- log2(seq(desiredFC[1], desiredFC[2], 0.025))
        	desiredFC <- 2^delta
        	z_alpha <- qnorm(1-alpha/2)
        	z_beta <- qnorm(power)
        	aa <- (delta/(z_alpha+z_beta))^2
        	numSample <- round(2 * (median.sigma.error + median.sigma.subject) / aa, 0)
        	CV <- round(2 * (median.sigma.error / numSample + median.sigma.subject / numSample) / desiredFC, 3)

	  		##
        	processout <- rbind(processout, 
        	                    "The number of sample is calculated. - okay")
        	write.table(processout, file=finalfile, row.names=FALSE)
        
        	out <- data.frame(desiredFC,numSample,FDR,power,CV)
       	 	return(out)
      	}
	} # when power is numeric
  	#} ## label-free
}


#############################################
## designSampleSizePlots
#############################################
#' @export
designSampleSizePlots <- function(data=data) {
  
  	if (length(unique(data$numSample)) > 1) {
  	    index <- "numSample"
  	}
  	if (length(unique(data$power)) > 1) {
  	    index <- "power"
  	}
  
  	if (length(unique(data$numSample)) == 1 & length(unique(data$power)) == 1) {
  	    index <- "numSample"
  	}
  
  	text.size <- 1.2	
  	axis.size <- 1.3	
  	lab.size <- 1.7
  
  	if (index == "numSample") {
    
        plot(data$desiredFC, data$numSample, 
             lwd=2, xlab="", ylab="", 
             cex.axis=axis.size, type="l", xaxt="n")
  	    axis(1, at=seq(min(data$desiredFC), max(data$desiredFC), 0.05), 
  	         labels=seq(min(data$desiredFC), max(data$desiredFC), 0.05),
  	         cex.axis=axis.size)
  	    #axis(3, at=seq(min(data$desiredFC), max(data$desiredFC), 0.05),
  	    #     labels=data$CV[which(data$desiredFC %in% seq(min(data$desiredFC), max(data$desiredFC), 0.05))],
  	    #     cex.axis=axis.size)
  	    #mtext("Coefficient of variation, CV", 3, line=2.5, cex=lab.size)
  	    mtext("Desired fold change", 1, line=3.5, cex=lab.size)
  	    mtext("Minimal number of biological replicates", 2, line=2.5, cex=lab.size)
  	    legend("topright",
  	           c(paste("FDR is", unique(data$FDR)),
  	             paste("Statistical power is", unique(data$power))),
  	           bty="n",
  	           cex=text.size)
  	}
  
    if (index == "power") {
    
    	plot(data$desiredFC, data$power, 
    	     lwd=2, xlab="", ylab="", 
    	     cex.axis=axis.size, type="l", xaxt="n")
        axis(1, at=seq(min(data$desiredFC), max(data$desiredFC), 0.05), 
             labels=seq(min(data$desiredFC), max(data$desiredFC), 0.05), 
             cex.axis=axis.size)
        mtext("Desired fold change", 1, line=3.5, cex=lab.size)
        mtext("Power", 2, line=2.5, cex=lab.size)
        legend("bottomright", 
               c(paste("Number of replicates is", unique(data$numSample)), 
                 paste("FDR is", unique(data$FDR))), 
               bty="n", 
               cex=text.size)
  	}
}



