
#############################################
# Whole plot testing
#############################################

#' @export groupComparison
#' @import lme4 
#' @import limma
#' @importFrom doSNOW registerDoSNOW
#' @importFrom snow makeCluster
#' @importFrom data.table rbindlist

groupComparison <- function(contrast.matrix=contrast.matrix,
                            data=data) {
  
    scopeOfBioReplication <- "expanded"
    scopeOfTechReplication <- "restricted"
    message.show <- FALSE ## whether show the 'testing xx protein' or not. It consume large memory and file size.

    ## save process output in each step
    allfiles <- list.files()
    filenaming <- "msstats"
    
    if (length(grep(filenaming, allfiles)) == 0) {
    
        finalfile <- "msstats.log"
        
        session <- sessionInfo()
        sink("sessionInfo.txt")
        print(session)
        sink()
        
        processout <- as.matrix(read.table("sessionInfo.txt", header=TRUE, sep="\t"))
    
    } else {
    
        num <- 0
        finalfile <- "msstats.log"
    
        while (is.element(finalfile, allfiles)) {
            num <- num + 1
            lastfilename <- finalfile ## in order to rea
            finalfile <- paste0(paste(filenaming, num, sep="-"), ".log")
        }
    
        finalfile <- lastfilename
        processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
    }
    
    processout <- rbind(processout,
                        as.matrix(c(" ", " ", "MSstats - groupComparison function", " "), ncol=1))
    
    ## check input is correct
    ## data format
    rawinput <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", 
                  "ProductCharge", "IsotopeLabelType", "Condition", 
                  "BioReplicate", "Run", "Intensity")
    
    if (length(setdiff(toupper(rawinput),toupper(colnames(data$ProcessedData)))) == 0) {
    
        processout <- rbind(processout,
                            "The required input - data : did not process from dataProcess function. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        	
        stop("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison.")
    }
    
    ## contrast. matrix
    if (ncol(contrast.matrix) != length(unique(data$ProcessedData$GROUP_ORIGINAL))) {
        processout <- rbind(processout,
                            "The required input - contrast.matrix: the number of column and the number of group are not the same. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please check contrast matrix. The number of group in data set is different with columns of contrast.matrix.")
    }
    
    ## check whether row.names of contrast.matrix.sub exists or not
    if (sum(is.null(row.names(contrast.matrix))) > 0) {
        processout <- rbind(processout,
                            "The required input - contrast.matrix: need row names of contrast.matrix . - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("No row.names of comparison exist.\n")
    }
    
    
    if (!(scopeOfBioReplication == "restricted" | scopeOfBioReplication == "expanded")) {
        processout <- rbind(processout,
                            "The required input - scopeOfBioReplication : 'scopeOfBioReplication' value is wrong. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("'scopeOfBioReplication' must be one of \"restricted\" or \"expanded\".")
    } 
    
    labeled <- ifelse(length(unique(data$ProcessedData$LABEL)) == 1, FALSE, TRUE)
    
    ## all input
    processout <- rbind(processout,
                        paste0("labeled = ", labeled))
    processout <- rbind(processout,
                        paste0("scopeOfBioReplication = ", scopeOfBioReplication))
    
    write.table(processout, file=finalfile, row.names=FALSE)
    
    ## check whether case-control(FALSE) or time-course(TRUE)
    repeated <- .checkRepeated(data$ProcessedData)
    
    if (repeated) { 
        processout <- rbind(processout, 
                            "Time course design of experiment - okay")
    } else {
        processout <- rbind(processout, 
                            "Case control design of experiment - okay")
    }
    write.table(processout, file=finalfile, row.names=FALSE)
    
    ## for final result report
    out <- NULL
    outsummary <- NULL
    outfitted <- NULL
    dataafterfit <- NULL
      
    ## check input for labeled
    
    ## start to analyze by protein ID
    message("\n Start to test and get inference in whole plot...")
     
    if (data$SummaryMethod == "logOfSum") {
        message("\n ** Use t-test for log sum summarization per subject.")
     	  processout <- rbind(processout, 
     	                    "** Use t-test for log sum summarization per subject.")
     	  write.table(processout, file=finalfile, row.names=FALSE)
    
        if (repeated) {
     		
         	processout <- rbind(processout,
         	                    "Only group comparsion is available for t-test.- stop")
         	write.table(processout, file=finalfile, row.names=FALSE)
        
         	stop("\n * Only group comparsion is available for t-test.")
     	  }
    }
    
    ## need original group information
    rqall <- data$RunlevelData
    rqall$Protein <- factor(rqall$Protein)
    processall <- data$ProcessedData
    
    origGroup <- unique(rqall$GROUP_ORIGINAL)
    groupinfo <- levels(data$ProcessedData$GROUP_ORIGINAL)
      
    pb <- txtProgressBar(max = nlevels(rqall$Protein), style = 3)
      
    for (i in 1:nlevels(rqall$Protein)) {
          
        sub <- rqall[rqall$Protein == levels(rqall$Protein)[i], ]
        colnames(sub)[colnames(sub) == "LogIntensities"] <- "ABUNDANCE"
        colnames(sub)[colnames(sub) == "Protein"] <- "PROTEIN"
    
        # it is important to remove NA first, before we have the correct structure for each factor
        sub <- sub[!is.na(sub$ABUNDANCE), ]
    
        ## 1. for logsum, t-test
        if (data$SummaryMethod=="logOfSum") {
    	
            ## subject level
        	sub$GROUP <- factor(sub$GROUP)
           	sub$SUBJECT <- factor(sub$SUBJECT)
        	sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
    
        	## testing and inference in whole plot
        	if (message.show) {
        		message(paste("Testing a comparison for protein ", 
        		              unique(sub$PROTEIN), "(", i, " of ", length(unique(rqall$Protein)), ")"))
        	}
    		
     	    temp <- try(.ttest.logsum(contrast.matrix, sub, origGroup), silent=TRUE)
     	
        } else { ## linear model
     	
        	sub$GROUP <- factor(sub$GROUP)
          sub$SUBJECT <- factor(sub$SUBJECT)
        	sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL, levels=groupinfo)	
        	sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)
        	sub$SUBJECT_NESTED <- factor(sub$SUBJECT_NESTED)
        	sub$RUN <- factor(sub$RUN)
        
        	singleSubject <- .checkSingleSubject(sub)
        	TechReplicate <- .checkTechReplicate(sub) ## use for label-free model
        
        	MissSubjectByGroup <- .checkRunbyFeature(sub)
        	UnequalSubject <- .checkUnequalSubject(sub)
                
        	## testing and inference in whole plot
        	if (message.show) {
        	  message(paste("Testing a comparison for protein ", 
        	                unique(sub$PROTEIN), "(", i, " of ", length(unique(rqall$Protein)), ")"))
        	}
        		
         	## fit the model 
         	temp <- try(.fit.model.single(contrast.matrix, sub, 
         	                              TechReplicate, singleSubject, repeated, origGroup, processout), 
         	            silent=TRUE)
        }
    
        if (class(temp) == "try-error") {
            #message("*** error : can't analyze ", levels(rqall$Protein)[i], " for comparison.")
      
            processout <- rbind(processout, 
                                c(paste0("error : can't analyze ", levels(rqall$Protein)[i], " for comparison.")))
            write.table(processout, file=finalfile, row.names=FALSE)
        	
            tempresult <- list(result=NULL, valueresid=NULL, valuefitted=NULL, fittedmodel=NULL)
          
            for(k in 1:nrow(contrast.matrix)) {	
                tempresult$result <- rbind(tempresult$result,
                                           data.frame("Protein"=levels(rqall$Protein)[i], 
                                                      "Label"=row.names(contrast.matrix)[k],
                                                      "logFC"=NA,
                                                      "SE"=NA,
                                                      "Tvalue"=NA,
                                                      "DF"=NA,
                                                      "pvalue"=NA,
                                                      "issue"=NA))
            }
        
        } else {
    	    tempresult <- temp
            processout <- temp$processout
        }
    	
        temptempresult <- tempresult$result
    	    	
        ## need to add information about % missingness and %imputation
        if ( is.element("NumImputedFeature", colnames(data$RunlevelData)) ) {
            subprocess <- processall[processall$PROTEIN == as.character(unique(sub$PROTEIN)), ]
        
        	temptempresult <- .count.missing.percentage(contrast.matrix, temptempresult, sub, subprocess)
        }
    	
        ## comparison result table
        temptempresult$Label <- as.character(temptempresult$Label)
        out <- rbindlist(list(out, temptempresult))
        
        ## for checking model assumptions
        ## add residual and fitted after fitting the model
        if (data$SummaryMethod != "logOfSum" & class(temp) == "try-error") {
            if (nrow(sub) != 0) {
              sub$residuals <- NA
              sub$fitted <- NA
            }
        } else if (data$SummaryMethod != "logOfSum" & any(is.na(tempresult$result$pvalue))) { 
        	## even though there is no error for .fit.model function, still output can have empty fitted and residuals.
            sub$residuals <- NA
        	sub$fitted <- NA
              
        } else if (data$SummaryMethod != "logOfSum"){
            sub$residuals <- tempresult$valueresid
            sub$fitted <- tempresult$valuefitted
        }
    
        dataafterfit <- rbindlist(list(dataafterfit,sub), fill=TRUE)
    
        ## save fitted model
        outfitted <- c(outfitted, list(tempresult$fittedmodel))
    	
        ## progress
        setTxtProgressBar(pb, i)
    	
    } ### end protein loop
      
    close(pb)
      
    ## Combine the processout_MC with prior processout
    processout <- rbind(processout, 
                        "Comparisons for all proteins are done.- okay")
    write.table(processout, file=finalfile, row.names=FALSE)
      
    ## finalize result
    ## need to FDR per comparison
    out.all <- NULL
    
    out$Label <- as.character(out$Label)
    
    for(i in 1:length(unique(out$Label))) {
        outsub <- out[out$Label == unique(out$Label)[i], ]
        outsub <- data.frame(outsub, adj.pvalue=p.adjust(outsub$pvalue, method="BH"))
        out.all <- rbindlist(list(out.all, outsub))
        #out.all <- rbind(out.all, outsub)
    }
    out.all$Label <- factor(out.all$Label)
    
    ##
    processout <- rbind(processout, 
                        "Adjust p-values per comparison are calculated - okay.")
    write.table(processout, file=finalfile, row.names=FALSE)
    
    temp <- data$ProcessedData[!is.na(data$ProcessedData[,"ABUNDANCE"]) & 
                                   !is.na(data$ProcessedData[,"INTENSITY"]), ]
    temp <- temp[temp$ABUNDANCE > 2, ]
      
    if (abs(log2(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"]) < 
        abs(log10(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"])) {
        colnames(out.all)[3] <- "log2FC"
    }
    
    if (abs(log2(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"]) >
        abs(log10(temp[1, "INTENSITY"]) - temp[1, "ABUNDANCE"])) {
        colnames(out.all)[3] <- "log10FC"
    }
    
    ## change the format as data.frame
    out.all <- data.frame(out.all)
    
    ## change order of columns, 
    if (any(is.element(colnames(out.all), "ImputationPercentage"))) {
        out.all <- out.all[, c(1:7, 11, 8, 9, 10)]		
    } else if (any(is.element(colnames(out.all), "MissingPercentage"))) {
        out.all <- out.all[, c(1:7, 10, 8, 9)]
    } else {
    ## logOfSum output
        out.all <- out.all[, c(1:7, 9)]
    }
    
    ## if just one condition is completely missing, replace adjust pvalue as zero
    if (any(is.element(colnames(out.all), "issue"))){
        out.all[!is.na(out.all$issue) & out.all$issue == "oneConditionMissing", "adj.pvalue"] <- 0
    }
      
    ##
    processout <- rbind(processout, 
                        "Group comparison is done. - okay")
    write.table(processout, file=finalfile, row.names=FALSE)
    
    finalout <- list("ComparisonResult"=out.all, 
                     "ModelQC"=dataafterfit, 
                     "fittedmodel"=outfitted)
    return(finalout)	
}
    
    
#############################################
## fit.model
#############################################

#############################################
## check if measurements are missing for entire group
## if yes, length of group and length of contrast won't agree
#############################################
.chechGroupComparisonAgreement <- function(sub1, contrast.matrix) {
    tempSub <- as.numeric(as.character(unique(sub1[, c("GROUP")])))
    positionMiss <- setdiff(seq(1, length(contrast.matrix)), tempSub)
    contrast.matrix.sub1 <- contrast.matrix[tempSub]
    # either one of the groups for the comparison of interest is not present 

    return(list("sign"=length(setdiff(contrast.matrix[tempSub], 0)) < 2, 
                "positionMiss"=positionMiss))
}
    
#############################################
## check repeated (case-control? or time-course?)
#############################################
.checkRepeated <- function(data) {

    data.light <- data[data$LABEL=="L", ]	
    subjectByGroup <- table(data.light$SUBJECT_ORIGINAL, data.light$GROUP_ORIGINAL)
    subjectAppearances <- apply(subjectByGroup, 1, function(x) sum(x > 0))
    crossedIndicator <- any(subjectAppearances > 1)
    
    return(crossedIndicator)
}
    
#############################################
## check single subject for both case-control and time-course?
#############################################
.checkSingleSubject <- function(data) {

    temp <- unique(data[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
    temp$GROUP_ORIGINAL <- factor(temp$GROUP_ORIGINAL)
    temp1 <- xtabs(~ GROUP_ORIGINAL, data=temp)
    singleSubject <- all(temp1 == "1")
    
    return(singleSubject)	
}
    
#############################################
## check .checkSingleFeature
#############################################
.checkSingleFeature <- function(data) {

    sigleFeature <- length(unique(data$FEATURE)) < 2
    
    return(sigleFeature)	
}
    
#############################################
## check .checkTechReplicate
#############################################

.checkTechReplicate <-  function(data) {

    temp <- unique(data[, c("SUBJECT_NESTED", "RUN")])
    temp$SUBJECT_NESTED <- factor(temp$SUBJECT_NESTED)
    temp1 <- xtabs(~ SUBJECT_NESTED, data=temp)
    TechReplicate <- all(temp1 != "1")
    
    return(TechReplicate)	
}
    
    
#############################################
## check .checkRunByFeature
#############################################
# it might not be right
.checkRunbyFeature <- function(data) {

    data.light <- data[data$LABEL=="L", ]	
    RunByFeature <- table(data.light$RUN, data.light$FEATURE)
    emptyRow <- apply(RunByFeature, 1, sum)
    noRunFeature <- any(emptyRow == 0)
    
    return(noRunFeature)	
}

    
#############################################
## check .checkMissGroupByFeature
#############################################
.checkMissGroupByFeature <- function(data) {

    temp <- unique(data[, c("GROUP", "FEATURE")])
    temp1 <- xtabs(~ GROUP, data=temp)
    
    return(any(temp1 != temp1[1]))
}

    
#############################################
## check .checkMissRunByFeature
#############################################
.checkMissRunByFeature <- function(data) {

    temp <- unique(data[data$LABEL == "L", c("RUN", "FEATURE")])
    temp1 <- xtabs(~ RUN, data=temp)
    
    return(any(temp1 != length(unique(data$FEATURE))))
}


#############################################
## check .checkMissFeature for label-free missingness
#############################################
.checkMissFeature <- function(data) {

    dataByPeptide <- tapply(as.numeric(data$ABUNDANCE), 
                          list(data$FEATURE, data$GROUP_ORIGINAL), 
                          function(x) sum(x > 0, na.rm=TRUE))
    missPeptideInd <- apply(dataByPeptide, 1, function(x) any(x == 0 | is.na(x)))
    missingPeptides <- names(missPeptideInd)[missPeptideInd == TRUE]
    
    return(missingPeptides)
}
    
    
#############################################
## check .checkUnequalSubject
#############################################

.checkUnequalSubject <- function(data) {

    temp <- unique(data[data$LABEL == "L", c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
    temp1 <- xtabs(~ GROUP_ORIGINAL, data=temp)
    
    return(any(temp1 != temp1[1]))
}
    
    
########################################################
.fit.model.single <- function(contrast.matrix,
                              data,
                              TechReplicate,
                              singleSubject,
                              repeated,
                              origGroup,
                              processout) {
    
    ## input : output of run quantification
    data2 <- data
    data2$GROUP <- factor(data2$GROUP)
    data2$SUBJECT <- factor(data2$SUBJECT)
    
    ## if there is only one condition between two conditions, make error message and next
    if (length(unique(data2$GROUP_ORIGINAL)) == 1) {
    	
        ## each comparison
        allout <- NULL
          
        for (k in 1:nrow(contrast.matrix)) {
            
            ## choose each comparison
            contrast.matrix.sub <- matrix(contrast.matrix[k, ], nrow=1)
            row.names(contrast.matrix.sub) <- row.names(contrast.matrix)[k]
                
            if (any(levels(origGroup)[contrast.matrix.sub != 0] == unique(data2$GROUP_ORIGINAL))) {
                    
                processout <- rbind(processout, 
                                    paste0("*** error : results of Protein ", 
                                           unique(data2$PROTEIN), " for comparison ", 
                                           row.names(contrast.matrix.sub), " are NA because there are measurements only in Group ", 
                                           unique(data2$GROUP_ORIGINAL), "."))
                      
                ## need to check Inf vs -Inf
                if (contrast.matrix.sub[contrast.matrix.sub != 0 & 
                                        (levels(origGroup) == unique(data2$GROUP_ORIGINAL)) ] > 0) {
                          
                    out <- data.frame("Protein"=unique(data2$PROTEIN),
                                      "Label"=row.names(contrast.matrix.sub),
                                      "logFC"=Inf,
                                      "SE"=NA,
                                      "Tvalue"=NA,
                                      "DF"=NA,
                                      "pvalue"=NA,
                                      "issue"='oneConditionMissing')
                    
                } else {
                    out <- data.frame("Protein"=unique(data2$PROTEIN),
                                      "Label"=row.names(contrast.matrix.sub),
                                      "logFC"=(-Inf),
                                      "SE"=NA,
                                      "Tvalue"=NA,
                                      "DF"=NA,
                                      "pvalue"=NA,
                                      "issue"='oneConditionMissing')
                }         	
                	
            } else {
                processout <- rbind(processout, 
                                    paste0("*** error : results of Protein ", unique(data2$PROTEIN), 
                                           " for comparison ", row.names(contrast.matrix.sub), 
                                           " are NA because there are no measurements in both conditions."))
                        
                out <- data.frame("Protein"=unique(data2$PROTEIN),
                                  "Label"=row.names(contrast.matrix.sub),
                                  "logFC"=NA,
                                  "SE"=NA,
                                  "Tvalue"=NA,
                                  "DF"=NA,
                                  "pvalue"=NA,
                                  "issue"='completeMissing')
                	
            }
             	       	        
            allout <- rbind(allout, out)
            	
        } ## end loop for comparion    
      	
        finalresid <- NULL
        finalfitted <- NULL
        fit.full <- NULL
        	
    } else {
    	
        ## when subject is fixed, it is ok using lm function.
        ## when single feature, consider technical replicates for time-course.
        
        ## case-control
        if (!repeated) {
        	if (!TechReplicate | singleSubject) {
                fit.full <- lm(ABUNDANCE ~ GROUP, data=data2)
            } else {
                fit.full <- suppressMessages(try(lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data=data2),
                                                 TRUE))
                df.full <- suppressMessages(try(lm(ABUNDANCE ~ GROUP + SUBJECT, data=data2)$df.residual,
                                                TRUE))
            }
        } else { ## time-course
            if (singleSubject) {
                fit.full <- lm(ABUNDANCE ~ GROUP, 
                               data=data2)
            } else { ## no single subject
                if (!TechReplicate) {
                      fit.full <- suppressMessages(try(lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data=data2), 
                                       TRUE))
                      df.full <- suppressMessages(try(lm(ABUNDANCE ~ GROUP + SUBJECT, data=data2)$df.residual, 
                                    TRUE))
                } else {
                      fit.full <- suppressMessages(try(lmer(ABUNDANCE ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT), 
                                       data=data2), 
                                       TRUE))## SUBJECT==SUBJECT_NESTED here
                      df.full <- suppressMessages(try(lm(ABUNDANCE ~ GROUP + SUBJECT + GROUP:SUBJECT, 
                                    data=data2)$df.residual,
                                    TRUE))
                }
            }	
        } ## time-course
      
        ## get parameter from model
        if (class(fit.full) == "lm") {
        	Para <- .getParameterFixed(fit.full)
        } else {
        	Para <- .getParameterRandom(fit.full, df.full)
        }
      
        ## each comparison
        allout <- NULL
        	
        ## get condition IDs which are completely missing.
        emptycondition <- setdiff(levels(origGroup), unique(data2$GROUP_ORIGINAL))
          
        for (k in 1:nrow(contrast.matrix)) {
            
            ## choose each comparison
            contrast.matrix.sub <- matrix(contrast.matrix[k, ], nrow=1)
            row.names(contrast.matrix.sub) <- row.names(contrast.matrix)[k]
  	
            if (length(emptycondition) != 0) { # if there are any completely missing in any condition,
            		
                ## one by one comparison is simple. However, for linear combination of condition can be complicated
                ## get + and - condition separately
                count.pos <- levels(origGroup)[contrast.matrix.sub != 0 & contrast.matrix.sub > 0]
                count.neg <- levels(origGroup)[contrast.matrix.sub != 0 & contrast.matrix.sub < 0]
                		
                ## then check whether any + or - part is completely missing
                count.diff.pos <- intersect(levels(origGroup)[contrast.matrix.sub != 0 & 
                                                                contrast.matrix.sub > 0], 
                                            emptycondition)
                count.diff.neg <- intersect(levels(origGroup)[contrast.matrix.sub != 0 & 
                                                                contrast.matrix.sub < 0], 
                                            emptycondition)
                		
                ## positive side
                if (length(count.diff.pos) != 0) {
                    flag.issue.pos <- TRUE ## TRUE : there are problematic conditions
                } else {
                    flag.issue.pos <- FALSE ## FALSE : no issue about completely missing
                }
                		
                ## negative side
                if (length(count.diff.neg) != 0) {
                  flag.issue.neg <- TRUE ## TRUE : there are problematic conditions
                } else {
                  flag.issue.neg <- FALSE ## FALSE : no issue about completely missing
                }
          
                ## message and output  			
                if (flag.issue.pos & flag.issue.neg) { ## both sides are completely missing
                			
                    processout <- rbind(processout, 
                                        paste0("*** error : results of Protein ", unique(data2$PROTEIN), 
                                               " for comparison ", row.names(contrast.matrix.sub), 
                                               " are NA because there are no measurements in both conditions."))
		    
                    out <- data.frame("Protein"=unique(data2$PROTEIN), 
                                      "Label"=row.names(contrast.matrix.sub),
                                      "logFC"=NA, 
                                      "SE"=NA, 
                                      "Tvalue"=NA, 
                                      "DF"=NA, 
                                      "pvalue"=NA,
                                      "issue"='completeMissing')	
                  		
                } else if (flag.issue.pos | flag.issue.neg) {
                			
                    if (flag.issue.pos) {
                        issue.side <- count.diff.pos
                				
                		processout <- rbind(processout, 
                		                    paste0("*** error : results of Protein ", unique(data2$PROTEIN), 
                		                           " for comparison ", row.names(contrast.matrix.sub), 
                		                           " are NA because there are measurements only in Group ", 
                		                           paste(issue.side, collapse = ", "), "."))
                    	          	
                        out <- data.frame("Protein"=unique(data2$PROTEIN),
                                          "Label"=row.names(contrast.matrix.sub),
                                          "logFC"=(-Inf), 
                                          "SE"=NA, 
                                          "Tvalue"=NA, 
                                          "DF"=NA, 
                                          "pvalue"=NA,
                                          "issue"='oneConditionMissing')
                	}
                
                    if (flag.issue.neg) {
                	    issue.side <- count.diff.neg
                				
                		processout <- rbind(processout, 
                		                    paste0("*** error : results of Protein ", unique(data2$PROTEIN), 
                		                           " for comparison ", row.names(contrast.matrix.sub),
                		                           " are NA because there are measurements only in Group ", 
                		                           paste(issue.side, collapse = ", "), "."))
                    	          	
                        out <- data.frame("Protein"=unique(data2$PROTEIN),
                                          "Label"=row.names(contrast.matrix.sub),
                                          "logFC"=Inf, 
                                          "SE"=NA, 
                                          "Tvalue"=NA, 
                                          "DF"=NA, 
                                          "pvalue"=NA,
                                          "issue"='oneConditionMissing')
                	}       			
                			
                } else { ## then same as regulat calculation
                    contrast <- .make.contrast.free.single(fit.full, contrast.matrix.sub, data2)
                    out <- .estimableFixedRandom(Para, contrast)
                  
                    ## any error for out, just NA
                    if (is.null(out)) {
                        out <- data.frame("Protein"=unique(data2$PROTEIN), 
                                          "Label"=row.names(contrast.matrix.sub), 
                                          "logFC"=NA, 
                                          "SE"=NA, 
                                          "Tvalue"=NA, 
                                          "DF"=NA, 
                                          "pvalue"=NA, 
                                          "issue"=NA)
                    } else {
                        out <- data.frame("Protein"=unique(data2$PROTEIN), 
                                          "Label"=row.names(contrast.matrix.sub), 
                                          out, 
                                          "issue"=NA)	
                    }
                }
            } else {
                contrast <- .make.contrast.free.single(fit.full, contrast.matrix.sub, data2)
                out <- .estimableFixedRandom(Para, contrast)
                  
                ## any error for out, just NA
                if (is.null(out)) {
                    out <- data.frame("Protein"=unique(data2$PROTEIN), 
                                      "Label"=row.names(contrast.matrix.sub), 
                                      "logFC"=NA, 
                                      "SE"=NA, 
                                      "Tvalue"=NA, 
                                      "DF"=NA, 
                                      "pvalue"=NA, 
                                      "issue"=NA)
                } else {
                    out <- data.frame("Protein"=unique(data2$PROTEIN), 
                                      "Label"=row.names(contrast.matrix.sub), 
                                      out, 
                                      "issue"=NA)	
                }
            }
            
            allout <- rbind(allout, out)
        } ## end loop for comparion
          
        if (class(fit.full) == "lm") {  ## lm model
            finalresid <- fit.full$residuals
            finalfitted <- fit.full$fitted.values
        } else {   ## lmer model
            finalresid <- resid(fit.full)
            finalfitted <- fitted(fit.full)
        }
    } ## more than 2 conditions in the dataset
    
    finalout <- list("result"=allout, 
                     "valueresid"=finalresid, 
                     "valuefitted"=finalfitted, 
                     "fittedmodel"=fit.full, 
                     "processout"=processout)	
    return(finalout)
    
} ## .fit.model.single
    
    
#############################################
.ttest.logsum <- function(contrast.matrix, 
                          data, 
                          origGroup) {

    ## each comparison
    allout <- NULL
      
    for (k in 1:nrow(contrast.matrix)) {
        
        ## choose each comparison
        contrast.matrix.sub <- matrix(contrast.matrix[k,], nrow=1)
        row.names(contrast.matrix.sub) <- row.names(contrast.matrix)[k]
            
        ## get two groups from contrast.matrix
        datasub <- data[which(data$GROUP_ORIGINAL %in% origGroup[contrast.matrix.sub != 0]), ]
              
        ## t test
        sumresult <- try(t.test(datasub$ABUNDANCE~datasub$GROUP_ORIGINAL, 
                                var.equal=TRUE), 
                         silent=TRUE)
        
        if (class(sumresult) == "try-error") {
            out <- data.frame("Protein"=unique(data$PROTEIN),
                              "Label"=row.names(contrast.matrix.sub), 
                              "logFC"=NA, 
                              "SE"=NA, 
                              "Tvalue"=NA, 
                              "DF"=NA, 
                              "pvalue"=NA, 
                              "issue"=NA)
        } else {
            out <- data.frame("Protein"=unique(data$PROTEIN), 
                              "Label"=paste0(names(sumresult$estimate)[1], " - ", names(sumresult$estimate)[2]), 
                              "logFC"=sumresult$estimate[1] - sumresult$estimate[2],
                              "SE"=(sumresult$estimate[1] - sumresult$estimate[2]) / sumresult$statistic, 
                              "Tvalue"=sumresult$statistic,
                              "DF"=sumresult$parameter, 
                              "pvalue"=sumresult$p.value,
                              "issue"=NA)
            rownames(out) <- NULL
        }
        
        allout <- rbind(allout, out)
        
    } # end loop for comparion
      
    finalout <- list("result"=allout, 
                     "valueresid"=NULL, 
                     "valuefitted"=NULL, 
                     "fittedmodel"=NULL)	
    return(finalout)
}
    
    
#############################################
.count.missing.percentage <- function(contrast.matrix, temptempresult, sub, subprocess) {

    totaln.cell <- aggregate(PROTEIN ~ GROUP_ORIGINAL, 
                             data=subprocess[subprocess$LABEL == "L", ], 
                             length) 
    ## just for count total measurement, use PROTEIN instead of ABUNDANCE in order to prevent to remove NA in ABUNDANCE
    colnames(totaln.cell)[colnames(totaln.cell) == "PROTEIN"] <- "totalN"
    
    totaln.measured <- aggregate(NumMeasuredFeature ~ GROUP_ORIGINAL, data=sub, sum, na.rm=TRUE)
    totaln.imputed <- aggregate(NumImputedFeature ~ GROUP_ORIGINAL, data=sub, sum, na.rm=TRUE)
    
    totaln <- merge(totaln.cell, totaln.measured, by="GROUP_ORIGINAL", all=TRUE)
    totaln <- merge(totaln, totaln.imputed, by="GROUP_ORIGINAL", all=TRUE)
    
    if (any(is.element(colnames(sub), "NumImputedFeature"))) {
        temptempresult$MissingPercentage <- NA
        temptempresult$ImputationPercentage <- NA
    } else {
        temptempresult$MissingPercentage <- NA
    }
    
    for (k in 1:nrow(contrast.matrix)) {
        
        ## choose each comparison
        contrast.matrix.sub <- matrix(contrast.matrix[k, ], nrow=1)
        row.names(contrast.matrix.sub) <- row.names(contrast.matrix)[k]
            
        condition.needed <- contrast.matrix.sub != 0
            
        MissingPercentage.new <- NA
        ImputationPercentage.new <- NA
            
        ## total # missing
        MissingPercentage.new <- 1 - 
            sum(totaln[condition.needed, "NumMeasuredFeature"], na.rm=TRUE) / 
            sum(totaln[condition.needed, "totalN"], na.rm=TRUE)
    	
        ## imputed intensity
        if (any(is.element(colnames(sub), "NumImputedFeature"))) {
            ImputationPercentage.new <- sum(
                totaln[condition.needed, "NumImputedFeature"], na.rm=TRUE) / 
                sum(totaln[condition.needed, "totalN"], na.rm=TRUE)
            
            temptempresult[temptempresult$Label == row.names(contrast.matrix.sub), 
                           "MissingPercentage"] <- MissingPercentage.new
            temptempresult[temptempresult$Label == row.names(contrast.matrix.sub), 
                           "ImputationPercentage"] <- ImputationPercentage.new
        } else {
            temptempresult[temptempresult$Label == row.names(contrast.matrix.sub), 
                           "MissingPercentage"] <- MissingPercentage.new
        }
    } # end loop for multiple comparisons
    return(temptempresult)
}
	
