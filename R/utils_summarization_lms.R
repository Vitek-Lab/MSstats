
##################################
########### make contrast ########
##################################

##================================
## .make.contrast.free: 
## label-free: equal or unequal subjects per group
## fixed or random subject
##================================

.make.contrast.free <- function(fit, contrast.matrix, sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
    } else {
        coef_name <- names(fixef(fit))
    }
    
    ## intercept in label-free
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(0, length(temp))
    names(intercept_c) <- temp
    if (length(temp) == 0) {
        intercept_c <- NULL
    }
    
    ## feature in label-free
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    feature_c <- rep(0, length(temp))
    names(feature_c) <- temp
    if (length(temp) == 0) {
        feature_c <- NULL
    }
    
    ## subject_nested in label-free
    temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
    if (length(temp) > 0) {
        temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\.")), nrow=2))
        temp2 <- as.vector(xtabs(~ temp1[, 1]))	
        tempdata <- fit$model
        group_levels <- unique(tempdata$GROUP)
        sub_contrast <- contrast.matrix[as.numeric(as.character(group_levels))]
        ## the base is alway be the first SUBJECT_NESTED
        ## (SUBJECT_NESTED1.1)
        temp3 <- temp2
        if (length(temp2) == length(sub_contrast)) {
            ## this statement goes first, otherwise length of temp3 becomes > 1
            temp3[1] <- temp2[1] + 1
        } else {
            temp3 <- c(1, temp3)
        }
        # subjectNested_c<-rep(contrast.matrix/(temp3),temp2) ## in case of unequal sample per group, wrong
        subjectNested_c <- rep(sub_contrast / (temp3), temp3)[-1]
        names(subjectNested_c) <- temp
    } else {
        ## length(temp) == 0
        subjectNested_c <- NULL
    }
    
    ## subject in label-free: for time-course
    # temp<-coef_name[grep("SUBJECT",coef_name)[!grep("SUBJECT",coef_name)%in%c(grep(":",coef_name),grep("NESTED",coef_name))]]
    # if(length(temp)>0){
    # 	subject_c<-rep(0,length(temp))
    # 	names(subject_c)<-temp
    # }
    # if(length(temp)==0) subject_c<-NULL
    
    temp <- coef_name[setdiff(grep("SUBJECT", coef_name), 
                              grep(":|NESTED", coef_name))]
	if (length(temp) > 0) {
	    # subject_c<-rep(0,length(temp))
	    # names(subject_c)<-temp
	    tempdata <- fit$model
	    group_levels <- unique(tempdata$GROUP)
	    labels <- paste("GROUP", group_levels, sep="")
	    patients <- NULL
	    for (i in 1:length(group_levels)) {
	        sub <- tempdata[tempdata$GROUP == group_levels[i], ]
	        sub_patients <- cbind(
	            GROUP=paste("GROUP", group_levels[i], sep=""), 
	            SUBJECT=paste("SUBJECT", as.character(unique(sub$SUBJECT)), sep=""), 
	            Value=contrast.matrix[as.numeric(as.character(group_levels[i]))])
	        patients <- data.frame(rbind(patients, sub_patients))
	    }
	    patient_count <- tapply(patients$SUBJECT, patients$GROUP, 
	                            function(x) length(unique(x)))
	    patient_seq <- rep(0, length(temp))
	    for (i in 1:length(as.character(patients$SUBJECT))) {
	        match <- any(temp == as.character(patients$SUBJECT)[i])
	        if (match & as.numeric(as.character(patients$Value[i])) != 0) {
	            res <- temp == as.character(patients$SUBJECT)[i]
	            index <- which(res == TRUE)
	            group <- as.character(patients[i, ]$GROUP)
	            count <- as.numeric(patient_count[names(patient_count) == group])
	            value <- as.numeric(as.character(patients[i, ]$Value))
	            patient_value <- c(rep(0, index - 1), 
	                               value / count, 
	                               rep(0, length(temp) - index))
	        } else {
	            patient_value <- rep(0, length(temp))
	        }
	        patient_seq <- patient_value + patient_seq
	    }
	    subject_c <- patient_seq
	    names(subject_c) <- temp
	} else {
	    ## length(temp) == 0
	    subject_c <- NULL
	}
    
    ## subject by group in label-free: only for time-course - SUBJECT and GROUP (order) even GROUP:SUBJECT in model
    # temp<-coef_name[intersect(grep("SUBJECT",coef_name),grep("GROUP",coef_name))]
    
    # tempSub<-unique(sub1[,c("GROUP","SUBJECT")])
    # tempSub1<-xtabs(~GROUP,data=tempSub)
    # tempSub2<-tempSub1[-1]
    # if(length(temp)>0){
    # 	temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
    # 	temp2<-as.vector(xtabs(~temp1[,2]))	## count per GROUP
    # 	#gs_c<-rep(as.vector(contrast.matrix[-1]/(tempSub2)),temp2[1]) ## assume no missing for group and subject
    # 	# when Group completely  missing
    # 	sub.matrix<-contrast.matrix[unique(tempSub$GROUP)]
    # 	gs_c<-rep(as.vector(sub.matrix[-1]/(tempSub2)),each=temp2[1])
    # 	names(gs_c)<-temp
    # }
    # if(length(temp)==0) gs_c<-NULL
    
    temp <- coef_name[intersect(grep("SUBJECT", coef_name), 
                                grep("GROUP", coef_name))]
    if (length(temp) > 0) {
        # subject_c<-rep(0,length(temp))
        # names(subject_c)<-temp
        tempdata <- fit$model
        group_levels <- unique(tempdata$GROUP)
        labels <- paste("GROUP", group_levels, sep="")
        patients <- NULL
        for (i in 1:length(group_levels)) {
            sub <- tempdata[tempdata$GROUP == group_levels[i], ]
            sub_patients <- cbind(
                GROUP=paste("GROUP", group_levels[i], sep=""), 
                SUBJECT=paste("SUBJECT", as.character(unique(sub$SUBJECT)), sep=""), 
                Value=contrast.matrix[as.numeric(as.character(group_levels[i]))])
            patients <- data.frame(rbind(patients, sub_patients))
        }
        patient_count <- tapply(patients$SUBJECT, patients$GROUP, 
                                function(x) length(unique(x)))
        interaction_seq <- rep(0, length(temp))
        interaction_labels <- paste(as.character(patients$GROUP), 
                                    as.character(patients$SUBJECT), sep=":")
        for (i in 1:length(as.character(patients$SUBJECT))) {
            match <- any(temp == interaction_labels[i])
            if (match & as.numeric(as.character(patients$Value[i])) != 0) {
                res <- temp == interaction_labels[i]
                index <- which(res == TRUE)
                group <- as.character(patients[i, ]$GROUP)
                count <- as.numeric(patient_count[names(patient_count) == group])
                value <- as.numeric(as.character(patients[i, ]$Value))
                interaction_value <- c(rep(0, index - 1), 
                                       value/count, 
                                       rep(0, length(temp) - index))
            } else {
                interaction_value <- rep(0, length(temp))
            }
            interaction_seq <- interaction_value + interaction_seq
        }
        gs_c <- interaction_seq
        names(gs_c) <- temp
    } else {
	    ## length(temp)==0
        gs_c <- NULL
	}
	
    ## group in label-free (different from labeled)
    temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
    ## when there are some groups which are all missing
    tempSub <- as.numeric(as.character(unique(sub1[, c("GROUP")])))
    tempcontrast <- contrast.matrix[tempSub]
    group_c <- tempcontrast[-1]  # for label-free, remove the first element
    names(group_c) <- temp
	if (length(temp) == 0) {
	    group_c <- NULL
	}
    
	## feature by group in label-free (different from labeled)
	temp <- coef_name[intersect(grep("GROUP", coef_name), 
	                            grep("FEATURE", coef_name))]
	tempSub <- unique(sub1[, c("GROUP", "FEATURE")])
	tempSub1 <- xtabs(~ GROUP, data=tempSub)
	tempSub2 <- tempSub1[-1]
	if (length(temp) > 0) {
	    temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\:")), nrow=2))
	    temp2 <- as.vector(xtabs(~ temp1[, 2]))
	    gf_c <- rep(contrast.matrix[-1] / (tempSub2), temp2)
	    names(gf_c) <- temp
	} else {
	    gf_c <- NULL
	}
	
	## be careful of the order
	contrast <- c(intercept_c, feature_c, subjectNested_c, subject_c, group_c, 
	              gs_c, gf_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}


##================================
## .make.contrast.based: 
## label-based: equal or unequal subjects per group
## fixed or random subject
##================================

.make.contrast.based <- function(fit, contrast.matrix, sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
	} else {
	    coef_name <- names(fixef(fit))
	}
    
    ## intercept in label-based
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(0, length(temp))
    names(intercept_c) <- temp
    if (length(temp) == 0) {
        intercept_c <- NULL
    }
    
	## feature in label-based
	temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
	feature_c <- rep(0, length(temp))
	names(feature_c) <- temp
	if (length(temp) == 0) {
	    feature_c <- NULL
	}
	
	## subject_nested in label-based
	temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
	if (length(temp) > 0) {
	    temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\.")), nrow=2))
	    temp2 <- as.vector(xtabs(~ temp1[, 1]))	
		## the base is alway be the first SUBJECT_NESTED
		## (SUBJECT_NESTED0.0)
		temp3 <- temp2
		## in label-free: temp3 <- temp2; temp3[1] <- temp2[1]+1
		subjectNested_c <- rep(contrast.matrix / (temp3), temp2)
		names(subjectNested_c) <- temp
	} else {
	    ## length(temp) == 0
	    subjectNested_c <- NULL
	}
	
	## subject in label-based
	temp <- coef_name[setdiff(grep("SUBJECT", coef_name), 
	                          grep(":|NESTED", coef_name))]
	if (length(temp) > 0) {
	    subject_c <- rep(0, length(temp))
	    names(subject_c) <- temp
	} else {
	    subject_c <- NULL
	}
	
	## subject by group in label-based
	temp <- coef_name[intersect(grep("SUBJECT", coef_name), 
	                            grep("GROUP", coef_name))]
	tempSub <- unique(sub1[ ,c("GROUP", "SUBJECT")])
	tempSub1 <- xtabs(~ GROUP, data=tempSub)
	tempSub2 <- tempSub1[-1]
	if (length(temp) > 0) {
	    temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\:")), nrow=2))
	    temp2 <- as.vector(xtabs(~ temp1[, 2]))	
	    gs_c <- rep(as.numeric(contrast.matrix) / tempSub2, temp2)
	    names(gs_c) <- temp
	} else {
	    gs_c <- NULL
	}
	
	## subject_original_nested in label-based
	temp <- coef_name[grep("SUBJECT_ORIGINAL_NESTED", coef_name)]
	if (length(temp) > 0) {
	    temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\.")), nrow=2))
	    temp2 <- as.vector(xtabs(~ temp1[, 1]))
	    ## the base is alway be the first SUBJECT_ORIGINAL_NESTED
		## (SUBJECT_ORIGINAL_NESTED1.1)
		temp3 <- temp2
		temp3[1] <- temp2[1] + 1
		subjectOriginalNested_c <- rep(contrast.matrix / temp3, temp2)
		names(subjectOriginalNested_c) <- temp
	} else {
	    subjectOriginalNested_c <- NULL
	}
	
	## group in label-based (different from label-free)
	temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
	## when there are some groups which are all missing
	tempSub <- as.numeric(as.character(unique(sub1[, "GROUP"])))
	group_c <- contrast.matrix[tempSub]
	## in label-free: group_c<-contrast.matrix[-1]
	names(group_c) <- temp
	if (length(temp) == 0) {
	    group_c <- NULL
	}
	
	## run in label-based
	temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
	run_c <- rep(0, length(temp))
	names(run_c) <- temp
	if (length(temp) == 0) {
	    run_c <- NULL
	}
	
	## feature by group in label-based (different from label-free)
	temp <- coef_name[intersect(grep("GROUP", coef_name), 
	                            grep("FEATURE", coef_name))]
	tempSub <- unique(sub1[, c("GROUP", "FEATURE")])
	tempSub1 <- xtabs(~ GROUP, data=tempSub)
	tempSub2 <- tempSub1[-1]
	if (length(temp) > 0) {
	    temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\:")), nrow=2))
		temp2 <- as.vector(xtabs(~ temp1[, 2]))
		gf_c <- rep(as.numeric(contrast.matrix) / tempSub2, temp2)
		# in label-free: gf_c<-rep(contrast.matrix[-1]/(tempSub2),temp2)
		names(gf_c) <- temp
	} else {
	    gf_c <- NULL
	}
	
	## run by feature in label-based
	temp <- coef_name[intersect(grep("RUN", coef_name), 
	                            grep("FEATURE", coef_name))]
	if (length(temp) > 0) {
	    rf_c <- rep(0, length(temp))
	    names(rf_c) <- temp
	} else {
	    rf_c<-NULL
	}

	contrast <- c(intercept_c, feature_c, subjectNested_c, subject_c, 
	              subjectOriginalNested_c, group_c, run_c, gs_c, gf_c, rf_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}


##########################################################################################

.make.contrast.run.quantification <- function(fit, contrast.matrix, sub1, 
                                              labeled) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
	} else {
	    coef_name <- names(fixef(fit))
	}
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    if (length(temp) > 0) {
        intercept_c <- rep(1, length(temp))
        names(intercept_c) <- temp
    } else {
        intercept_c <- NULL
    }
    
    ## feature
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("FEATURE", "RUN")])
        tempSub1 <- xtabs(~ RUN + FEATURE, data=tempSub)
        tempSub2 <- tempSub1[contrast.matrix == 1, ]
        feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
        names(feature_c) <- temp
    } else {
        feature_c <- NULL
    }
    
    ## run: different with other quantification - first try
    temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
    if (!labeled) { 
        ## label-free
        if (length(temp) > 0) {
            run_c <- contrast.matrix[-1]
            names(run_c) <- temp
        } else {
			run_c <- NULL
		}
	} else { 
	    ## label-based
	    if (length(temp) > 0) {
	        run_c <- rep(1 / nlevels(sub1$RUN), length(temp))
	        names(run_c) <- temp
	    } else {
	        run_c <- NULL
	    }
	}
    
    ## ref
    temp <- coef_name[grep("ref", coef_name)]
    if (length(temp) > 0) {
        if (nlevels(sub1$LABEL) == 2) {
            ref_levels <- levels(sub1$ref)
			ref_c <- contrast.matrix[1:(length(ref_levels) - 1)]
        }
        names(ref_c) <- temp
	} else {
	    ref_c <- NULL
	}
    
    ## run by feature
    temp <- coef_name[intersect(grep("RUN", coef_name), 
                                grep("FEATURE", coef_name))]
    tempSub <- dim(unique(sub1[, c("RUN", "FEATURE")]))[1]
    if (length(temp) > 0) {
        rf_c <- rep(1 / tempSub, length(temp))
        names(rf_c) <- temp
    } else {
        rf_c <- NULL
    }
    
    ## subject_nested
    temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
	if (length(temp) > 0) {
	    if (nlevels(sub1$LABEL) == 2) {
	        subjectNested_c <- contrast.matrix
	    } else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    subjectNested_c <- contrast.matrix[-1]
		}
		names(subjectNested_c) <- temp
	} else {
	    subjectNested_c <- NULL
	}

	contrast <- c(intercept_c, feature_c, run_c, ref_c, rf_c, subjectNested_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##########################################################################################
.make.contrast.run.quantification.reference <- function(fit, contrast.matrix, 
                                                        sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
    } else {
        coef_name <- names(fixef(fit))
    }
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(1, length(temp))
    names(intercept_c) <- temp
    if (length(temp) == 0) {
        intercept_c <- NULL
    }
    
    ## feature
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("FEATURE", "RUN")])
        tempSub1 <- xtabs(~ RUN + FEATURE, data=tempSub)
        tempSub2 <- tempSub1[contrast.matrix == 1, ]
        feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
        names(feature_c) <- temp
    } else {
        feature_c <- NULL
    }
    
    ## run
    temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        run_c <- rep(1 / nlevels(sub1$RUN), length(temp))
        names(run_c) <- temp
    } else {
        run_c <- NULL
    }
    
    ## ref
    temp <- coef_name[grep("ref", coef_name)]
    if (length(temp) > 0) {
        ref_c <- rep(0, length(temp))
        names(ref_c) <- temp
    } else {
        ref_c <- NULL
    }

	contrast <- c(intercept_c, feature_c, run_c, ref_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##================================
## .make.contrast.group.quantification: 
## label-based/label-free; single/multiple features
## all fixed subject and run
##================================

.make.contrast.group.quantification <- function(fit, contrast.matrix, sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
    } else {
        coef_name <- names(fixef(fit))
    }
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(1, length(temp))
    names(intercept_c) <- temp
    if (length(temp) == 0) {
        intercept_c <- NULL
    }
    
    ## feature
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    
    # if(length(temp)>0){
    # 	if(class(fit)=="lm"){
    # 		feature_c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
    # 	}else{
    # 		### need to fix fit$omdel
    # 		tempfeature<-eval(getCall(fit)$data)
    # 		tempfeature$FEATURE<-factor(tempfeature$FEATURE)
    # 		feature_c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
    # 	}
    # 	names(feature_c)<-temp
    # }else{
    # 	feature_c<-NULL
    # }
    
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("FEATURE", "GROUP")])
        tempSub1 <- xtabs(~ GROUP + FEATURE, data=tempSub)
        if (nlevels(sub1$LABEL) == 2) {
            tempSub1 <- tempSub1[-1, ]
            tempSub2 <- tempSub1[contrast.matrix == 1, ]
        } else if (nlevels(sub1$LABEL) == 1) {
            ## label-free
            tempSub2 <- tempSub1[contrast.matrix == 1, ]
        }
        feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
        names(feature_c) <- temp
    } else {
        feature_c <- NULL
    }
    
    ## subject_nested
    temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
    if (length(temp) > 0) {
        temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\.")), nrow=2))
        temp2 <- as.vector(xtabs(~ temp1[, 1]))
        ## the base is alway be the first SUBJECT_NESTED
        ## (SUBJECT_NESTED0.0)
        if (nlevels(sub1$LABEL) == 2) {
            temp3 <- temp2
            subjectNested_c <- rep(contrast.matrix / temp3, temp2)
        } else if (nlevels(sub1$LABEL) == 1) {
            ## label-free
            temp3 <- temp2
            if (length(temp2) == length(contrast.matrix)) {
                ## this statement goes first, otherwise length of temp3 >1
                temp3[1] <- temp2[1] + 1
            } else {
                temp3 <- c(1, temp3)
            }
            subjectNested_c <- rep(contrast.matrix / temp3, temp3)[-1]
        }
        names(subjectNested_c) <- temp
    } else {
        subjectNested_c <- NULL
    }
    
    ## subject_original
    temp <- coef_name[setdiff(grep("SUBJECT_ORIGINAL", coef_name), 
                              grep(":", coef_name))]
	if (length(temp) > 0) {
	    
	    # if(class(fit)=="lm"){
	    # 	subjectOrig_c<-rep(1/nlevels(fit$model$SUBJECT_ORIGINAL),length(temp))
	    # }else{
	    # 	### need to fix fit$omdel
	    # 	subjectOrig_c<-rep(1/nlevels(eval(getCall(fit)$data)$SUBJECT_ORIGINAL),length(temp))
	    # }
	    
	    subjectOrig_c <- rep(1 / nlevels(sub1$SUBJECT_ORIGINAL), length(temp))
	    names(subjectOrig_c) <- temp
	} else {
	    subjectOrig_c <- NULL
	}
    
    ## subject_original: group
    temp <- coef_name[intersect(grep("SUBJECT_ORIGINAL", coef_name), 
                                grep("GROUP", coef_name))]
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("SUBJECT_ORIGINAL", "GROUP")])
        tempSub1 <- xtabs(~ GROUP, data=tempSub)
		tempSub2 <- tempSub1[-1]
		temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\:")), nrow=2))
		temp2 <- as.vector(xtabs(~ temp1[, 2]))	
		if (nlevels(sub1$LABEL) == 2) {
		    subjectOrigGroup_c <- 
		        rep(as.numeric(contrast.matrix) / tempSub2, temp2)
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    subjectOrigGroup_c <- rep(contrast.matrix[-1] / tempSub2, temp2)
		}
		names(subjectOrigGroup_c) <- temp
    } else {
        subjectOrigGroup_c <- NULL
    }
    
    ## group
    temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        if (nlevels(sub1$LABEL) == 2) {
            group_c <- contrast.matrix
        } else if (nlevels(sub1$LABEL) == 1) {
            ## label-free
            group_c <- contrast.matrix[-1]
        }
        names(group_c) <- temp
        
        ## if some group's coef is NA, need to remove
        #tempname<-rownames(summary(fit)$coefficients)
        #tempname1<-tempname[grep("GROUP",tempname)[!grep("GROUP",tempname)%in%grep(":",tempname)]]
        #group_c<-group_c[names(group_c)==tempname1]
        
    } else {
        group_c <- NULL
    }
    
    ## run
    temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        
        # if(class(fit)=="lm"){
        # 	run_c<-rep(1/nlevels(fit$model$RUN),length(temp))
        # }else{
        # 	### need to fix fit$omdel
        # 	run_c<-rep(1/nlevels(eval(getCall(fit)$data)$RUN),length(temp))
        # }
        
        run_c <- rep(1 / nlevels(sub1$RUN), length(temp))
        names(run_c) <- temp
    } else {
        run_c <- NULL
    }
    
    ## feature by group
    temp <- coef_name[intersect(grep("GROUP", coef_name), 
                                grep("FEATURE", coef_name))]
    if (length(temp) > 0) {
		tempSub <- unique(sub1[, c("GROUP", "FEATURE")])
		tempSub1 <- xtabs(~ GROUP, data=tempSub)
		tempSub2 <- tempSub1[-1]
		temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\:")), nrow=2))
		temp2 <- as.vector(xtabs(~ temp1[, 2]))
		if (nlevels(sub1$LABEL) == 2) {
		    gf_c <- rep(as.numeric(contrast.matrix) / tempSub2, temp2)
		} else if (nlevels(sub1$LABEL)==1) {
		    ## label-free
		    gf_c <- rep(contrast.matrix[-1] / tempSub2, temp2)
		}
		names(gf_c) <- temp
    } else {
        gf_c <- NULL
    }
    
    ## run by feature
    temp <- coef_name[intersect(grep("RUN", coef_name), 
                                grep("FEATURE", coef_name))]
    tempSub <- dim(unique(sub1[, c("RUN", "FEATURE")]))[1]
    if (length(temp) > 0) {
        rf_c <- rep(1 / tempSub, length(temp))
        names(rf_c) <- temp
    } else {
        rf_c <- NULL
    }
    
    contrast <- c(intercept_c, feature_c, subjectNested_c, subjectOrig_c, 
                  subjectOrigGroup_c, group_c, run_c, gf_c, rf_c)
    if (class(fit) == "lm") {
        contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
    
	return(contrast1)
}

##================================
## .make.contrast.group.quantification.reference: 
## label-based/label-free; single/multiple features
## all fixed subject and run
##================================

.make.contrast.group.quantification.reference <- function(fit, contrast.matrix, 
                                                          sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
    } else {
        coef_name <- names(fixef(fit))
    }
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(1, length(temp))
    names(intercept_c) <- temp
    if (length(temp) == 0) {
        intercept_c <- NULL
    }
    
    ## feature
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    
    # if(length(temp)>0){
    # 	if(class(fit)=="lm"){
    # 		feature_c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
    # 	}else{
    # 		### need to fix fit$omdel
    # 		tempfeature<-eval(getCall(fit)$data)
    # 		tempfeature$FEATURE<-factor(tempfeature$FEATURE)
    # 		feature_c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
    # 	}
    # 	names(feature_c)<-temp
    # }else{
    # 	feature_c<-NULL
    # }
    
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("FEATURE", "GROUP")])
        tempSub1 <- xtabs(~ GROUP + FEATURE, data=tempSub)
        if (nlevels(sub1$LABEL) == 2) {
            tempSub1 <- tempSub1[-1, ]
            tempSub2 <- tempSub1[contrast.matrix == 1, ]
        } else if (nlevels(sub1$LABEL) == 1) {
            ## label-free
            tempSub2 <- tempSub1[contrast.matrix == 1, ]
        }
        feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
        names(feature_c) <- temp
    } else {
        feature_c <- NULL
    }
    
    ## subject_nested
    temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
    if (length(temp) > 0) {
        subjectNested_c <- rep(0, length(temp))
        names(subjectNested_c) <- temp
    } else {
        subjectNested_c <- NULL
    }
    
    ## subject_original
    temp <- coef_name[grep("SUBJECT_ORIGINAL", coef_name)]
    if (length(temp) > 0) {
        subjectOrig_c <- rep(0, length(temp))
        names(subjectOrig_c) <- temp
    } else {
        subjectOrig_c <- NULL
    }
    
    ## subject_original:group
    temp <- coef_name[intersect(grep("SUBJECT_ORIGINAL", coef_name), 
                                grep("GROUP", coef_name))]
    if (length(temp) > 0) {
        subjectOrigGroup_c <- rep(0, length(temp))
        names(subjectOrigGroup_c) <- temp
    } else {
        subjectOrigGroup_c <- NULL
    }
    
    ## group
    temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        group_c <- rep(0, length(temp))
        names(group_c) <- temp
    } else {
        group_c <- NULL
    }

    ## run
    temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        
        # if(class(fit)=="lm"){
        # 	run_c<-rep(1/nlevels(fit$model$RUN),length(temp))
        # }else{
        # 	### need to fix fit$omdel
        # 	temprun<-eval(getCall(fit)$data)
        # 	tempfeature$FEATURE<-factor(tempfeature$FEATURE)
        # 	feature_c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
        # 	run_c<-rep(1/nlevels(eval(getCall(fit)$d)$RUN),length(temp))
        # }
        
        run_c <- rep(1 / nlevels(sub1$RUN), length(temp))
        names(run_c) <- temp
    } else {
        run_c <- NULL
    }
    
    ## feature by group
    temp <- coef_name[intersect(grep("GROUP", coef_name), 
                                grep("FEATURE", coef_name))]
    if (length(temp) > 0) {
        gf_c <- rep(0, length(temp))
        names(gf_c) <- temp
    } else {
        gf_c <- NULL
    }
    
    ## run by feature
    temp <- coef_name[intersect(grep("RUN", coef_name), 
                                grep("FEATURE", coef_name))]
    tempSub <- dim(unique(sub1[, c("RUN", "FEATURE")]))[1]
    if (length(temp) > 0) {
        rf_c <- rep(1 / tempSub, length(temp))
        names(rf_c) <- temp
    } else {
        rf_c <- NULL
    }

	contrast <- c(intercept_c, feature_c, subjectNested_c, subjectOrig_c, 
	              subjectOrigGroup_c, group_c, run_c, gf_c, rf_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}

	return(contrast1)
}


##================================
## .make.contrast.subject.quantification.single: 
## label-based/label-free; single features
## all fixed subject and run
##================================

.make.contrast.subject.quantification.single <- function(fit, contrast.matrix, 
                                                         sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
    } else {
        coef_name <- names(fixef(fit))
    }
    
    #contrastList<-unique(sub1[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
    #contrastGroup<-rep(0,nlevels(sub1$GROUP_ORIGINAL))
    #contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP_ORIGINAL"])]<-1
    
    ## for label-based	
	if (nlevels(sub1$LABEL) == 2) {
	    ## remove GROUP==0
	    contrastList <- unique(sub1[, c("GROUP", "SUBJECT_ORIGINAL")])[-1, ]
	    contrastList$GROUP <- factor(contrastList$GROUP)
		contrastGroup <- rep(0, nlevels(sub1$GROUP) - 1)
		contrastGroup[as.numeric(contrastList[contrast.matrix == 1, 
		                                      "GROUP"])] <- 1
	} else {
	    ## for label-free
	    contrastList <- unique(sub1[, c("GROUP", "SUBJECT_ORIGINAL")])
	    contrastGroup <- rep(0, nlevels(sub1$GROUP))
	    contrastGroup[as.numeric(contrastList[contrast.matrix == 1, 
	                                          "GROUP"])] <- 1
	}
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(1, length(temp))
    names(intercept_c) <- temp
    if (length(temp) == 0) {
        intercept_c <- NULL
    }
    
    ## feature
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    
    # if(length(temp)>0){
    # 	if(class(fit)=="lm"){
    # 		feature_c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
    # 	}else{
    # 		### need to fix fit$omdel
    # 		tempfeature<-eval(getCall(fit)$data)
    # 		tempfeature$FEATURE<-factor(tempfeature$FEATURE)
    # 		feature_c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
    # 	}
    # 	names(feature_c)<-temp
    # }else{
    # 	feature_c<-NULL
    # }
    
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("FEATURE", "SUBJECT_NESTED")])
        tempSub1 <- xtabs(~ SUBJECT_NESTED + FEATURE, data=tempSub)
        if (nlevels(sub1$LABEL) == 2) {
            tempSub1 <- tempSub1[-1, ]
			tempSub2 <- tempSub1[contrast.matrix == 1, ]
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
			tempSub2 <- tempSub1[contrast.matrix == 1, ]
		}
        feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
		names(feature_c) <- temp
	} else {
	    feature_c <- NULL
	}
    
    ## subject_nested
    temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
    if (length(temp) > 0) {
        if (nlevels(sub1$LABEL) == 2) {
            subjectNested_c <- contrast.matrix
        } else if (nlevels(sub1$LABEL) == 1) {
            ## label-free
            subjectNested_c <- contrast.matrix[-1]
		}
		names(subjectNested_c) <- temp
	} else {
	    subjectNested_c <- NULL
	}
    
    ## group
    temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        if (nlevels(sub1$LABEL) == 2) {
            group_c <- contrastGroup
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    group_c <- contrastGroup[-1]
		}
        names(group_c) <- temp
    } else {
        group_c <- NULL
    }
    
    ## run
    temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        
        #run_c<-rep(1/nlevels(fit$model$RUN),length(temp))
        
        ## when no technical replicate: subject_nested = run
		run_c <- contrast.matrix[-1]
		## however, with technical replicate: subject_nested != run, need others
		names(run_c) <- temp
    } else {
        run_c <- NULL
    }
    
    ## feature by group
    temp <- coef_name[intersect(grep("GROUP", coef_name), 
                                grep("FEATURE", coef_name))]
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("GROUP", "FEATURE")])
        tempSub1 <- xtabs(~ GROUP, data=tempSub)
		tempSub2 <- tempSub1[-1]
		temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\:")), nrow=2))
		temp2 <- as.vector(xtabs(~ temp1[, 2]))
		if (nlevels(sub1$LABEL) == 2) {
		    gf_c <- rep(as.numeric(contrastGroup) / tempSub2, temp2)
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    gf_c <- rep(contrastGroup[-1] / tempSub2, temp2)
		}
		names(gf_c) <- temp
    } else {
        gf_c <- NULL
    }
    
    ## run by feature
    temp <- coef_name[intersect(grep("RUN", coef_name), 
                                grep("FEATURE", coef_name))]
    tempSub <- dim(unique(sub1[, c("RUN", "FEATURE")]))[1]
    if (length(temp) > 0) {
        rf_c <- rep(1 / tempSub, length(temp))
        names(rf_c) <- temp
    } else {
        rf_c <- NULL
    }
    
    contrast <- c(intercept_c, feature_c, subjectNested_c, group_c, run_c, 
                  gf_c, rf_c)
    if (class(fit) == "lm") {
        contrast1 <- contrast[!is.na(coef(fit))]
    } else {
        contrast1 <- contrast[!is.na(fixef(fit))]
    }

	return(contrast1)
}


##================================
## .make.contrast.subject.quantification: 
## label-based/label-free; multiple features
## all fixed subject and run
##================================

.make.contrast.subject.quantification <- function(fit, contrast.matrix, sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
    } else {
        coef_name <- names(fixef(fit))
    }
    
    ## when there are missing value in endogenous, there are error, because the 
    ## number of group_original and fitted group are different
    #contrastList<-unique(sub1[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
    #contrastGroup<-rep(0,nlevels(sub1$GROUP_ORIGINAL))
    #contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP_ORIGINAL"])]<-1
    
    ## for label-based
    if (nlevels(sub1$LABEL) == 2) {
        contrastList <- unique(sub1[, c("GROUP", "SUBJECT_ORIGINAL", 
                                        "SUBJECT_NESTED")])
        ## remove GROUP==0
        contrastList <- contrastList[contrastList$GROUP != "0", ]
        ## remove '0' group
		contrastList$GROUP <- factor(contrastList$GROUP)
		## remove '0' group
		contrastList$SUBJECT_ORIGINAL <- factor(contrastList$SUBJECT_ORIGINAL)
		contrastGroup <- rep(0, nlevels(sub1$GROUP) - 1)
		contrastGroup[as.numeric(contrastList[contrast.matrix == 1, 
		                                      "GROUP"])] <- 1
		contrastSubjectOriginal <- rep(0, nlevels(sub1$SUBJECT_ORIGINAL))
		contrastSubjectOriginal[as.numeric(
		    contrastList[contrast.matrix == 1, "SUBJECT_ORIGINAL"])] <- 1
    } else {
        ## for label-free
        contrastList <- unique(sub1[, c("GROUP", "SUBJECT_ORIGINAL")])
        contrastGroup <- rep(0, nlevels(sub1$GROUP))
        contrastGroup[as.numeric(contrastList[contrast.matrix == 1, 
                                              "GROUP"])] <- 1
    }
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    if (length(temp) > 0) {
        intercept_c <- rep(1, length(temp))
        names(intercept_c) <- temp
    } else {
        intercept_c <- NULL
    }
    
    ## feature
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    
    # if(length(temp)>0){
    # 	if(class(fit)=="lm"){
    # 		feature_c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
    # 	}else{
    # 		### need to fix fit$omdel
    # 		tempfeature<-eval(getCall(fit)$data)
    # 		tempfeature$FEATURE<-factor(tempfeature$FEATURE)
    # 		feature_c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
    # 	}
    # 	names(feature_c)<-temp
    # }else{
    # 	feature_c<-NULL
    # }
    
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("FEATURE", "SUBJECT_NESTED")])
		tempSub1 <- xtabs(~ SUBJECT_NESTED + FEATURE, data=tempSub)
		if (nlevels(sub1$LABEL) == 2) {
		    tempSub1 <- tempSub1[-1, ]
			tempSub2 <- tempSub1[contrast.matrix == 1, ]
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
			tempSub2 <- tempSub1[contrast.matrix == 1, ]
		}
		feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
		names(feature_c) <- temp
    } else {
        feature_c <- NULL
    }
    
    ## subject_nested
    temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
    if (length(temp) > 0) {
        if (nlevels(sub1$LABEL) == 2) {
            subjectNested_c <- contrast.matrix
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    subjectNested_c <- contrast.matrix[-1]
		}
        names(subjectNested_c) <- temp
	} else {
	    subjectNested_c <- NULL
	}
    
    ## subject_original
    temp <- coef_name[setdiff(grep("SUBJECT_ORIGINAL", coef_name), 
                              grep(":", coef_name))]
    if (length(temp) > 0) {
        subjectOrig_c <- contrastSubjectOriginal[-1]
        names(subjectOrig_c) <- temp
    } else {
        subjectOrig_c <- NULL
    }
    
    ## group
    temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        if (nlevels(sub1$LABEL) == 2) {
            group_c <- contrastGroup
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    group_c <- contrastGroup[-1]
		}
        names(group_c) <- temp
    } else {
        group_c <- NULL
    }
    
    ## run
    temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        
        # if(class(fit)=="lm"){
        # 	run_c<-rep(1/nlevels(fit$model$RUN),length(temp))
        # }else{
        # 	### need to fix fit$omdel
        # 	run_c<-rep(1/nlevels(eval(getCall(fit)$data)$RUN),length(temp))
        # }
        
        run_c <- rep(1 / nlevels(sub1$RUN), length(temp))
        names(run_c) <- temp
    } else {
        run_c<-NULL
    }
    
    ## feature by group
    temp <- coef_name[intersect(grep("GROUP", coef_name), 
                                grep("FEATURE", coef_name))]
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("GROUP", "FEATURE")])
		tempSub1 <- xtabs(~ GROUP, data=tempSub)
		tempSub2 <- tempSub1[-1]
		temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\:")), nrow=2))
		temp2 <- as.vector(xtabs(~ temp1[, 2]))
		if (nlevels(sub1$LABEL) == 2) {
		    gf_c <- rep(as.numeric(contrastGroup) / tempSub2, temp2)
		} else if (nlevels(sub1$LABEL) == 1) {
		    # label-free
			gf_c <- rep(contrastGroup[-1] / tempSub2, temp2)
		}
		names(gf_c) <- temp
    } else {
        gf_c <- NULL
    }
    
    ## run by feature
    temp <- coef_name[intersect(grep("RUN", coef_name), 
                                grep("FEATURE", coef_name))]
    tempSub <- dim(unique(sub1[, c("RUN", "FEATURE")]))[1]
    if (length(temp) > 0) {
        rf_c <- rep(1 / tempSub, length(temp))
        names(rf_c) <- temp
    } else {
        rf_c <- NULL
    }

	contrast <- c(intercept_c, feature_c, subjectNested_c, subjectOrig_c, 
	              group_c, run_c, gf_c, rf_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##================================
## .make.contrast.subject.quantification.reference: 
## label-based/label-free; single/multiple features
## all fixed subject and run
##================================

.make.contrast.subject.quantification.reference <- function(fit, 
                                                            contrast.matrix, 
                                                            sub1) {
    
    if (class(fit) == "lm") {
        coef_name <- names(coef(fit))
    } else {
        coef_name <- names(fixef(fit))
    }
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(1, length(temp))
	names(intercept_c) <- temp
	if (length(temp) == 0) {
	    intercept_c <- NULL
	}
	
	## feature
	temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
	
	# if(length(temp)>0){
	# 	if(class(fit)=="lm"){
	# 		feature_c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
	# 	}else{
	# 		### need to fix fit$omdel
	# 		tempfeature<-eval(getCall(fit)$data)
	# 		tempfeature$FEATURE<-factor(tempfeature$FEATURE)
	# 		feature_c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
	# 	}
	# 	names(feature_c)<-temp
	# }else{
	# 	feature_c<-NULL
	# }
	
	if (length(temp) > 0) {
	    tempSub <- unique(sub1[, c("FEATURE", "SUBJECT_NESTED")])
	    tempSub1 <- xtabs(~ SUBJECT_NESTED + FEATURE, data=tempSub)
	    if (nlevels(sub1$LABEL) == 2) {
	        tempSub1 <- tempSub1[-1, ]
			tempSub2 <- tempSub1[contrast.matrix == 1, ]
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    tempSub2 <- tempSub1[contrast.matrix == 1, ]
		}
	    feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
		names(feature_c) <- temp
	} else {
	    feature_c <- NULL
	}
	
	## subject_nested
	temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
	if (length(temp) > 0) {
	    subjectNested_c <- rep(0, length(temp))
	    names(subjectNested_c) <- temp
	} else {
	    subjectNested_c <- NULL
	}
	
	## subject_original
	temp <- coef_name[setdiff(grep("SUBJECT_ORIGINAL", coef_name), 
	                          grep(":", coef_name))]
	if (length(temp) > 0) {
	    subjectOrig_c <- rep(0, length(temp))
	    names(subjectOrig_c) <- temp
	} else {
	    subjectOrig_c <- NULL
	}
	
	## group
	temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
	if (length(temp) > 0) {
	    group_c <- rep(0, length(temp))
	    names(group_c) <- temp
	} else {
	    group_c <- NULL
	}
	
	## run
	temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
	if (length(temp) > 0) {
	    
	    # if(class(fit)=="lm"){
	    # 	run_c<-rep(1/nlevels(fit$model$RUN),length(temp))
	    # }else{
	    # 	### need to fix fit$omdel
	    # 	run_c<-rep(1/nlevels(eval(getCall(fit)$data)$RUN),length(temp))
	    # }
	    
	    run_c <- rep(1 / nlevels(sub1$RUN), length(temp))
	    names(run_c) <- temp
	} else {
	    run_c <- NULL
	}
	
	## feature by group
	temp <- coef_name[intersect(grep("GROUP", coef_name), 
	                            grep("FEATURE", coef_name))]
	if (length(temp) > 0) {
	    gf_c <- rep(0, length(temp))	
		names(gf_c) <- temp
	} else {
	    gf_c <- NULL
	}
	
	## run by feature
	temp <- coef_name[intersect(grep("RUN", coef_name), 
	                            grep("FEATURE", coef_name))]
	tempSub <-dim(unique(sub1[, c("RUN", "FEATURE")]))[1]
	if (length(temp) > 0) {
	    rf_c <- rep(1 / tempSub, length(temp))
		names(rf_c) <- temp
	} else {
	    rf_c <- NULL
	}
	
	contrast <- c(intercept_c, feature_c, subjectNested_c, subjectOrig_c, 
	              run_c, gf_c, rf_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}

##================================
## label-free, single
##================================
.make.contrast.free.single <- function(fit, contrast.matrix, sub1) {
    
    if (class(fit) == "lm") {
	    coef_name <- names(coef(fit))
	} else {
	    coef_name <- names(fixef(fit))
	}
    
    ## change contrast.matrix without Group1
    # cons=1
    # if(contrast.matrix[1]==0) contrast.free<-contrast.matrix[-1]
    # if(contrast.matrix[1]<0){ 
    # 	contrast.free<-contrast.matrix[-1]/abs(contrast.matrix[1])
    # 	cons=abs(contrast.matrix[1])
    # 	}
    # if(contrast.matrix[1]>0){ 
    # 	contrast.free<-contrast.matrix[-1]*(-1)/contrast.matrix[1]
    # 	cons=contrast.matrix[1]*(-1)
    # 	}
    # if(class(fit)=="lm"){
    # 	coef_name<-names(coef(fit))
    # }else{
    # 	coef_name<-names(fixef(fit))
    # }
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    intercept_c <- rep(0, length(temp))
    names(intercept_c) <- temp
    if (length(temp) == 0) {
        intercept_c <- NULL
    }
    
    ## subject
    temp <- coef_name[setdiff(grep("SUBJECT", coef_name), 
                              grep(":|NESTED", coef_name))]
    if (length(temp) > 0) {
        
        # subject_c<-rep(0,length(temp))
        # names(subject_c)<-temp
        
        tempdata <- fit$model
        group_levels <- levels(tempdata$GROUP)
        labels <- paste("GROUP", group_levels, sep="")
        patients <- NULL
        for (i in 1:length(group_levels)) {
            sub <- tempdata[tempdata$GROUP == group_levels[i], ]
            sub_patients <- cbind(
                GROUP=paste("GROUP", group_levels[i], sep=""), 
                SUBJECT=paste("SUBJECT", as.character(group_levels(sub$SUBJECT)), sep=""), 
                Value=contrast.matrix[as.numeric(as.character(group_levels[i]))])
            patients <- data.frame(rbind(patients, sub_patients))
        }
        patient_count <- tapply(patients$SUBJECT, patients$GROUP, 
                                function(x) length(unique(x)))
        patient_seq <- rep(0, length(temp))
        for (i in 1:length(as.character(patients$SUBJECT))) {
            match <- any(temp == as.character(patients$SUBJECT)[i])
            if (match & as.numeric(as.character(patients$Value[i])) != 0) {
                res <- temp == as.character(patients$SUBJECT)[i]
                index <- which(res == TRUE)
                group <- as.character(patients[i, ]$GROUP)
                count <- as.numeric(patient_count[names(patient_count) == group])
                value <- as.numeric(as.character(patients[i, ]$Value))
                patient_value <- c(rep(0, index-1), 
                                   value / count, 
                                   rep(0, length(temp) - index))
            } else {
                patient_value <- rep(0, length(temp))
            }
            patient_seq <- patient_value + patient_seq
        }
        subject_c <- patient_seq
        names(subject_c) <- temp
    } else if (length(temp) == 0) {
        subject_c <- NULL
    }
    
    ## group: different from labeled
    temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
    ## when there are some groups which are all missing
	tempSub <- as.numeric(as.character(levels(sub1[, c("GROUP")])))
	tempcontrast <- contrast.matrix[tempSub]
	## for label-free, need to remove first
	group_c <- tempcontrast[-1] 
	names(group_c) <- temp
	if (length(temp) == 0) {
	    group_c<-NULL
	}
	
	## subject_nested
	temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
	if (length(temp) > 0) {
	    temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\.")), nrow=2))
		temp2 <- as.vector(xtabs(~ temp1[, 1]))
		tempdata <- fit$model
		group_levels <- levels(tempdata$GROUP)
		sub_contrast <- contrast.matrix[as.numeric(as.character(group_levels))]
		
		## the base is alway be the first SUBJECT_NESTED
		## (SUBJECT_NESTED1.1)
		temp3 <- temp2
		if (length(temp2) == length(sub_contrast)) {
		    ## this line first, otherwise length of temp3 >1
		    temp3[1] <- temp2[1] + 1
		} else {
		    temp3 <- c(1, temp3)
		}
		
		## in case of unequal sample per group, wrong
		# subjectNested_c<-rep(contrast.matrix/(temp3),temp2) 
		
		subjectNested_c <- rep(sub_contrast / temp3, temp3)[-1]
		names(subjectNested_c) <- temp
	} else if (length(temp) == 0) {
	    subjectNested_c <- NULL
	}
	
	## subject by group: only for time-course - SUBJECT and GROUP (order) even GROUP:SUBJECT in model
	temp <- coef_name[intersect(grep("SUBJECT", coef_name), 
	                            grep("GROUP", coef_name))]
	if (length(temp) > 0) {
	    
	    # subject_c<-rep(0,length(temp))
	    # names(subject_c)<-temp
	    
	    tempdata <- fit$model
	    group_levels <- levels(tempdata$GROUP)
	    labels <- paste("GROUP", group_levels, sep="")
	    patients <- NULL
	    for (i in 1:length(group_levels)) {
	        sub <- tempdata[tempdata$GROUP == group_levels[i], ]
	        sub_patients <- cbind(
	            GROUP=paste("GROUP", group_levels[i], sep=""), 
	            SUBJECT=paste("SUBJECT", as.character(levels(sub$SUBJECT)), sep=""), 
	            Value=contrast.matrix[as.numeric(as.character(group_levels[i]))])
	        patients <- data.frame(rbind(patients, sub_patients))
	    }
	    patient_count <- tapply(patients$SUBJECT, patients$GROUP, 
	                            function(x) length(unique(x)))
	    interaction_seq <- rep(0, length(temp))
	    interaction_labels <- paste(as.character(patients$GROUP), 
	                                as.character(patients$SUBJECT), sep=":")
	    for (i in 1:length(as.character(patients$SUBJECT))) {
	        match <- any(temp == interaction_labels[i])
	        if (match & as.numeric(as.character(patients$Value[i])) != 0) {
	            res <- temp == interaction_labels[i]
	            index <- which(res == TRUE)
	            group <- as.character(patients[i, ]$GROUP)
	            count <- as.numeric(patient_count[names(patient_count) == group])
	            value <- as.numeric(as.character(patients[i, ]$Value))
	            interaction_value <- c(rep(0, index - 1), 
	                                   value / count, 
	                                   rep(0, length(temp) - index))
	        } else {
	            interaction_value <- rep(0, length(temp))
	        }
	        interaction_seq <- interaction_value + interaction_seq
	    }
	    gs_c <- interaction_seq
	    names(gs_c) <- temp
	} else if (length(temp) == 0) {
	    gs_c <- NULL
	}
	
	## combine all
	contrast <- c(intercept_c, group_c, subjectNested_c, subject_c, gs_c)
	if (class(fit) == "lm") {
	    contrast1 <- contrast[!is.na(coef(fit))]
	} else {
	    contrast1 <- contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##################################
##################################
########### estimate   ###########
##################################
##################################	
.getParameterFixed <- function(obj) {
    temp1 <- summary.lm(obj)
    cf <- temp1$coefficients
    vcv <- temp1$cov.unscaled * temp1$sigma ^ 2
    ## TODO (): for unbalanced case, variance is weighted by degree of freedom	
	df <- obj$df.residual
	parameter <- list(cf=cf, vcv=vcv, df=df)
	
	return(parameter)
}

.getParameterRandom <- function(obj, df.full) {
	cf <- as.matrix(fixef(obj))
	vcv <- as.matrix(vcov(obj))
	df <- df.full
	parameter <- list(cf=cf, vcv=vcv, df=df)
	
	return(parameter)
}

.estimableFixedRandom <- function(parameter, cm) {
	cm <- matrix(cm, nrow=1)
	ct <- cm %*% parameter$cf[, 1]
	vc <- sqrt(diag(cm %*% parameter$vcv %*% t(cm)))
	prob <- 2 * (1 - pt(abs(ct / vc), parameter$df))
	result <- cbind(est=ct, stderr=vc, t=ct / vc, df=parameter$df, prob=prob)
	colnames(result) <- c("logFC", "SE", "Tvalue", "DF", "pvalue")
	
	return(result)
}

##########################################################################################

.estimableFixedQuantification <- function(cf, cm) {
    cm <- matrix(cm, nrow=1)
    ct <- cm %*% cf[, 1]
    result <- cbind(est=ct)
    colnames(result) <- c("log-intensities")
    
    return(result)
}

##########################################################################################

.estimableRandomQuantification <- function(cf, cm) {
    cm <- matrix(cm, nrow=1)
    ct <- cm %*% cf
    result <- cbind(est=ct)
    colnames(result) <- c("log-intensities")
    
    return(result)
}

##########################################################################################

.estimableFixedQuantificationSurvival <- function(cf, cm) {
    cm <- matrix(cm, nrow=1)
    ct <- cm %*% cf
    result <- cbind(est=ct)
    colnames(result) <- c("log-intensities")
    
    return(result)
}


##########################################################################################

.make.contrast.run.quantification.Survival <- function(fit, contrast.matrix, 
                                                       sub1, labeled) {
    
    coef_name <- names(fit$coefficients)
    
    ## intercept
    temp <- coef_name[grep("Intercept", coef_name)]
    if (length(temp) > 0) {
        intercept_c <- rep(1, length(temp))
        names(intercept_c) <- temp
    } else {
        intercept_c <- NULL
    }
    
    ## feature
    temp <- coef_name[setdiff(grep("FEATURE", coef_name), grep(":", coef_name))]
    if (length(temp) > 0) {
        tempSub <- unique(sub1[, c("FEATURE", "RUN")])
		tempSub1 <- xtabs(~ RUN + FEATURE, data=tempSub)
		tempSub2 <- tempSub1[contrast.matrix == 1, ]
		feature_c <- as.numeric(tempSub2[-1]) / sum(tempSub2)
		names(feature_c) <- temp
    } else {
        feature_c <- NULL
    }
    
    ## run: different with other quantification - first try
    if (!labeled) {
        ## label-free
        temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
        if (length(temp) > 0) {
            run_c <- contrast.matrix[-1]
            names(run_c) <- temp
        } else {
            run_c <- NULL
        }
    } else {
        ## label-based
        temp <- coef_name[setdiff(grep("RUN", coef_name), grep(":", coef_name))]
        if (length(temp) > 0) {
            run_c <- rep(1 / nlevels(sub1$RUN), length(temp))
            names(run_c) <- temp
        } else {
            run_c <- NULL
        }
    }
    
    ## ref
    temp <- coef_name[grep("ref", coef_name)]
    if (length(temp) > 0) {
        ref_levels <- levels(sub1$ref)
        ref_c <- contrast.matrix[1:(length(ref_levels) - 1)]
        names(ref_c) <- temp
    } else {
        ref_c <- NULL
    }
    
    ## run by feature
    temp <- coef_name[intersect(grep("RUN", coef_name), 
                                grep("FEATURE", coef_name))]
    tempSub <- dim(unique(sub1[, c("RUN", "FEATURE")]))[1]
    if (length(temp) > 0) {
        rf_c <- rep(1 / tempSub, length(temp))
		names(rf_c) <- temp
    } else {
        rf_c <- NULL
    }
    
    ## subject_nested
    temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
    if (length(temp) > 0) {
        if (nlevels(sub1$LABEL) == 2) {
            subjectNested_c <- contrast.matrix
		} else if (nlevels(sub1$LABEL) == 1) {
		    ## label-free
		    subjectNested_c <- contrast.matrix[-1]
		}
        names(subjectNested_c) <- temp
    } else {
        subjectNested_c <- NULL
    }
    
    contrast <- c(intercept_c, feature_c, run_c, ref_c, rf_c, subjectNested_c)
    contrast1 <- contrast[!is.na(fit$coefficients)]
    
    return(contrast1)
}



##########################################################################################

.iter.wls.fit.model <- function(data, fit, nrepeats) {
    for (i in 1:nrepeats) {
        if (i == 1) {
            ## lm or lmer
            if (class(fit) == "lm") {
                abs.resids <- data.frame(abs.resids=abs(fit$residuals))
                fitted <- data.frame(fitted=fit$fitted.values)
			} else {
			    abs.resids <- data.frame(abs.resids=abs(resid(fit)))
			    fitted <- data.frame(fitted=fitted(fit))
			}
			data <- data.frame(data, "abs.resids"=abs.resids, "fitted"=fitted)
		}
		fit.loess <- loess(abs.resids ~ fitted, data=data)
		loess.fitted <- data.frame(loess.fitted=fitted(fit.loess))
		data <- data.frame(data, "loess.fitted"=loess.fitted)
		
		## loess fitted valuaes are predicted sd
		data$weight <- 1 / (data$loess.fitted ^ 2)
		data <- data[, -which(colnames(data) %in% "abs.resids")]
		
		## re-fit using weight
		if (class(fit) == "lm") {
		    wls.fit <- lm(formula(fit), data=data, weights=weight)
		} else {
		    wls.fit <- lmer(formula(fit), data=data, weights=weight)
		}
		
		## lm or lmer
		if (class(wls.fit) == "lm") {
		    abs.resids <- data.frame(abs.resids=abs(wls.fit$residuals))
		} else {
		    abs.resids <- data.frame(abs.resids=abs(resid(wls.fit)))
		}
		data <- data.frame(data, "abs.resids"=abs.resids)
		
		data <- data[, -which(colnames(data) %in% c("loess.fitted", "weight"))]
    }
    
    return(wls.fit)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
#' @importFrom grid viewport grid.newpage pushViewport grid.layout
#' 
.multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    
    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    
    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                         ncol=cols, nrow=ceiling(numPlots / cols))
    }
    if (numPlots == 1) {
        print(plots[[1]])
    } else {
        ## Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        ## Make each plot, in the correct location
        for (i in 1:numPlots) {
            ## Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind=TRUE))
            print(plots[[i]], vp=viewport(layout.pos.row=matchidx$row,
                                          layout.pos.col=matchidx$col))
        }
    }
}

