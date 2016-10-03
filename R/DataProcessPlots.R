
#############################################
## dataProcessPlots
#############################################
#' @export
#' @import ggplot2 
#' @importFrom graphics axis image legend mtext par plot.new title
#' @importFrom grDevices dev.off hcl pdf

dataProcessPlots <- function(data=data,
							type=type,
							featureName="Transition",
							ylimUp=FALSE,
							ylimDown=FALSE,
							scale=FALSE,
							interval="CI",
							x.axis.size=10,
							y.axis.size=10,
							text.size=4,
							text.angle=0,
							legend.size=7,
							dot.size.profile=2,
							dot.size.condition=3,
							width=10,
							height=10, 
							which.Protein="all", 
							originalPlot=TRUE, 
							summaryPlot=TRUE,
							save_condition_plot_result=FALSE, 
							address="") {
	
	datafeature <- data$ProcessedData
	datarun <- data$RunlevelData
  
	datafeature$PROTEIN <- factor(datafeature$PROTEIN)	
	datarun$Protein <- factor(datarun$Protein)	
  
	if (!is.element("SUBJECT_NESTED", colnames(datafeature))) {
		stop("Input for dataProcessPlots function should be processed by dataProcess function previously. Please use 'dataProcess' function first.")
	}
  
	if (length(setdiff(toupper(type), c(toupper("ProfilePlot"), toupper("QCPlot"), toupper("ConditionPlot")))) != 0) {
		stop(paste("Input for type=", type, ". However,'type' should be one of \"ProfilePlot\", \"QCPlot\",\"ConditionPlot\".", sep=""))
	}
  

	## Profile plot ##
	## --------------- 
	if (toupper(type) == "PROFILEPLOT") {
    
		## choose Proteins or not
    if (which.Protein != "all") {
			## check which.Protein is name of Protein
			if (is.character(which.Protein)) {
        		temp.name <- which.Protein
        
        		## message if name of Protein is wrong.
        		if (length(setdiff(temp.name,unique(datafeature$PROTEIN)))>0) {
          			stop(paste("Please check protein name. Data set does not have this protein. -", paste(temp.name, collapse=", "), sep=" "))
          		}
      		}
      
			## check which.Protein is order number of Protein
			if (is.numeric(which.Protein)) {
				temp.name <- levels(datafeature$PROTEIN)[which.Protein]
        
				## message if name of Protein is wrong.
				if (length(levels(datafeature$PROTEIN)) < max(which.Protein)) {
					stop(paste("Please check your selection of proteins. There are ", length(levels(datafeature$PROTEIN))," proteins in this dataset.", sep=" "))
				}
			}
      
				## use only assigned proteins
				datafeature <- datafeature[which(datafeature$PROTEIN %in% temp.name), ]
				datafeature$PROTEIN <- factor(datafeature$PROTEIN)
      
				datarun <- datarun[which(datarun$Protein %in% temp.name), ]
				datarun$PROTEIN <- factor(datarun$Protein)
		}
    
    ## assign upper or lower limit
    # MC, 2016/04/21, default upper limit is maximum log2(intensity) after normalization+3, then round-up
    y.limup <- ceiling ( max ( datafeature$ABUNDANCE, na.rm=TRUE ) + 3 )
    
    if (is.numeric(ylimUp)) {
    	y.limup <- ylimUp 
    }
    
    y.limdown=-1
    if (is.numeric(ylimDown)) {
    	y.limdown <- ylimDown 
    }
    
		datafeature <- datafeature[with(datafeature, order(GROUP_ORIGINAL, SUBJECT_ORIGINAL, LABEL)), ]
		datafeature$RUN <- factor(datafeature$RUN, levels=unique(datafeature$RUN), labels=seq(1, length(unique(datafeature$RUN))))
		datafeature$RUN <- as.numeric(datafeature$RUN)
		tempGroupName <- unique(datafeature[, c("GROUP_ORIGINAL", "RUN")])
    
		groupAxis <- as.numeric(xtabs(~GROUP_ORIGINAL, tempGroupName))
		cumGroupAxis <- cumsum(groupAxis)
		lineNameAxis <- cumGroupAxis[-nlevels(datafeature$GROUP_ORIGINAL)]
    
		groupName <- data.frame(RUN=c(0, lineNameAxis) + groupAxis / 2 + 0.5, ABUNDANCE=rep(y.limup-1, length(groupAxis)), Name=levels(datafeature$GROUP_ORIGINAL))
    
		if (length(unique(datafeature$LABEL)) == 2) {
			datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Reference", "Endogenous"))	
		} else {
			if (unique(datafeature$LABEL) == "L") {
				datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Endogenous"))	
      		}
      		if (unique(datafeature$LABEL) == "H") {
        		datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Reference"))
      		}
    	}
    
		## need to fill in incomplete rows for Runlevel data
		haverun <- FALSE
    
		if (sum(is.element(colnames(datarun), "RUN")) != 0) {
			datamat = dcast( Protein ~ RUN, data=datarun, value.var='LogIntensities', keep=TRUE) 
    
    		datarun = melt(datamat, id.vars=c('Protein'))
     		colnames(datarun)[colnames(datarun) %in% c("variable", "value")] <- c('RUN', 'ABUNDANCE')
		
			haverun <- TRUE
    	}
    
    	## remove the column called 'SuggestToFilter' if there.
    	if (any(is.element(colnames(datafeature),"SuggestToFilter"))) {
    		datafeature$SuggestToFilter <- NULL
    	}
    	
    	## remove the column called 'Fiter.Repro' if there.
    	if (any(is.element(colnames(datafeature),"Filter.Repro"))) {
    		datafeature$Filter.Repro <- NULL
    	}
        
   		## save the plots as pdf or not
    	## If there are the file with the same name, add next numbering at the end of file name
    	
    	## y-axis labeling
      temp <- datafeature[!is.na(datafeature[,"ABUNDANCE"]) & !is.na(datafeature[,"INTENSITY"]), ]
      temp <- temp[1, ]
      temptest <- abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
      if (temptest) {
          yaxis.name <- 'Log2-intensities'
      } else {
          yaxis.name <- 'Log10-intensities'
      }	
        					
    	if (originalPlot) {
    		if (address!=FALSE) {
      			allfiles <- list.files()
      
      			num <- 0
      			filenaming <- paste(address, "ProfilePlot", sep="")
      			finalfile <- paste(address, "ProfilePlot.pdf", sep="")
      
      			while(is.element(finalfile, allfiles)) {
        			num <- num + 1
        			finalfile <- paste(paste(filenaming, num, sep="-"), ".pdf", sep="")
      			}	
      
      			pdf(finalfile, width=width, height=height)
    		}
    
    		for (i in 1:nlevels(datafeature$PROTEIN)) {	
    		  
    		  sub <- datafeature[datafeature$PROTEIN == levels(datafeature$PROTEIN)[i], ]
      		sub$FEATURE <- factor(as.character(sub$FEATURE))	
      		sub$SUBJECT <- factor(sub$SUBJECT)	
      		sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
      		sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)
      		sub$PEPTIDE <- factor(as.character(sub$PEPTIDE))
      
				  ## if all measurements are NA,
				  if (nrow(sub) == sum(is.na(sub$ABUNDANCE))) {
        			message(paste("Can't the Profile plot for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datafeature$PROTEIN)), ") because all measurements are NAs."))
        			next()
     	 		}
      
      		## seq for peptide and transition
      		b <- unique(sub[, c("PEPTIDE", "FEATURE")])
      		b <- b[with(b, order(PEPTIDE, FEATURE)), ] ## add because if there are missing value, orders are different.
      
      		temp1 <- xtabs(~b[, 1])
      		ss <- NULL
      		s <- NULL
      
      		for (j in 1:length(temp1)) {
        		temp3 <- rep(j, temp1[j])
        		s <- c(s, temp3)
        		temp2 <- seq(1, temp1[j])
        		ss <- c(ss, temp2)	
      		}
      			
				  ## for annotation of condition
      		groupNametemp <- data.frame(groupName, FEATURE=unique(sub$FEATURE)[1], PEPTIDE=unique(sub$PEPTIDE)[1])

				  if (toupper(featureName) == "TRANSITION") {
				    
				    if ( any(is.element(colnames(sub), "censored")) ) {
              
				      sub$censored <- factor(sub$censored, levels=c('FALSE', 'TRUE'))
				      
       	 			## 1st plot for original plot
        			ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='FEATURE', linetype='FEATURE'), data=sub)+
        			  facet_grid(~LABEL)+
        			  geom_line(size=0.5)+
        			  geom_point(aes_string(x='RUN', y='ABUNDANCE', color='FEATURE', linetype='FEATURE', shape='censored'), data=sub, 
        			             size=dot.size.profile)+
        			  scale_colour_manual(values=s)+
        			  scale_linetype_manual(values=ss)+
        			  scale_shape_manual(values=c(16, 1), labels=c("Detected data", "Censored missing data"))+
        			  scale_x_continuous('MS runs', breaks=cumGroupAxis)+
        			  scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
        			  geom_vline(xintercept=lineNameAxis + 0.5, colour="grey", linetype="longdash")+
        			  labs(title=unique(sub$PROTEIN))+
        			  geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
        			  theme(
        			    panel.background=element_rect(fill='white', colour="black"),
        			    legend.key=element_rect(fill='white', colour='white'),
        			    panel.grid.minor = element_blank(),
        			    strip.background=element_rect(fill='gray95'),
        			    strip.text.x=element_text(colour=c("#00B0F6"), size=14),
        			    axis.text.x=element_text(size=x.axis.size, colour="black"),
        			    axis.text.y=element_text(size=y.axis.size, colour="black"),
        			    axis.ticks=element_line(colour="black"),
        			    axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
        			    axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
        			    title=element_text(size=x.axis.size+8, vjust=1.5),
        			    legend.position="top",
        			    legend.text=element_text(size=legend.size)
        			    )+
        			  guides(color=guide_legend(title=paste("# peptide:", nlevels(sub$PEPTIDE)), 
        			                            title.theme = element_text(size=13, angle=0),
        			                            ncol=3),
        			         linetype=guide_legend(title=paste("# peptide:", nlevels(sub$PEPTIDE)), 
        			                               title.theme = element_text(size=13, angle=0),
        			                               ncol=3),
        			         shape=guide_legend(title=NULL,
        			                            label.theme = element_text(size=11, angle=0)))
        			  
				    } else {
				      ## 1st plot for original plot
				      ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='FEATURE',linetype='FEATURE'), data=sub)+
				        facet_grid(~LABEL)+
				        geom_point(size=dot.size.profile)+
				        geom_line(size=0.5)+
				        scale_colour_manual(values=s)+
				        scale_linetype_manual(values=ss)+
				        scale_x_continuous('MS runs', breaks=cumGroupAxis)+
				        scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
				        geom_vline(xintercept=lineNameAxis + 0.5, colour="grey", linetype="longdash")+
				        labs(title=unique(sub$PROTEIN))+
				        geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
				        theme(
				          panel.background=element_rect(fill='white', colour="black"),
				          legend.key=element_rect(fill='white', colour='white'),
				          panel.grid.minor = element_blank(),
				          strip.background=element_rect(fill='gray95'),
				          strip.text.x=element_text(colour=c("#00B0F6"), size=14),
				          axis.text.x=element_text(size=x.axis.size, colour="black"),
				          axis.text.y=element_text(size=y.axis.size, colour="black"),
				          axis.ticks=element_line(colour="black"),
				          axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
				          axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
				          title=element_text(size=x.axis.size+8, vjust=1.5),
				          legend.position="top",
				          legend.text=element_text(size=legend.size)
				        )+
				        guides(color=guide_legend(title=paste("# peptide:", nlevels(sub$PEPTIDE)), 
				                                  title.theme = element_text(size=13, angle=0),
				                                  ncol=3), 
				               linetype=guide_legend(title=paste("# peptide:", nlevels(sub$PEPTIDE)), 
				                                     title.theme = element_text(size=13, angle=0),
				                                     ncol=3))
				    }
        
					  print(ptemp)
					
					  message(paste("Drew the Profile plot for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datafeature$PROTEIN)), ")"))
				  }
      
      		if (toupper(featureName) == "PEPTIDE") {
      		  
      		  if ( any(is.element(colnames(sub), "censored")) ) {
      		    
      		    sub$censored <- factor(sub$censored, levels=c('FALSE', 'TRUE'))
      		    
      		    ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE', linetype='FEATURE'), data=sub)+
      		      facet_grid(~LABEL)+
      		      geom_line(size=0.5)+
      		      geom_point(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE', linetype='FEATURE', shape='censored'), data=sub, 
      		                 size=dot.size.profile)+
      		      scale_colour_manual(values=unique(s))+ ## unique(s) ??
      		      scale_linetype_manual(values=ss, guide="none")+
      		      scale_shape_manual(values=c(16, 1), labels=c("Detected data", "Censored missing data"))+
      		      scale_x_continuous('MS runs', breaks=cumGroupAxis)+
      		      scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
      		      geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
      		      labs(title=unique(sub$PROTEIN))+
      		      geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
      		      theme(
      		        panel.background=element_rect(fill='white', colour="black"),
      		        legend.key=element_rect(fill='white', colour='white'),
      		        panel.grid.minor = element_blank(),
      		        strip.background=element_rect(fill='gray95'),	
      		        strip.text.x=element_text(colour=c("#00B0F6"), size=14),
      		        axis.text.x=element_text(size=x.axis.size, colour="black"),
      		        axis.text.y=element_text(size=y.axis.size, colour="black"),
      		        axis.ticks=element_line(colour="black"),
      		        axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
      		        axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
      		        title=element_text(size=x.axis.size+8, vjust=1.5),
      		        legend.position="top",
      		        legend.text=element_text(size=legend.size)
      		      )+
      		      guides(color=guide_legend(title=paste("# peptide:", nlevels(sub$PEPTIDE)), 
      		                                title.theme = element_text(size=13, angle=0),
      		                                ncol=3), 
      		             shape=guide_legend(title=NULL,
      		                                label.theme = element_text(size=11, angle=0)))
      		    
      		  } else {
      		    
      		    ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE', linetype='FEATURE'), data=sub)+
      		      facet_grid(~LABEL)+
      		      geom_point(size=dot.size.profile)+
      		      geom_line(size=0.5)+
      		      scale_colour_manual(values=unique(s))+
      		      scale_linetype_manual(values=ss, guide="none")+
      		      scale_x_continuous('MS runs', breaks=cumGroupAxis)+
      		      scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
      		      geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
      		      labs(title=unique(sub$PROTEIN))+
      		      geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
      		      theme(
      		        panel.background=element_rect(fill='white', colour="black"),
      		        legend.key=element_rect(fill='white', colour='white'),
      		        panel.grid.minor = element_blank(),
      		        strip.background=element_rect(fill='gray95'),	
      		        strip.text.x=element_text(colour=c("#00B0F6"), size=14),
      		        axis.text.x=element_text(size=x.axis.size, colour="black"),
      		        axis.text.y=element_text(size=y.axis.size, colour="black"),
      		        axis.ticks=element_line(colour="black"),
      		        axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
      		        axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
      		        title=element_text(size=x.axis.size+8, vjust=1.5),
      		        legend.position="top",
      		        legend.text=element_text(size=legend.size)
      		      )+
      		      guides(color=guide_legend(title=paste("# peptide:", nlevels(sub$PEPTIDE)), 
      		                                title.theme = element_text(size=13, angle=0),
      		                                ncol=3))
      		    
      		  }
        
					  print(ptemp)

        		message(paste("Drew the Profile plot for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datafeature$PROTEIN)), ")"))
        
      		}
      
      		if (toupper(featureName) == "NA") {
      		  
      		  if ( any(is.element(colnames(sub), "censored")) ) {
      		    
      		    sub$censored <- factor(sub$censored, levels=c('FALSE', 'TRUE'))
      		    
      		    ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE', linetype='FEATURE'), data=sub)+
      		      facet_grid(~LABEL)+
      		      geom_line(size=0.5)+
      		      geom_point(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE', linetype='FEATURE', shape='censored'), data=sub, 
      		                 size=dot.size.profile)+
      		      scale_colour_manual(values=unique(s), guide="none")+ 
      		      scale_linetype_manual(values=ss, guide="none")+
      		      scale_shape_manual(values=c(16, 1), labels=c("Detected data", "Censored missing data"))+
      		      scale_x_continuous('MS runs', breaks=cumGroupAxis)+
      		      scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
      		      geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
      		      labs(title=unique(sub$PROTEIN))+
      		      geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
      		      theme(
      		        panel.background=element_rect(fill='white', colour="black"),
      		        legend.key=element_rect(fill='white', colour='white'),
      		        panel.grid.minor = element_blank(),
      		        strip.background=element_rect(fill='gray95'),	
      		        strip.text.x=element_text(colour=c("#00B0F6"), size=14),
      		        axis.text.x=element_text(size=x.axis.size, colour="black"),
      		        axis.text.y=element_text(size=y.axis.size, colour="black"),
      		        axis.ticks=element_line(colour="black"),
      		        axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
      		        axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
      		        title=element_text(size=x.axis.size+8, vjust=1.5),
      		        legend.position="top",
      		        legend.text=element_text(size=legend.size)
      		      )+
      		      guides(shape=guide_legend(title=NULL,
      		                                label.theme = element_text(size=11, angle=0)))
      		    
      		  } else {
      		    
      		    ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE', linetype='FEATURE'), data=sub)+
      		      facet_grid(~LABEL)+
      		      geom_point(size=dot.size.profile)+
      		      geom_line(size=0.5)+
      		      scale_colour_manual(values=unique(s), guide="none")+
      		      scale_linetype_manual(values=ss, guide="none")+
      		      scale_x_continuous('MS runs', breaks=cumGroupAxis)+
      		      scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
      		      geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
      		      labs(title=unique(sub$PROTEIN))+
      		      geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
      		      theme(
      		        panel.background=element_rect(fill='white', colour="black"),
      		        legend.key=element_rect(fill='white', colour='white'),
      		        panel.grid.minor = element_blank(),
      		        strip.background=element_rect(fill='gray95'),	
      		        strip.text.x=element_text(colour=c("#00B0F6"), size=14),
      		        axis.text.x=element_text(size=x.axis.size, colour="black"),
      		        axis.text.y=element_text(size=y.axis.size, colour="black"),
      		        axis.ticks=element_line(colour="black"),
      		        axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
      		        axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
      		        title=element_text(size=x.axis.size+8, vjust=1.5),
      		        legend.position="top",
      		        legend.text=element_text(size=legend.size)
      		      )
      		    
      		  }
        
       			print(ptemp)
        	
        		message(paste("Drew the Profile plot for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datafeature$PROTEIN)), ")"))
        
     		 	}
    		} # end-loop for each protein
    		
   	 		if (address!=FALSE) {
   	 			dev.off()
   	 		} 
   	 		
    	} # end original plot
    
		  ## 2st plot for original plot : summary ##
    	## ---------------------------------------
    	
    	if (summaryPlot) {
			  if (address!=FALSE) {
     	 		allfiles <- list.files()
      
     	 		num <- 0
      			filenaming <- paste(address, "ProfilePlot_wSummarization", sep="")
      			finalfile <- paste(address, "ProfilePlot_wSummarization.pdf", sep="")
      
      			while(is.element(finalfile, allfiles)) {
        			num <- num + 1
        			finalfile <- paste(paste(filenaming, num, sep="-"), ".pdf", sep="")
      			}	
      
     	 		pdf(finalfile, width=width, height=height)
    		}

    		for (i in 1:nlevels(datafeature$PROTEIN)) {	
    			
     			sub <- datafeature[datafeature$PROTEIN == levels(datafeature$PROTEIN)[i], ]
      		sub$FEATURE <- factor(as.character(sub$FEATURE))	
      		sub$SUBJECT <- factor(sub$SUBJECT)	
      		sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
      		sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)
      		sub$PEPTIDE <- factor(as.character(sub$PEPTIDE))
      
      		## if all measurements are NA,
     		 	if (nrow(sub) == sum(is.na(sub$ABUNDANCE))) {
       			 	message(paste("Can't the Profile plot for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datafeature$PROTEIN)), ") because all measurements are NAs."))
        			next()
     		 	}
      
      		## seq for peptide and transition
     		 	b <- unique(sub[, c("PEPTIDE", "FEATURE")])
     		 	b <- b[with(b, order(PEPTIDE, FEATURE)), ] ## add because if there are missing value, orders are different.
      
     		 	temp1 <- xtabs(~b[, 1])
      		ss <- NULL
      		s <- NULL
      
      		for( j in 1:length(temp1) ) {
        		temp3 <- rep(j, temp1[j])
       			s <- c(s, temp3)
       			temp2 <- seq(1, temp1[j])
        		ss <- c(ss, temp2)	
     			}
     			
     			## for annotation of condition
      		groupNametemp <- data.frame(groupName, FEATURE=unique(sub$FEATURE)[1], analysis="Run summary")
      
      		if (haverun) {
        		subrun <- datarun[datarun$Protein == levels(datafeature$PROTEIN)[i], ]
			
					  if ( nrow(subrun) != 0 ) {

						  quantrun <- sub[1,]
						  quantrun[, 2:ncol(quantrun)] <- NA
						  quantrun <- quantrun[rep(seq_len(nrow(subrun))), ]
						  
						  quantrun$PROTEIN <- subrun$Protein 
						  quantrun$PEPTIDE <- "Run summary"
						  quantrun$TRANSITION <- "Run summary" 
						  quantrun$FEATURE <- "Run summary" 
						  quantrun$LABEL <- "Endogenous"
						  quantrun$RUN <- subrun$RUN
						  quantrun$ABUNDANCE <- subrun$ABUNDANCE
						  quantrun$METHOD <- 1
					    
            } else { # if there is only one Run measured across all runs, no Run information for linear with censored

						  quantrun <- datafeature[1,]
						  quantrun[, 2:ncol(quantrun)] <- NA

						  quantrun$PROTEIN <- levels(datafeature$PROTEIN)[i]
						  quantrun$PEPTIDE <- "Run summary"
						  quantrun$TRANSITION <- "Run summary" 
						  quantrun$FEATURE <- "Run summary" 
						  quantrun$LABEL <- "Endogenous"
						  quantrun$RUN <- unique(datafeature$RUN)[1]
						  quantrun$ABUNDANCE <- NA
						  quantrun$METHOD <- 1
						  
					  }
			
					  if (any(is.element(colnames(sub), "censored"))) {
						  quantrun$censored <- FALSE
					  }

					  quantrun$analysis <- "Run summary"
					  sub$analysis <- "Processed feature-level data"
					
					  ## if 'Filter' column after feature selection, remove this column in order to match columns with run quantification
					  filter_column <- is.element(colnames(sub), "Filter")
					  if (any(filter_column)) {
						  sub<-sub[, !filter_column]
					  }
					
					  final <- rbind(sub, quantrun)
					  final$analysis <- factor(final$analysis)
					  final$FEATURE <- factor(final$FEATURE)
					  final$RUN <- as.numeric(final$RUN)
					  
					  if (any(is.element(colnames(sub), "censored"))) {
					    
					    final$censored <- factor(final$censored, levels=c('FALSE', 'TRUE'))
					    
					    ptempall <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='analysis', linetype='FEATURE', size='analysis'), data=final)+
					      facet_grid(~LABEL)+
					      geom_line(size=0.5)+
					      geom_point(aes_string(x='RUN', y='ABUNDANCE', color='analysis', linetype='FEATURE', size='analysis', shape='censored'), data=final)+
					      scale_colour_manual(values=c("lightgray", "darkred"))+
					      scale_shape_manual(values=c(16, 1), labels=c("Detected data", "Censored missing data"))+
					      scale_size_manual(values=c(1.7, 2), guide="none")+
					      scale_linetype_manual(values=c(rep(1, times=length(unique(final$FEATURE))-1), 2), guide="none")+
					      scale_x_continuous('MS runs',breaks=cumGroupAxis)+
					      scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
					      geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
					      labs(title=unique(final$PROTEIN))+
					      geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
					      theme(
					        panel.background=element_rect(fill='white', colour="black"),
					        legend.key=element_rect(fill='white', colour='white'),
					        panel.grid.minor = element_blank(),
					        strip.background=element_rect(fill='gray95'),
					        strip.text.x=element_text(colour=c("#00B0F6"), size=14),
					        axis.text.x=element_text(size=x.axis.size, colour="black"),
					        axis.text.y=element_text(size=y.axis.size, colour="black"),
					        axis.ticks=element_line(colour="black"),
					        axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
					        axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
					        title=element_text(size=x.axis.size+8, vjust=1.5),
					        legend.position="top",
					        legend.text=element_text(size=legend.size),
					        legend.title=element_blank()
					      )+
					      guides(color=guide_legend(order=1, title=NULL, label.theme = element_text(size=10, angle=0)),
					             shape=guide_legend(order=2, title=NULL, label.theme = element_text(size=10, angle=0)))
					    
					  } else {
					    
					    ptempall <- ggplot(aes_string(x='RUN', y='ABUNDANCE', color='analysis', linetype='FEATURE', size='analysis'), data=final)+
					      facet_grid(~LABEL)+
					      geom_point(size=dot.size.profile)+
					      geom_line(size=0.5)+
					      scale_colour_manual(values=c("lightgray", "darkred"))+
					      scale_size_manual(values=c(1.7, 2), guide="none")+
					      scale_linetype_manual(values=c(rep(1, times=length(unique(final$FEATURE))-1), 2), guide="none")+
					      scale_x_continuous('MS runs',breaks=cumGroupAxis)+
					      scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
					      geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
					      labs(title=unique(final$PROTEIN))+
					      geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
					      theme(
					        panel.background=element_rect(fill='white', colour="black"),
					        legend.key=element_rect(fill='white', colour='white'),
					        panel.grid.minor = element_blank(),
					        strip.background=element_rect(fill='gray95'),
					        strip.text.x=element_text(colour=c("#00B0F6"), size=14),
					        axis.text.x=element_text(size=x.axis.size, colour="black"),
					        axis.text.y=element_text(size=y.axis.size, colour="black"),
					        axis.ticks=element_line(colour="black"),
					        axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
					        axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
					        title=element_text(size=x.axis.size+8, vjust=1.5),
					        legend.position="top",
					        legend.text=element_text(size=legend.size),
					        legend.title=element_blank()
					      )+
					      guides(color=guide_legend(order=1, title=NULL, label.theme = element_text(size=10, angle=0)))
					    
					   ## draw point again because some red summary dots could be hiden
					   ptempall <- ptempall+geom_point(data=final, aes(x=RUN, y=ABUNDANCE, size=analysis, color=analysis))
					    
					 }
				
					print(ptempall)
			
					message(paste("Drew the Profile plot with summarization for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datafeature$PROTEIN)), ")"))
			
				}

    		} # end-loop for each protein
    		if (address!=FALSE) {
    			dev.off()
    		}
  		}  
	} # end Profile plot	
  

	## QC plot (Quality control plot) ##
	## ---------------------------------
	if (toupper(type) == "QCPLOT") {
    
        ## y-axis labeling
       	temp <- datafeature[!is.na(datafeature[,"ABUNDANCE"]) & !is.na(datafeature[,"INTENSITY"]), ]
       	temp <- temp[1, ]
        temptest <- abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
        if (temptest) {
          	yaxis.name <- 'Log2-intensities'
        } else {
          	yaxis.name <- 'Log10-intensities'
        }	
        
		## save the plots as pdf or not
    	## If there are the file with the same name, add next numbering at the end of file name		
    	if (address!=FALSE) {
      		allfiles <- list.files()
      
      		num <- 0
      		filenaming <- paste(address,"QCPlot", sep="")
      		finalfile <- paste(address,"QCPlot.pdf", sep="")
      
      		while(is.element(finalfile, allfiles)) {
        		num <- num + 1
        		finalfile <- paste(paste(filenaming, num, sep="-"), ".pdf", sep="")
      		}	
      
      		pdf(finalfile, width=width, height=height)
    	}

    	## assign upper or lower limit
	  	# MC, 2016/04/21, default upper limit is maximum log2(intensity) after normalization+3, then round-up
	  	y.limup <- ceiling ( max ( datafeature$ABUNDANCE, na.rm=TRUE ) + 3 )
    
    	if (is.numeric(ylimUp)) {
    		y.limup <- ylimUp 
    	}
    
    	y.limdown=-1
    	if (is.numeric(ylimDown)) {
    		y.limdown <- ylimDown 
    	}
    	
    	## relabel the Run (make it sorted by group first)
    	datafeature <- datafeature[with(datafeature, order(GROUP_ORIGINAL, SUBJECT_ORIGINAL)), ]
    	datafeature$RUN <- factor(datafeature$RUN, levels=unique(datafeature$RUN), labels=seq(1, length(unique(datafeature$RUN))))
    
    	if (length(unique(datafeature$LABEL)) == 2) {
      		datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Reference", "Endogenous"))	
      		label.color <- c("darkseagreen1", "lightblue")
    	} else {
      		if (unique(datafeature$LABEL) == "L") {
        		datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Endogenous"))
        		label.color <- c("lightblue")	
      		}
      		if (unique(datafeature$LABEL) == "H") {
        		datafeature$LABEL <- factor(datafeature$LABEL, labels=c("Reference"))
        		label.color <- c("darkseagreen1")
      		}
    	}	
    
    	tempGroupName <- unique(datafeature[, c("GROUP_ORIGINAL", "RUN")])
    	datafeature <- datafeature[with(datafeature, order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL)),]
    
    	groupAxis <- as.numeric(xtabs(~GROUP_ORIGINAL, tempGroupName))
    	cumGroupAxis <- cumsum(groupAxis)
    	lineNameAxis <- cumGroupAxis[-nlevels(datafeature$GROUP_ORIGINAL)]
    
    	groupName <- data.frame(RUN=c(0, lineNameAxis)+groupAxis / 2 + 0.5, ABUNDANCE=rep(y.limup-1, length(groupAxis)), Name=levels(datafeature$GROUP_ORIGINAL))
    
    	## all protein
    	ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE'), data=datafeature)+
    			facet_grid(~LABEL)+
    			geom_boxplot(aes_string(fill='LABEL'), outlier.shape=1, outlier.size=1.5)+
    			scale_fill_manual(values=label.color, guide="none")+
    			scale_x_discrete('MS runs', breaks=cumGroupAxis)+
    			scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
    			geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
    			labs(title="All proteins")+
    			geom_text(data=groupName, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
    			theme(
      				panel.background=element_rect(fill='white', colour="black"),
      				legend.key=element_rect(fill='white', colour='white'),
      				panel.grid.minor = element_blank(),
      				strip.background=element_rect(fill='gray95'),	
      				strip.text.x=element_text(colour=c("#00B0F6"), size=14),
      				axis.text.x=element_text(size=x.axis.size,colour="black"),
      				axis.text.y=element_text(size=y.axis.size,colour="black"),
      				axis.ticks=element_line(colour="black"),
      				axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
      				axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
      				title=element_text(size=x.axis.size+8, vjust=1.5)
      			)
    
    	print(ptemp)
    
    	message("Drew the Quality Contol plot(boxplot) for all proteins.")
    
    	## each protein
		## choose Proteins or not
    	if (which.Protein != "all") {
      		## check which.Protein is name of Protein
      		if (is.character(which.Protein)) {
        
        		temp.name <- which.Protein
        
        		## message if name of Protein is wrong.
        		if (length(setdiff(temp.name, unique(datafeature$PROTEIN))) > 0) {
          			dev.off()
          			stop(paste("Please check protein name. Data set does not have this protein. -", paste(temp.name, collapse=", "), sep=" "))
        		}
      		}
      
      		## check which.Protein is order number of Protein
      		if (is.numeric(which.Protein)) {
       			temp.name <- levels(datafeature$PROTEIN)[which.Protein]
        
        		## message if name of Protein is wrong.
        		if (length(levels(datafeature$PROTEIN))<max(which.Protein)) {
          			dev.off()
          			stop(paste("Please check your selection of proteins. There are ", length(levels(datafeature$PROTEIN))," proteins in this dataset.",sep=" "))
        		}
      		}
      
      		## use only assigned proteins
      		datafeature <- datafeature[which(datafeature$PROTEIN %in% temp.name), ]
      		datafeature$PROTEIN <- factor(datafeature$PROTEIN)
    	}
    
    	for (i in 1:nlevels(datafeature$PROTEIN)) {	
      		sub <- datafeature[datafeature$PROTEIN==levels(datafeature$PROTEIN)[i], ]
      		subTemp <- sub[!is.na(sub$ABUNDANCE),]
      		sub <- sub[with(sub, order(LABEL,RUN)),]
      
     		## if all measurements are NA,
      		if (nrow(sub)==sum(is.na(sub$ABUNDANCE))) {
        		message(paste("Can't the Quality Control plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),") because all measurements are NAs."))
        		next()
      		}
      
      		ptemp <- ggplot(aes_string(x='RUN', y='ABUNDANCE'), data=sub)+
      				facet_grid(~LABEL)+
      				geom_boxplot(aes_string(fill='LABEL'), outlier.shape=1, outlier.size=1.5)+
      				scale_fill_manual(values=label.color, guide="none")+
      				scale_x_discrete('MS runs', breaks=cumGroupAxis)+
      				scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
      				geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
      				labs(title=unique(sub$PROTEIN))+
      				geom_text(data=groupName, aes(x=RUN, y=ABUNDANCE, label=Name), size=text.size, angle=text.angle, color="black")+
      				theme(
        				panel.background=element_rect(fill='white', colour="black"),
        				legend.key=element_rect(fill='white', colour='white'),
        				panel.grid.minor = element_blank(),
        				strip.background=element_rect(fill='gray95'),	
        				strip.text.x=element_text(colour=c("#00B0F6"), size=14),
        				axis.text.x=element_text(size=x.axis.size, colour="black"),
        				axis.text.y=element_text(size=y.axis.size, colour="black"),
        				axis.ticks=element_line(colour="black"),
        				axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
        				axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
        				title=element_text(size=x.axis.size+8, vjust=1.5)
        			)
      
      		print(ptemp)
      
      		message(paste("Drew the Quality Contol plot(boxplot) for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datafeature$PROTEIN)), ")"))
      
    	} # end-loop
    	if (address!=FALSE) {
    		dev.off()
    	}
  	} # end QC plot	
  
  

  ## Condition plot ##
  ## -----------------
	if (toupper(type) == "CONDITIONPLOT") {
    	
    	colnames(datarun)[colnames(datarun) == "Protein"] <- "PROTEIN"
    	colnames(datarun)[colnames(datarun) == "LogIntensities"] <- "ABUNDANCE"

    	## choose Proteins or not
    	if (which.Protein != "all") {
      		## check which.Protein is name of Protein
      		if (is.character(which.Protein)) {
        
        		temp.name <- which.Protein
        
        		## message if name of Protein is wrong.
        		if (length(setdiff(temp.name, unique(datarun$PROTEIN))) > 0) {
        			stop(paste("Please check protein name. Dataset does not have this protein. -", paste(temp.name, collapse=", "),sep=" "))
        		}
          			
      		}
      
      		## check which.Protein is order number of Protein
      		if (is.numeric(which.Protein)) {
        
        		temp.name <- levels(datarun$PROTEIN)[which.Protein]
        
        		## message if name of Protein is wrong.
        		if (length(levels(datarun$PROTEIN))<max(which.Protein)) {
          			stop(paste("Please check your selection of proteins. There are ", length(levels(datarun$PROTEIN))," proteins in this dataset.",sep=" "))
          		}
      		}
      
      		## use only assigned proteins
      		datarun <- datarun[which(datarun$PROTEIN %in% temp.name), ]
      		datarun$PROTEIN <- factor(datarun$PROTEIN)
    	}
    
    	## save the plots as pdf or not
    	## If there are the file with the same name, add next numbering at the end of file name		
    	if (address!=FALSE) {
      		allfiles <- list.files()
      
      		num <- 0
      		filenaming <- paste(address, "ConditionPlot", sep="")
      		finalfile <- paste(address, "ConditionPlot.pdf", sep="")
      
      		while(is.element(finalfile, allfiles)) {
        		num <- num + 1
        		finalfile <- paste(paste(filenaming, num, sep="-"), ".pdf", sep="")
      		}	
      
      		pdf(finalfile, width=width, height=height)
    	}
    	
    	## save all results
    	resultall <- NULL
    	
    	## y-axis labeling, find log 2 or log 10
      temp <- datafeature[!is.na(datafeature[, "ABUNDANCE"]) & !is.na(datafeature[, "INTENSITY"]),]
      temp <- temp[1,]
      temptest <- abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"]) < abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
      if (temptest) {
        yaxis.name <- 'Log2-intensities'
      } else {
        yaxis.name <- 'Log10-intensities'
      }
    
      for (i in 1:nlevels(datarun$PROTEIN)) {	
      		
      	suball <- NULL
      		
        sub <- datarun[datarun$PROTEIN == levels(datarun$PROTEIN)[i], ]
        sub <- na.omit(sub)	
      	sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
        sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)	
        
        ## if all measurements are NA,
        if (nrow(sub) == sum(is.na(sub$ABUNDANCE))) {
          	message(paste("Can't the Condition plot for ", unique(sub$PROTEIN), "(", i, " of ",length(unique(datarun$PROTEIN)), ") because all measurements are NAs."))
          	next()
        }
        
        	## statistics
        	sub.mean <- aggregate(ABUNDANCE ~ GROUP_ORIGINAL, data=sub, mean, na.rm=TRUE)
        	sub.sd <- aggregate(ABUNDANCE ~ GROUP_ORIGINAL, data=sub, sd)
        	sub.len <- aggregate(ABUNDANCE ~ GROUP_ORIGINAL, data=sub, length)
        		
        	## make the table for result
        	colnames(sub.mean)[colnames(sub.mean) == "ABUNDANCE"] <- "Mean"
        	colnames(sub.sd)[colnames(sub.sd) == "ABUNDANCE"] <- "SD"
        	colnames(sub.len)[colnames(sub.len) == "ABUNDANCE"] <- "numMeasurement"

			suball <- merge(sub.mean, sub.sd, by="GROUP_ORIGINAL")
			suball <- merge(suball, sub.len, by="GROUP_ORIGINAL")

        	if (interval=="CI") {
        		suball$ciw <- qt(0.975, suball$numMeasurement) * suball$SD / sqrt(suball$numMeasurement)
        	}
        	if (interval=="SD") {
        		suball$ciw <- suball$SD
        	}
        
        	if (sum(is.na(suball$ciw)) >= 1) {
        		suball$ciw[is.na(suball$ciw)] <- 0
        	}
        
        	## assign upper or lower limit
        	y.limup <- ceiling(max(suball$Mean + suball$ciw))
        	if (is.numeric(ylimUp)) {
        		y.limup <- ylimUp 
        	}
        
        	y.limdown <- floor(min(suball$Mean - suball$ciw))
        	if (is.numeric(ylimDown)) {
        		y.limdown <- ylimDown 
        	}
        	
        	## re-order (1, 10, 2, 3, -> 1, 2, 3, ... , 10)
        	suball <- suball[order(suball$GROUP_ORIGINAL), ]
        	suball <- data.frame(Protein=unique(sub$PROTEIN), suball)
        	resultall <- rbind(resultall, suball)
        	
        	if (!scale) {  ## scale: false
          
          		## reformat as data.frame
          		#tempsummary <- data.frame(Label=unique(sub$GROUP_ORIGINAL), mean=as.vector(sub.mean), ciw=as.vector(ciw))
          		tempsummary <- suball
          		colnames(tempsummary)[colnames(tempsummary) == "GROUP_ORIGINAL"] <- "Label"
          
          		ptemp <- ggplot(aes_string(x='Label', y='Mean'), data=tempsummary)+
          				geom_errorbar(aes(ymax = Mean + ciw, ymin= Mean - ciw), data=tempsummary, width=0.1, colour="red")+
          				geom_point(size = dot.size.condition, colour = "darkred")+
          				scale_x_discrete('Condition')+
          				scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
          				geom_hline(yintercept = 0, linetype = "twodash", colour = "darkgrey", size = 0.6)+
          				labs(title=unique(sub$PROTEIN))+
          				theme(
            				panel.background=element_rect(fill='white', colour="black"),
            				panel.grid.major.y = element_line(colour="grey95"),
            				panel.grid.minor.y = element_blank(),
            				axis.text.x=element_text(size=x.axis.size, colour="black", angle=text.angle),
            				axis.text.y=element_text(size=y.axis.size, colour="black"),
            				axis.ticks=element_line(colour="black"),
            				axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
            				axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
            				title=element_text(size=x.axis.size+8, vjust=1.5)
            			)
          
        	} else {
        		## scale : true
          		## extract numeric value, don't use levels (because T1,T10,T3,...)
				## reformat as data.frame
          		
          		tempsummary <- suball
          		colnames(tempsummary)[colnames(tempsummary) == "GROUP_ORIGINAL"] <- "Label"
          		tempsummary$Label <- as.numeric(gsub("\\D", "", unique(tempsummary$Label)))
          
          		ptemp <- ggplot(aes_string(x='Label', y='Mean'), data=tempsummary)+
          				geom_errorbar(aes(ymax = Mean + ciw, ymin = Mean - ciw), data=tempsummary, width=0.1, colour="red")+
          				geom_point(size=dot.size.condition, colour="darkred")+
          				scale_x_continuous('Condition', breaks=tempsummary$Label, labels=tempsummary$Label)+
          				scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
          				geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+
          				labs(title=unique(sub$PROTEIN))+
          				theme(
            				panel.background=element_rect(fill='white', colour="black"),
            				panel.grid.major.y = element_line(colour="grey95"),
            				panel.grid.minor.y = element_blank(),
            				axis.text.x=element_text(size=x.axis.size, colour="black", angle=text.angle),
            				axis.text.y=element_text(size=y.axis.size, colour="black"),
            				axis.ticks=element_line(colour="black"),
            				axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
            				axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
            				title=element_text(size=x.axis.size+8, vjust=1.5)
            			)
          
        	}
        
        	print(ptemp)
        
        	message(paste("Drew the condition plot for ", unique(sub$PROTEIN), "(", i, " of ", length(unique(datarun$PROTEIN)), ")"))
        
      	} # end-loop
    
    	if (address != FALSE) {
    		dev.off()
    	}
    	
    	## save the table for condition plot
    	if(save_condition_plot_result){
    		colnames(resultall)[colnames(resultall) == "GROUP_ORIGINAL"] <- 'Condition'
    		
    		if (interval == "CI") {
        		colnames(resultall)[colnames(resultall) == "ciw"] <- '95% CI'
        	}
        	if (interval == "SD") {
        		colnames(resultall)[colnames(resultall) == "ciw"] <- 'SD'
        	}
        	
        	if (address!=FALSE) {
      			allfiles <- list.files()
      
      			num <- 0
      			filenaming <- paste(address, "ConditionPlot_value", sep="")
      			finalfile <- paste(address, "ConditionPlot_value.csv", sep="")
      
      			while(is.element(finalfile, allfiles)) {
        			num <- num + 1
        			finalfile <- paste(paste(filenaming, num, sep="-"), ".csv", sep="")
      			}	
      
      			write.csv(resultall, file=finalfile, row.names=FALSE)
    		}

    	}
    	

  	} # end Condition plot
}


