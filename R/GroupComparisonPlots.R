
#############################################
## groupComparisonPlots
#############################################

#' @export
#' @importFrom gplots heatmap.2
#' @importFrom stats hclust
#' @importFrom ggrepel geom_text_repel
#' @importFrom marray maPalette

groupComparisonPlots <- function(data=data,
				type=type,
				sig=0.05,
				FCcutoff=FALSE,
				logBase.pvalue=10,
				ylimUp=FALSE,
				ylimDown=FALSE,
				xlimUp=FALSE,
				x.axis.size=10,
				y.axis.size=10,
				dot.size=3,
				text.size=4, 
				legend.size=13,
				ProteinName=TRUE,
				numProtein=100, 
				clustering="both", 
				width=10, 
				height=10, 
				which.Comparison="all", 
				address="") {
  
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
    	processout <- as.matrix(read.table(finalfile, header=T, sep="\t"))
  	}	
  
  	processout <- rbind(processout, as.matrix(c(" ", " ", "MSstats - groupComparisonPlots function", " "), ncol=1))
  
  
  	## make upper letter
  	type <- toupper(type)
  
  	if (length(setdiff(type, c("HEATMAP", "VOLCANOPLOT", "COMPARISONPLOT"))) != 0) {
    
    	processout <- rbind(processout, c(paste("Input for type=", type, ". However,'type' should be one of \"Heatmap\", \"VolcanoPlot\",\"ComparisonPlot\".", sep="")))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop(paste("Input for type=", type, ". However,'type' should be one of \"Heatmap\", \"VolcanoPlot\",\"ComparisonPlot\".", sep=""))
  	}
  
  	## check logBase.pvalue is 2,10 or not
  	if (logBase.pvalue != 2 & logBase.pvalue != 10) {
    	processout <- rbind(processout, c("ERROR : (-) Logarithm transformation for adjusted p-values : log2 or log10 only - stop"))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop("Only -log2 or -log10 for logarithm transformation for adjusted p-values are posssible.\n")
  	}
  
  	## choose comparison to draw plots
  
  	if (which.Comparison != "all") {
    	## check which.comparison is name of comparison
    	if (is.character(which.Comparison)) {
      
      		temp.name <- which.Comparison
      
      		## message if name of comparison is wrong.
      		if (length(setdiff(temp.name, unique(data$Label))) > 0) {
        
        		processout <- rbind(processout, paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse=", "), sep=" "))
      			write.table(processout, file=finalfile, row.names=FALSE)
      
      			stop(paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse=", "), sep=" "))
      		}
    	}
     
    	## check which.comparison is order number of comparison
    	if (is.numeric(which.Comparison)) {
      
    		temp.name <- levels(data$Label)[which.Comparison]
      
      		## message if name of comparison is wrong.
      		if (length(levels(data$Label))<max(which.Comparison)) {
        		stop(paste("Please check your selection of comparisons. There are ", length(levels(data$Label)), " comparisons in this result.", sep=" "))
        	}
    	}  
    
   	 	## use only assigned proteins
    	data <- data[which(data$Label %in% temp.name), ]
    
   	 	data$Protein <- factor(data$Protein)
    	data$Label <- factor(data$Label)
  	} else {
    
    	data$Protein <- factor(data$Protein)
    	data$Label <- factor(data$Label)
  	}
  
  
  	#######################
  	## Heatmap
  	#######################
  
  	if (type == "HEATMAP") {
    
    	## check whether there is only one protein or not,
    	if (length(unique(data$Protein)) <= 1) {
      		stop("At least two proteins are needed for heatmaps.")
    	}
    
    	## check whether there is only one comparison or not,
    	if (length(unique(data$Label)) <= 1) {
      		stop("At least two comparisons are needed for heatmaps.")
    	}
    
    	if (logBase.pvalue == 2) {
      		y.limUp  <- 30
    	}
    
    	if (logBase.pvalue == 10) {
      		y.limUp  <- 10
    	}
    
    	if (is.numeric(ylimUp)) y.limUp <- ylimUp 
    
   		## when NA, change it
    	#data$adj.pvalue[is.na(data$adj.pvalue)] <- 1 ## missing will be grey
    
    	if (logBase.pvalue == 2) {
      		data$adj.pvalue[data$adj.pvalue<2^(-y.limUp)] <- 2^(-y.limUp)
    	}
    
    	if (logBase.pvalue == 10) {
      		data$adj.pvalue[data$adj.pvalue<10^(-y.limUp)] <- 10^(-y.limUp)
    	}
    
    
    	## if FCcutoff is assigned, make p-value insignificant.
   		if (is.numeric(FCcutoff)) {
      		if (colnames(data)[3] == "log2FC") {
        		data$adj.pvalue[data[, 3] < log2(FCcutoff) & data[, 3] > (-log2(FCcutoff)) ] <- 1
      		}
      		if (colnames(data)[3] == "log10FC") {
        		data$adj.pvalue[data[, 3] < log10(FCcutoff) & data[, 3] > (-log10(FCcutoff))] <- 1
      		}
    	}
    
    	final <- NULL
    
        ## based on p-value
    	for (i in 1:nlevels(data$Label)) {
      
      		sub <- data[data$Label == levels(data$Label)[i], ]
      
      		if (logBase.pvalue==2) {
      		  temp <-  -log2(sub$adj.pvalue)*sign(sub[,3])
      		}
      
      		if (logBase.pvalue==10) {
        		temp <-  -log10(sub$adj.pvalue)*sign(sub[,3])
      		}
      
      		final <- data.frame(cbind(final,temp))
    	}
    
    	obj <- final
    	data$Protein <- factor(data$Protein)
    	rownames(obj) <- levels(data$Protein)
    	colnames(obj) <- levels(data$Label)
    
    	## remove if whole rows or columns are NA
    	obj <- obj[rowSums(!is.na(obj)) != 0, colSums(!is.na(obj)) != 0]
    
    	## clustering for order
    	tempobj <- obj
    	tempobj[is.na(tempobj)] <- 50
    
    	if (toupper(clustering) == 'PROTEIN') {
      		obj <- obj[hclust(dist(tempobj), method="ward.D")$order, ]
   	 	}
    	if (toupper(clustering) == 'COMPARISON') {
      		obj <- obj[, hclust(dist(t(tempobj)), method="ward.D")$order]
    	}
   	    if (toupper(clustering) == 'BOTH') {
      		obj <- obj[hclust(dist(tempobj), method="ward.D")$order, hclust(dist(t(tempobj)), method="ward.D")$order]
    	}
    	if (toupper(clustering) == 'NONE') {
      		obj <- obj
    	}
    
    	rm(tempobj)
    
    	## change the order
    	#obj$id <- seq(1:nrow(obj))
    	#obj <- obj[order(obj$id,decreasing=TRUE),]
    	#obj <- subset(obj, select=-c(id))
    
    	## color scale
    	blue.red.18  <-  maPalette(low = "blue", high = "red", mid = "black", k = 12)
    	my.colors  <- blue.red.18
    	#my.colors[my.colors=="#FFFFFF"] <- "gold"
    	my.colors <- c(my.colors,"grey") ## for NA
    
    	## color scale is fixed with log 10 based. then change break for log 2 based
    	up <- 10 
    		
    	temp <- 10^(-sort(ceiling(seq(2, up, length=10)[c(1, 2, 3, 5, 10)]), decreasing = TRUE))
    	breaks <- c(temp,sig)
    	
    	if (logBase.pvalue == 10) {
    		
    		neg.breaks  <-  log(breaks, 10)
    		my.breaks   <-  c(neg.breaks, 0, -neg.breaks[6:1], 101)	
    		
		} else if(logBase.pvalue == 2) {
			
			neg.breaks  <-  log(breaks, 2)
    		my.breaks   <-  c(neg.breaks, 0, -neg.breaks[6:1], 101)	
    		
		}
   
    	## draw color key
    	blocks <- c(-breaks, 1, breaks[6:1])
    	x.at <- seq(-0.05, 1.05, length.out=13)
    
    	## maximum number of proteins per heatmap
    	namepro <- rownames(obj)
    	totalpro <- length(namepro)
    	numheatmap <- totalpro %/% numProtein +1
    
    
    	## If there are the file with the same name, add next numbering at the end of file name
   	 	if (address != FALSE) {
      		allfiles <- list.files()
      
      		num <- 0
      		filenaming <- paste(address, "Heatmap", sep="")
      		finalfile <- paste(address, "Heatmap.pdf", sep="")
      
      		while(is.element(finalfile, allfiles)) {
        		num <- num + 1
        		finalfile <- paste(paste(filenaming, num, sep="-"), ".pdf", sep="")
      		}	
      
      		pdf(finalfile, width=width, height=height)
    	}
    
    	par(mar=c(3,3,3,3), mfrow=c(3,1),oma=c(3,0,3,0))
    	plot.new()
   		image(z = matrix(seq(1:(length(my.colors) - 1)), ncol = 1), col = my.colors[-length(my.colors)], xaxt = "n", yaxt = "n")
   	 	mtext("Color Key", side=3,line=1, cex=3)
    	mtext("(sign) Adjusted p-value", side=1, line=3,at=0.5, cex=1.7)
    	mtext(blocks, side=1, line=1, at=x.at, cex=1)
    
    	## draw heatmap
    
    	## loop for numProtein
    	for(j in 1:numheatmap) {
      
      		## get the number proteins needed
      		if (j != numheatmap) {
        		tempobj <- obj[((j-1) * numProtein + 1):(j * numProtein), ]
      		} else {
        		tempobj <- obj[((j-1) * numProtein + 1):nrow(obj), ]
      		}
      
      		par(oma=c(3,0,0,4))
     		heatmap.2(as.matrix(tempobj),
                col=my.colors,
                Rowv=FALSE,
                Colv=FALSE,
                dendrogram="none",
                breaks=my.breaks,
                trace="none",
                na.color="grey", ## missing value will be grey
                cexCol=(x.axis.size/10),cexRow=(y.axis.size/10), # assign text.size as option
                key=FALSE,
                lhei=c(0.1,0.9),
                lwid=c(0.1,0.9)
      		)   
      
    	} 
    	## end loop for heatmap
    
    	if (address!=FALSE) dev.off()
  	}
  
  
 	#######################
  	## VolcanoPlot
  	#######################
  	if (type == "VOLCANOPLOT") {
    
    	## If there are the file with the same name, add next numbering at the end of file name		
    	if (address != FALSE) {
      		allfiles <- list.files()
      
      		num <- 0
      		filenaming <- paste(address, "VolcanoPlot", sep="")
      		finalfile <- paste(address, "VolcanoPlot.pdf", sep="")
      
      		while(is.element(finalfile, allfiles)) {
        		num <- num + 1
        		finalfile <- paste(paste(filenaming, num, sep="-"), ".pdf", sep="")
      		}	
      
      		pdf(finalfile, width=width, height=height)
    	}
    
    	if (logBase.pvalue == 2) {
      		y.limUp  <- 30
    	}
    
    	if (logBase.pvalue == 10) {
     	 	y.limUp  <- 10
    	}
    
    	if (is.numeric(ylimUp)) y.limUp <- ylimUp 
    
    	## remove the result, NA
    	data <- data[!is.na(data$adj.pvalue),]
    
    	## group for coloring dots
    	if (!FCcutoff) {  
      		data[data$adj.pvalue >= sig, "colgroup"] <- "black"
      		data[data$adj.pvalue < sig & data[, 3] > 0, "colgroup"] <- "red"
      		data[data$adj.pvalue < sig & data[, 3] < 0, "colgroup"] <- "blue" 
    	}
    
    	if (is.numeric(FCcutoff)) {
      		data$colgroup <- "black"
      
      		if (colnames(data)[3] == "log2FC") {
        		data[data$adj.pvalue < sig & data[, 3] > log2(FCcutoff), "colgroup"] <- "red"
        		data[data$adj.pvalue < sig & data[, 3] < (-log2(FCcutoff)), "colgroup"] <- "blue"
      		}
      
      		if (colnames(data)[3] == "log10FC") {
        		data[data$adj.pvalue < sig & data[, 3] > log10(FCcutoff), "colgroup"] <- "red"
        		data[data$adj.pvalue < sig & data[, 3] < (-log10(FCcutoff)), "colgroup"] <- "blue"
      		}
    	}
    
     	data$colgroup <- factor(data$colgroup, levels=c("black", "blue", "red"))
    
    	## for multiple volcano plots, 
    	for(i in 1:nlevels(data$Label)) {
      
      		sub <- data[data$Label == levels(data$Label)[i], ]
      
      		if (logBase.pvalue == 2) {
        		sub$adj.pvalue[sub$adj.pvalue < 2^(-y.limUp)] <- 2^(-y.limUp)
      		}
      
      		if (logBase.pvalue == 10) {
        		sub$adj.pvalue[sub$adj.pvalue < 10^(-y.limUp)] <- 10^(-y.limUp)
      		}
      
      		sub <- as.data.frame(sub)
      
      		## ylimUp
      		if (logBase.pvalue == 2) {
        		y.limup <- ceiling(max(-log2(sub[!is.na(sub$adj.pvalue), "adj.pvalue"])))
        		if (y.limup < (-log2(sig))) {
        			y.limup <- (-log2(sig) + 1) ## for too small y.lim
        		}
      		}
      
     	 	if (logBase.pvalue == 10) {
        		y.limup <- ceiling(max(-log10(sub[!is.na(sub$adj.pvalue), "adj.pvalue"])))
        		if (y.limup < (-log10(sig))) {
        			y.limup <- (-log10(sig) + 1) ## for too small y.lim
        		}
      		}
       
      		## ylimDown
      		y.limdown <- 0 ## default is zero
      		if (is.numeric(ylimDown)) {
      			y.limdown <- ylimDown
      		}
      
      		## x.lim
      		x.lim <- ceiling(max(abs(sub[!is.na(sub[, 3]) & abs(sub[, 3]) != Inf , 3]))) ## log2FC or log10FC
      		if (x.lim < 3) {
      			x.lim <- 3
      		}
      		if (is.numeric(xlimUp)) {
      			x.lim <- xlimUp
      		}
      
      		## for assigning x in ggplot2
      		subtemp <- sub
      		colnames(subtemp)[3] <- "logFC"
      
      		if (logBase.pvalue == 2) {
        		subtemp$log2adjp <- (-log2(subtemp$adj.pvalue))
      		}
      
      		if (logBase.pvalue == 10) {
        		subtemp$log10adjp <- (-log10(subtemp$adj.pvalue))
        	}
        	
        	## for x limit for inf or -inf
        	subtemp$newlogFC <- subtemp$logFC
        	subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing" & subtemp$logFC == Inf, "newlogFC"] <- (x.lim - 0.2)
        	subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing" & subtemp$logFC == (-Inf), "newlogFC"] <- (x.lim - 0.2) *(-1)
        	
        	## add (*) in Protein name for Inf or -Inf
        	subtemp$Protein <- as.character(subtemp$Protein)
        	subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing", "Protein"] <- paste("*", subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing", "Protein"], sep="")

      
      		## Plotting
      		if (logBase.pvalue == 2) {
        		ptemp <- ggplot(aes_string(x='logFC', y='log2adjp', color='colgroup', label='Protein'), data=subtemp)+
        				geom_point(size=dot.size)+
        				scale_colour_manual(values=c("gray65", "blue", "red"), limits=c("black", "blue", "red"), breaks=c("black", "blue", "red"), labels=c("No regulation", "Down-regulated", "Up-regulated"))+
        				scale_y_continuous('-Log2 (adjusted p-value)', limits=c(y.limdown, y.limup))+
        				labs(title=unique(sub$Label))
      		}
      
      		if (logBase.pvalue == 10) {
        		ptemp <- ggplot(aes_string(x='logFC', y='log10adjp', color='colgroup', label='Protein'), data=subtemp)+
        				geom_point(size=dot.size)+
        				scale_colour_manual(values=c("gray65", "blue", "red"), limits=c("black", "blue", "red"), breaks=c("black", "blue", "red"), labels=c("No regulation", "Down-regulated", "Up-regulated"))+
        				scale_y_continuous('-Log10 (adjusted p-value)', limits=c(y.limdown, y.limup))+
        				labs(title=unique(sub$Label))
        	}
      
      
     		## x-axis labeling
      		if (colnames(sub)[3] == "log2FC") {
      			ptemp <- ptemp+scale_x_continuous('Log2 fold change', limits=c(-x.lim, x.lim))
      		}
      		if (colnames(sub)[3] == "log10FC") {
      			ptemp <- ptemp+scale_x_continuous('Log10 fold change', limits=c(-x.lim, x.lim))
      		}
      
     		## add protein name
      		if (ProteinName) {
      			if(length(unique(subtemp$colgroup)) == 1 & any(unique(subtemp$colgroup) == 'black')){
      				message(paste("The volcano plot for ", unique(subtemp$Label), " does not show the protein names because none of them is significant.", sep=""))
      				
      			} else {
          			
          			ptemp <- ptemp + geom_text_repel(data=subtemp[subtemp$colgroup != "black", ], aes(label=Protein), size=text.size, col='black')
          		}
      		} 
      
      		## For legend of linetype for cutoffs
      		## first assign line type
      		ltypes <- c("type1"="twodash", "type2"="dotted")
      
      		## cutoff lines, FDR only
     		if (!FCcutoff) { 
        		if (logBase.pvalue == 2) {
          			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=20), log2adjp=(-log2(sig)), line='twodash')
          
          			pfinal <- ptemp + geom_line(data=sigcut, aes_string(x='logFC', y='log2adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            			    scale_linetype_manual(values=c('twodash'=6), labels=c(paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
            			    		
       	 		}
        
        		if (logBase.pvalue == 10) {
          			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=20), log10adjp=(-log10(sig)), line='twodash')
          
         	 		pfinal <- ptemp + geom_line(data=sigcut, aes_string(x='logFC', y='log10adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            			    scale_linetype_manual(values=c('twodash'=6), labels=c(paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
				}				
      		}
      
      		## cutoff lines, FDR and Fold change cutoff
      		if (is.numeric(FCcutoff)) {
        		if (colnames(sub)[3] == "log2FC") {
          			if (logBase.pvalue == 2) {
            
           				## three different lines
            			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=10), log2adjp=(-log2(sig)), line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', logFC=log2(FCcutoff), log2adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', logFC=(-log2(FCcutoff)), log2adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log2adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log2adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log2adjp', linetype='line'), colour="darkgrey", size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}
          
          			if (logBase.pvalue == 10) {
            
            			## three different lines
            			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=10), log10adjp=(-log10(sig)), line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', logFC=log2(FCcutoff), log10adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', logFC=(-log2(FCcutoff)), log10adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log10adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log10adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log10adjp', linetype='line'), colour="darkgrey", size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}        
        		}
        
        		if (colnames(sub)[3] == "log10FC") {
          			if (logBase.pvalue == 2) {
            
            			## three different lines
            			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=10), log2adjp=(-log2(sig)), line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', logFC=log10(FCcutoff), log2adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', logFC=(-log10(FCcutoff)), log2adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log2adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log2adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log2adjp', linetype='line'), colour="darkgrey", size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}
          
          			if (logBase.pvalue == 10) {
            
            			## three different lines
            			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=10), log10adjp=(-log10(sig)), line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', logFC=log10(FCcutoff), log10adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', logFC=(-log10(FCcutoff)), log10adjp=seq(y.limdown, y.limup, length.out=10), line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log10adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log10adjp', linetype='line'), colour="darkgrey", size=0.6, show.legend=T)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log10adjp', linetype='line'), colour="darkgrey", size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}    
        		}
     		 }
      
      		pfinal <- pfinal+theme(
        		panel.background = element_rect(fill='white', colour="black"),
        		panel.grid.minor = element_blank(),
        		axis.text.x = element_text(size=x.axis.size, colour="black"),
        		axis.text.y = element_text(size=y.axis.size, colour="black"),
        		axis.ticks = element_line(colour="black"),
        		axis.title.x = element_text(size=x.axis.size+5, vjust=-0.4),
        		axis.title.y = element_text(size=y.axis.size+5, vjust=0.3),
        		title = element_text(size=x.axis.size+8, vjust=1.5),
        		legend.position="bottom",
        		legend.key = element_rect(fill='white', colour='white'),
        		legend.text = element_text(size=legend.size),
        		legend.title = element_blank()
        		)
      
      		print(pfinal)
    	} ## end-loop
    
    	if (address!=FALSE) dev.off()
	}	
  
  	#######################
  	## Comparison Plot
  	#######################
  	if (type == "COMPARISONPLOT") {
    
    	datatemp <- data[!is.na(data$adj.pvalue), ]
    	datatemp$Protein <- factor(datatemp$Protein)
    
    	## If there are the file with the same name, add next numbering at the end of file name		
    	if (address!=FALSE) {
     		allfiles <- list.files()
      
      		num <- 0
      		filenaming <- paste(address, "ComparisonPlot", sep="")
      		finalfile <- paste(address, "ComparisonPlot.pdf", sep="")
      
      		while(is.element(finalfile, allfiles)) {
        		num <- num+1
        		finalfile <- paste(paste(filenaming, num, sep="-"), ".pdf", sep="")
      		}	
      
      		pdf(finalfile, width=width, height=height)
   		}
    
   	 	for (i in 1:nlevels(datatemp$Protein)) {
      
      		sub <- datatemp[datatemp$Protein==levels(datatemp$Protein)[i], ] 		
      		#sub$ciw <- qt(1-sig/2,sub$DF)*sub$SE
      		## adjust for multiple comparison within protein
      		sub$ciw <- qt(1-sig/(2*nrow(sub)), sub$DF)*sub$SE
      
      		sub <- as.data.frame(sub)
      		
      		## for assigning x in ggplot2
      		colnames(sub)[3] <- "logFC"
      
      		## ylimUp
      		y.limup <- ceiling(max(sub$logFC + sub$ciw))
      		if (is.numeric(ylimUp)) {
      			y.limup <- ylimUp 
      		}
      
      		## ylimDown
      		y.limdown <- floor(min(sub$logFC - sub$ciw))
      		if (is.numeric(ylimDown)) {
      			y.limdown <- ylimDown 
      		}
      
      		ptemp <- ggplot(aes_string(x='Label', y='logFC'), data=sub)+
      			geom_errorbar(aes(ymax = logFC + ciw, ymin=logFC - ciw), data=sub, width=0.1, colour="red")+
      			geom_point(size=dot.size, colour="darkred")+
      			scale_x_discrete('Comparison')+
      			geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+
      			labs(title=levels(datatemp$Protein)[i])+
      			theme(
        			panel.background=element_rect(fill='white', colour="black"),
        			panel.grid.major.y = element_line(colour="grey95"),
        			panel.grid.minor.y = element_blank(),
        			axis.text.x = element_text(size=x.axis.size, colour="black"),
        			axis.text.y = element_text(size=y.axis.size, colour="black"),
        			axis.ticks = element_line(colour="black"),
        			axis.title.x = element_text(size=x.axis.size+5, vjust=-0.4),
        			axis.title.y = element_text(size=y.axis.size+5, vjust=0.3),
        			title = element_text(size=x.axis.size+8, vjust=1.5)
        			)
      
      		if (colnames(data)[3] == "log2FC") {
      			ptemp <- ptemp+scale_y_continuous("Log2-Fold Change", limits=c(y.limdown, y.limup))
      		}
     	 	if (colnames(data)[3] == "log10FC") {
     	 		ptemp <- ptemp+scale_y_continuous("Log10-Fold Change", limits=c(y.limdown, y.limup))
     	 	}
      
      		print(ptemp)
		
		} ## end-loop
    
   	 	if (address!=FALSE) {
    		dev.off()
    	}
  	} ## end Comparison plot
}
