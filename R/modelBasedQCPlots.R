#' @export
modelBasedQCPlots <- function(data,
                              type,
                              axis.size=10,
                              dot.size=3,
                              text.size=7,
                              legend.size=7,
                              width=10, 
                              height=10,
                              address="") {
  
  if (length(setdiff(toupper(type),c("QQPLOTS","RESIDUALPLOTS")))!=0) {
    stop(paste("Input for type=",type,". However,'type' should be one of \"QQPlots\", \"ResidualPlots\" .",sep=""))
  }
  
  datafit <- data$fittedmodel
  proid <- levels(data$ComparisonResult$Protein)
  
  #############################################
  ### normality, QQ plot
  #############################################
  if (toupper(type)=="QQPLOTS") {
    
    #### save the plots as pdf or not
    # If there are the file with the same name, add next numbering at the end of file name	
    if (address!=FALSE) {
      allfiles <- list.files()
        
      num <- 0
      filenaming <- paste(address,"QQPlot",sep="")
      finalfile <- paste(address,"QQPlot.pdf",sep="")
        
      while(is.element(finalfile,allfiles)) {
        num <- num+1
        inalfile <- paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
        
      pdf(finalfile, width=width, height=height)
    }
      
    for (i in 1:length(datafit)) {	
        
      sub <- datafit[[i]]
      
      if ( class(sub)=="lm" ) {  ### lm model
        sub.residuals <- sub$residuals
      } else {   ### lmer model
        sub.residuals <- resid(sub)
      }
      sub.residuals.table <- data.frame(residual=sub.residuals)
        
      ## get slope and intercept for qline
      y <- quantile(sub.residuals.table$residual[!is.na(sub.residuals.table$residual)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y)/diff(x)
      int <- y[1L]-slope*x[1L]
        
      ptemp <- ggplot(sub.residuals.table, aes(sample=residual))+
        geom_point(stat="qq",alpha=0.8, shape=20, size=dot.size)+
        scale_shape(solid=FALSE)+ 
        geom_abline(slope = slope, intercept = int,colour="red")+
        scale_y_continuous('Sample Quantiles')+
        scale_x_continuous('Theoretical Quantiles')+
        labs(title=paste("Normal Q-Q Plot (", proid[i], ")"))+
        theme(
          panel.background=element_rect(fill='white', colour="black"),
          panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor =element_blank(),
          axis.text.x=element_text(size=axis.size,colour="black"),
          axis.text.y=element_text(size=axis.size,colour="black"),
          axis.ticks=element_line(colour="black"),
          axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
          axis.title.y=element_text(size=axis.size+5,vjust=0.3),
          title=element_text(size=axis.size+8,vjust=1.5),
          legend.position="none"
        )
        
      print(ptemp)
        
      message(paste("Drew the QQ plot for ", proid[i], "(",i," of ",length(proid),")"))
        
    } ## end loop
      
    if (address!=FALSE) dev.off()
    
  } ### end QQplots
  
  
  #############################################
  #### Residual plot
  #############################################
  if (toupper(type)=="RESIDUALPLOTS") {
    
    #### save the plots as pdf or not
    # If there are the file with the same name, add next numbering at the end of file name	
    if (address!=FALSE) {
      allfiles <- list.files()
      
      num <- 0
      filenaming <- paste(address,"ResidualPlot",sep="")
      finalfile <- paste(address,"ResidualPlot.pdf",sep="")
      
      while(is.element(finalfile,allfiles)) {
        num <- num+1
        finalfile <- paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    for (i in 1:length(datafit)) {	
      
      sub <- datafit[[i]]
      
      if ( class(sub)=="lm" ) {  ### lm model
        sub.residuals <- sub$residuals
        sub.fitted <- sub$fitted.values
      } else {   ### lmer model
        sub.residuals <- resid(sub)
        sub.fitted <- fitted(sub)
      }
      
      sub.residuals.table <- data.frame(residual=sub.residuals, fitted=sub.fitted)
      
      ptemp <- ggplot(aes_string(x='fitted', y='residual'), data=sub.residuals.table)+
        geom_point(size=dot.size, alpha=0.5)+
        geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+
        scale_y_continuous('Residuals')+
        scale_x_continuous('Predicted Abundance')+
        labs(title=proid[i])+
        theme(
          panel.background=element_rect(fill='white', colour="black"),
          panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank(),
          axis.text.x=element_text(size=axis.size,colour="black"),
          axis.text.y=element_text(size=axis.size,colour="black"),
          axis.ticks=element_line(colour="black"),
          axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
          axis.title.y=element_text(size=axis.size+5,vjust=0.3),
          title=element_text(size=axis.size+8,vjust=1.5),
          legend.position="none")
    
      print(ptemp)
      
      message(paste("Drew the residual plot for ", proid[i], "(",i," of ",length(proid),")"))
      
    }
    
    if (address!=FALSE) dev.off()	
  } ## end residualplots
  
}


