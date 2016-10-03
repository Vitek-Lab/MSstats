
#The goal of this function is to calculate the LOB/LOD of the data provided in the data frame.
#The function returns a new data frame containing the value of the LOB/LOD

#' @export
plot_quantlim <- function(spikeindata, quantlim_out, alpha, dir_output, xlim_plot){
  
  

  
  if(is.null(quantlim_out)){
    print("Assay fit was incorrectly calculated by linear_quantlim or nonlinear_quantlim and cannot be plotted")
    return(NULL)
  }
  
  
  #percentile of the prediction interval considered
  if(missing(alpha)){
    alpha = 5/100;
  }
  
  if( alpha >= 1 | alpha  <= 0){
    print("incorrect specified value for alpha,  0 < alpha < 1")
    return(NULL)
  }
  
  expdata = spikeindata
  

  
  
  datain = quantlim_out
  #Define some colors here for the plots:
  black1 <- '#000000';  orange1 <- "#E69F00"; blue1 <- "#56B4E9"; green1 <- "#009E73"; yellow1 <- "#F0E442"; blue2 <- "#0072B2"; red1 <- "#D55E00"; pink1 <-  "#CC79A7";
  cbbPalette <- c(black1, orange1, blue1, green1, yellow1, blue2 ,red1 , pink1)
  
  #Rename variables for the function:
  names(datain)[names(datain) == 'CONCENTRATION'] <- 'C'
  names(expdata)[names(expdata) == 'CONCENTRATION'] <- 'C'
  names(expdata)[names(expdata) == 'INTENSITY'] <- 'I'
  names(datain)[names(datain) == 'MEAN'] <- 'mean'
  names(datain)[names(datain) == 'LOW'] <- 'low'
  names(datain)[names(datain) == 'UP'] <- 'up'
  
  #Remove NA and infinite numbers from spike in data:
  expdata <- expdata[!is.na(expdata$I) & !is.na(expdata$C),]  
  expdata <- expdata[!is.infinite(expdata$I) & !is.infinite(expdata$C),]  
  
  
  #Extract actual data points for plotting:
  Cdata = expdata$C
  Idata = expdata$I 
  
  tmp_blank <- expdata[expdata$C == 0,]
  n_blank = length(unique(tmp_blank$I))
  noise = mean(tmp_blank$I)
  var_noise = var(tmp_blank$I)
  
  
  fac = qt(1-alpha,n_blank - 1)*sqrt(1+1/n_blank)
  
  #upper bound of noise prediction interval
  up_noise = noise   +  fac* sqrt(var_noise)
  
  rel_size = 2.5; rel_size_2 = 1.8; lw = 1; pw = 2.5;
  
  xaxis_orig_2 = datain$C
  tmp_all = datain
  LOQ_pred = datain$LOD[1]
  LOD_pred = datain$LOB[1]
  
  lower_Q_pred = datain$low
  upper_Q_pred = datain$up
  mean_bilinear <- datain$mean
  
  if(LOD_pred >= 0){
    y_LOD_pred = up_noise
  }
  if(LOQ_pred >= 0){
    y_LOQ_pred = up_noise
  }  
  
  
  filename = paste(dir_output,'/',datain$NAME[1],'_',datain$METHOD[1],'_overall.pdf',sep='')
  pdf(filename)
  if(LOQ_pred > xaxis_orig_2[3]){
    C_max = xaxis_orig_2[min(which(abs(LOQ_pred - xaxis_orig_2) == min(abs(LOQ_pred - xaxis_orig_2))) +1, length(xaxis_orig_2))]      }else{
      C_max = xaxis_orig_2[which(abs(mean(xaxis_orig_2) - xaxis_orig_2) == min(abs(mean(xaxis_orig_2) - xaxis_orig_2)))]
    }
  
  low_p <- paste(alpha*100,"%",sep = "")
  high_p <- paste(100-alpha*100,"%", sep = "")
  upp_noise <- paste(high_p," upper bound of noise")
  low_pred <-  paste(low_p , 'percentile of predictions') 
  
  p1 <- ggplot() + .theme_complete_bw()
  p1 <- p1 + geom_point(data=data.frame(Cdata,Idata) , aes(x=Cdata,y=Idata),size =pw*1.5)
  p1 <- p1 +  geom_line(data=data.frame(x=xaxis_orig_2,y=mean_bilinear,col='mean prediction', lt = 'mean'),aes(x=x,y=y, color=col), size = lw) #, color = 'black' #as.factor(c('mean bootstrap'))
  p1 <- p1 +   geom_ribbon(data=data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred), aes(x=x, ymin=ymin, ymax=ymax), fill = red1, alpha = 0.3)
  p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =lower_Q_pred, col = low_pred, lt = 'Int') , aes(x=x,y=ymin, color = col), size = lw)
  p1 <- p1 +  geom_line(data=data.frame(x=xaxis_orig_2, ymax = rep(up_noise,length(xaxis_orig_2)),
                                                                                                 col = "95% upper bound of noise"),   aes(x=x, y=ymax, color = col),  size  = lw)
  p1 <- p1 + scale_alpha_continuous(guide = 'none') 
  p1 <- p1 + xlab('Spiked Concentration') + ylab('Estimated Concentration')
  p1 <- p1 + theme(axis.text.x = element_text(size = rel(rel_size))) + theme(axis.text.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.x = element_text(size = rel(rel_size))) + theme(axis.title.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.y=element_text(vjust=0.7))
  
  p1 <- p1 + scale_color_manual(values = c(orange1, blue1, red1),
                                labels = c(low_pred,
                                           upp_noise,
                                           "mean prediction"))
  
  #p1 <- p1 + scale_colour_manual(values=c('mean prediction'= red1,'upper 95% prediction' = orange1, '5% percentile of predictions' = orange1, '95% upper bound of noise' = blue1))  + theme(legend.title = element_blank()) 
  #colour_scales <- setNames(c('mean prediction','upper 95% prediction','5% percentile of predictions'),c('dasssta','messssan',"rrrr"))  
  #p1 <- p1 + scale_colour_manual(values = colour_scales)
    
  p1 <- p1 + theme(legend.title = element_blank()) + theme(legend.position = c(0.05, 0.5), legend.justification = c(0, 0), legend.text=element_text(size=rel(rel_size_2)))
  LOD_y = mean_bilinear[which(abs(xaxis_orig_2 - LOD_pred) == min(abs(xaxis_orig_2 - LOD_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOD_pred, y= y_LOD_pred, shape='LOD'),aes(x=x,y=y,shape = shape, guide=FALSE),colour="purple",size=5)
  LOQ_y = lower_Q_pred[which(abs(up_noise - lower_Q_pred) == min(abs(up_noise - lower_Q_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOQ_pred, y= y_LOQ_pred, shape='LOQ'),aes(x=x,y=y,shape = shape, guide=FALSE),colour=orange1,size=5)
  LOD_string <- paste('LOB=',round(LOD_pred, digits=1),sep='');   LOQ_string <- paste('LOD=',round(LOQ_pred, digits=1),sep=''); 
  p1 <- p1 + guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),     shape = guide_legend(order = 2)) + guides(shape=FALSE)
  p1 <- p1 + ggtitle(paste(datain$NAME,'\n', LOD_string,', ',LOQ_string,sep="")) + theme(plot.title = element_text(size = 20))
  print(p1);
  dev.off()
  

  #produce a second plot showing a zoomed view:
  if(missing(xlim_plot)){ #missing argument for the x limit in the function, pick a x limit that is close to the LOD/LOQ:
    
    if(LOQ_pred > 0){
      xlim = LOQ_pred*3.0      
    } else{
      xlim = unique(Cdata)[4]
    }
  
    
  } else{
    
    xlim = xlim_plot
    
  }
  
  

  filename = paste(dir_output,'/',datain$NAME[1],'_',datain$METHOD[1],'_zoom.pdf',sep='')
  pdf(filename)

  

  Idata =  subset(Idata, Cdata < xlim)
  Cdata =  subset(Cdata, Cdata < xlim)
    

  
  lower_Q_pred = subset(lower_Q_pred, xaxis_orig_2 < xlim)
  upper_Q_pred = subset(upper_Q_pred, xaxis_orig_2 < xlim)  
  mean_bilinear = subset(mean_bilinear, xaxis_orig_2 < xlim)
  xaxis_orig_2 = subset(xaxis_orig_2, xaxis_orig_2 < xlim)
  
  p1 <- ggplot() + .theme_complete_bw()
  p1 <- p1 + geom_point(data=data.frame(Cdata,Idata) , aes(x=Cdata,y=Idata),size =pw*1.5)
  p1 <- p1 +  geom_line(data=data.frame(x=xaxis_orig_2,y=mean_bilinear,col='mean prediction', lt = 'mean'),aes(x=x,y=y, color=col), size = lw) #, color = 'black' #as.factor(c('mean bootstrap'))
  p1 <- p1 +   geom_ribbon(data=data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred), aes(x=x, ymin=ymin, ymax=ymax), fill = red1, alpha = 0.3)
  p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =lower_Q_pred, col = 'lower 95% prediction', lt = 'Int') , aes(x=x,y=ymin, color = col), size = lw)
  p1 <- p1 +  geom_line(data=data.frame(x=xaxis_orig_2, ymax = rep(up_noise,length(xaxis_orig_2)),
                                        col = "95% upper bound of noise"),   aes(x=x, y=ymax, color = col),  size  = lw)
  p1 <- p1 + scale_alpha_continuous(guide = 'none') 
  p1 <- p1 + xlab('Spiked Concentration') + ylab('Estimated Concentration')
  p1 <- p1 + theme(axis.text.x = element_text(size = rel(rel_size))) + theme(axis.text.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.x = element_text(size = rel(rel_size))) + theme(axis.title.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.y=element_text(vjust=0.7))
  
  p1 <- p1 + scale_color_manual(values = c(blue1,orange1, red1),
                                labels = c(upp_noise,low_pred,
                                           "mean prediction"))
  
  
  #p1 <- p1 + scale_colour_manual(values=c('mean prediction'= red1,'upper 95% prediction' = orange1, 'lower 95% prediction' = orange1, '95% upper bound of noise' = blue1))  + theme(legend.title = element_blank()) 
  #p1 <- p1 + scale_linetype_manual(values=c('dashed','dashed','solid','solid'))
  p1 <- p1 + theme(legend.title = element_blank()) + theme(legend.position = c(0.05, 0.5), legend.justification = c(0, 0), legend.text=element_text(size=rel(rel_size_2)))
  LOD_y = mean_bilinear[which(abs(xaxis_orig_2 - LOD_pred) == min(abs(xaxis_orig_2 - LOD_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOD_pred, y= y_LOD_pred, shape='LOD'),aes(x=x,y=y,shape = shape, guide=FALSE),colour="purple",size=5)
  LOQ_y = lower_Q_pred[which(abs(up_noise - lower_Q_pred) == min(abs(up_noise - lower_Q_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOQ_pred, y= y_LOQ_pred, shape='LOQ'),aes(x=x,y=y,shape = shape, guide=FALSE),colour=orange1,size=5)
  LOD_string <- paste('LOB=',round(LOD_pred, digits=1),sep='');   LOQ_string <- paste('LOD=',round(LOQ_pred, digits=1),sep=''); 
  p1 <- p1 + guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),     shape = guide_legend(order = 2)) + guides(shape=FALSE)
  p1 <- p1 + ggtitle(paste(datain$NAME,'\n', LOD_string,', ',LOQ_string,sep="")) + theme(plot.title = element_text(size = 20))
  print(p1);
  dev.off()
  

}
