
#The goal of this function is to calculate the LoD/LoQ of the data provided in the data frame.
#The function returns a new data frame containing the value of the LoD/LoQ

plot_quantlim <- function(spikeindata,quantlim_out,dir_output,xlim_plot){
  
  
  expdata = spikeindata
  datain = quantlim_out
  #Define some colors here for the plots:
  black1 <- '#000000';  orange1 <- "#E69F00"; blue1 <- "#56B4E9"; green1 <- "#009E73"; yellow1 <- "#F0E442"; blue2 <- "#0072B2"; red1 <- "#D55E00"; pink1 <-  "#CC79A7";
  cbbPalette <- c(black1, orange1, blue1, green1, yellow1,blue2 ,red1 , pink1)
  
  #Rename variables for the function:
  names(datain)[names(datain) == 'CONCENTRATION'] <- 'C'
  names(expdata)[names(expdata) == 'CONCENTRATION'] <- 'C'
  names(expdata)[names(expdata) == 'INTENSITY'] <- 'I'
  names(datain)[names(datain) == 'MEAN'] <- 'mean'
  names(datain)[names(datain) == 'LOW'] <- 'low'
  names(datain)[names(datain) == 'UP'] <- 'up'
  
  #Extract actual data points for plotting:
  Cdata = expdata$C
  Idata = expdata$I 
  
  tmp_blank <- expdata[expdata$C == 0,]
  n_blank = length(unique(tmp_blank$I))
  noise = mean(tmp_blank$I)
  var_noise = var(tmp_blank$I)
  fac_low = qt(1-0.25/2,n_blank - 1)*sqrt(1+1/n_blank)  
  fac = qt(1-0.05/2,n_blank - 1)*sqrt(1+1/n_blank)
  fac_high = qt(1-0.01/2,n_blank - 1)*sqrt(1+1/n_blank)
  
  #2.575
  
  low_noise = noise  - fac* sqrt(var_noise)
  
  up_noise_low = noise + fac_low * sqrt(var_noise)
  up_noise = noise   +  fac* sqrt(var_noise)
  up_noise_high = noise + fac_high*sqrt(var_noise)
  
  rel_size = 2.5; rel_size_2 = 1.8; lw = 1; pw = 2.5;
  
  xaxis_orig_2 = datain$C
  tmp_all = datain
  LOQ_pred = datain$LOQ[1]
  LOD_pred = datain$LOD[1]
  
  lower_Q_pred = datain$low
  upper_Q_pred = datain$up
  mean_bilinear <- datain$mean
  
  if(LOD_pred > 0){
    y_LOD_pred = up_noise
  }
  if(LOQ_pred > 0){
    y_LOQ_pred = up_noise
  }  
  
  
  filename = paste(dir_output,'/',datain$NAME[1],'_',datain$METHOD[1],'_overall.pdf',sep='')
  pdf(filename)
  if(LOQ_pred > xaxis_orig_2[3]){
    C_max = xaxis_orig_2[min(which(abs(LOQ_pred - xaxis_orig_2) == min(abs(LOQ_pred - xaxis_orig_2))) +1, length(xaxis_orig_2))]      }else{
      C_max = xaxis_orig_2[which(abs(mean(xaxis_orig_2) - xaxis_orig_2) == min(abs(mean(xaxis_orig_2) - xaxis_orig_2)))]
    }
  
  
  p1 <- ggplot() + .theme_complete_bw()
  p1 <- p1 + geom_point(data=data.frame(Cdata,Idata) , aes(x=Cdata,y=Idata),size =pw*1.5)
  p1 <- p1 +  geom_line(data=data.frame(x=xaxis_orig_2,y=mean_bilinear,col='mean prediction', lt = 'mean'),aes(x=x,y=y, color=col), size = lw) #, color = 'black' #as.factor(c('mean bootstrap'))
  #p1 <- p1 +   scale_colour_manual(values=c("uuuu"="grey","mean bootstrap"=ggplotColours(n=2)[1])) + scale_alpha(guide = 'none')
  p1 <- p1 +   geom_ribbon(data=data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred), aes(x=x, ymin=ymin, ymax=ymax), fill = .ggplotColours(n=2)[1],alpha=0.1)
  p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred, col = '95% prediction', lt = 'Int') , aes(x=x,y=ymin, linetype = col), color = red1, size = lw)  #, color = 'red' #size=2.5,
  p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred, col = '95% prediction', lt = 'Int'), aes(x=x,y=ymax, linetype = col), color = red1, size  = lw) #, color = 'red' , size=2.5
  
  
  p1 <- p1 +  geom_ribbon(data=data.frame(x=xaxis_orig_2,y=rep(noise,length(xaxis_orig_2)), ymin = rep(low_noise,length(xaxis_orig_2)) , ymax = rep(up_noise,length(xaxis_orig_2)) ),
                          aes(x=x, ymin=ymin, ymax=ymax), fill = .ggplotColours(n=2)[2],alpha=0.3)
  
  #p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =rep(low_noise,length(xaxis_orig_2)) , ymax =upper_Q, col = '95% noise') %>% filter(x <= C_max), aes(x=x,y=ymin,color= col, linetype =col ), size = lw )
  #p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =lower_Q , ymax =rep(up_noise,length(xaxis_orig_2)), col = '95% noise') %>% filter(x <= C_max), aes(x=x,y=ymax, color= col, linetype =col ), size = lw )
  
  p1 <- p1 + scale_alpha_continuous(guide = 'none') 
  p1 <- p1 + xlab('Spiked Concentration') + ylab('Estimated Concentration')
  
  p1 <- p1 + theme(axis.text.x = element_text(size = rel(rel_size))) + theme(axis.text.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.x = element_text(size = rel(rel_size))) + theme(axis.title.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.y=element_text(vjust=0.7))
  
  p1 <- p1 + scale_colour_manual(values=c('Bootstrap samples'='grey','mean prediction'= red1,'95% prediction' = red1, '95% noise' = blue1))  + theme(legend.title = element_blank()) 
  p1 <- p1 + scale_linetype_manual(values=c('dashed','dashed','solid'))
  p1 <- p1 + theme(legend.title = element_blank()) + theme(legend.position = c(0.05, 0.5), legend.justification = c(0, 0), legend.text=element_text(size=rel(rel_size_2)))
  
  LOD_y = mean_bilinear[which(abs(xaxis_orig_2 - LOD_pred) == min(abs(xaxis_orig_2 - LOD_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOD_pred, y= y_LOD_pred, shape='LOD'),aes(x=x,y=y,shape = shape, guide=FALSE),colour="purple",size=5)
  
  LOQ_y = lower_Q_pred[which(abs(up_noise - lower_Q_pred) == min(abs(up_noise - lower_Q_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOQ_pred, y= y_LOQ_pred, shape='LOQ'),aes(x=x,y=y,shape = shape, guide=FALSE),colour=orange1,size=5)
  LOD_string <- paste('LOD=',round(LOD_pred, digits=1),sep='');   LOQ_string <- paste('LOQ=',round(LOQ_pred, digits=1),sep=''); 
 # p1 <- p1 +  scale_shape_manual( values=c(16,17), guide=FALSE)# ,labels=c(LOD_string, LOQ_string)) 
  
  p1 <- p1 + guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),     shape = guide_legend(order = 2)) + guides(shape=FALSE)
  #p1 <- p1 + xlim(0, LOQ_pred) + ylim(0,LOQ_y)

  #p1 <- p1 + annotate("text", x = 0.5*max(Cdata), y = 0.1*max(Idata), label = LOD_string, colour = "purple", size=5)
  #p1 <- p1 + annotate("text", x = 0.6*max(Cdata), y = 0.2*max(Idata), label = LOQ_string, colour = orange1, size=5)
  
  p1 <- p1 + ggtitle(paste(datain$NAME,': ', LOD_string,', ',LOQ_string,sep="")) + theme(plot.title = element_text(size = 20))
  
  print(p1);
  dev.off()
  

  #produce a second plot showing a zoomed view:
  if(missing(xlim_plot)){ #missing argument for the x limit in the function, pick a x limit that is close to the LOD/LOQ:
    
    if(LOQ_pred > 0){
      xlim = LOQ_pred*2      
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
  #p1 <- p1 +   scale_colour_manual(values=c("uuuu"="grey","mean bootstrap"=ggplotColours(n=2)[1])) + scale_alpha(guide = 'none')
  p1 <- p1 +   geom_ribbon(data=data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred), aes(x=x, ymin=ymin, ymax=ymax), fill = .ggplotColours(n=2)[1],alpha=0.1)

  p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred, col = '95% prediction', lt = 'Int') , aes(x=x,y=ymin, linetype = col), color = red1, size = lw)  #, color = 'red' #size=2.5,
  p1 <- p1 + geom_line(data = data.frame(x=xaxis_orig_2,ymin =lower_Q_pred , ymax =upper_Q_pred, col = '95% prediction', lt = 'Int'), aes(x=x,y=ymax, linetype = col), color = red1, size  = lw) #, color = 'red' , size=2.5
  
  
  p1 <- p1 +  geom_ribbon(data=data.frame(x=xaxis_orig_2,y=rep(noise,length(xaxis_orig_2)), ymin = rep(low_noise,length(xaxis_orig_2)) , ymax = rep(up_noise,length(xaxis_orig_2)) ),
                          aes(x=x, ymin=ymin, ymax=ymax), fill = .ggplotColours(n=2)[2],alpha=0.3)
  

  p1 <- p1 + scale_alpha_continuous(guide = 'none') 
  p1 <- p1 + xlab('Spiked Concentration') + ylab('Estimated Concentration')
  
  p1 <- p1 + theme(axis.text.x = element_text(size = rel(rel_size))) + theme(axis.text.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.x = element_text(size = rel(rel_size))) + theme(axis.title.y = element_text(size = rel(rel_size)))
  p1 <- p1 + theme(axis.title.y=element_text(vjust=0.7))
  
  p1 <- p1 + scale_colour_manual(values=c('Bootstrap samples'='grey','mean prediction'= red1,'95% prediction' = red1, '95% noise' = blue1))  + theme(legend.title = element_blank()) 
  p1 <- p1 + scale_linetype_manual(values=c('dashed','dashed','solid'))
  p1 <- p1 + theme(legend.title = element_blank()) + theme(legend.position = c(0.05, 0.5), legend.justification = c(0, 0), legend.text=element_text(size=rel(rel_size_2)))
  
  LOD_y = mean_bilinear[which(abs(xaxis_orig_2 - LOD_pred) == min(abs(xaxis_orig_2 - LOD_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOD_pred, y= y_LOD_pred, shape='LOD'),aes(x=x,y=y,shape = shape, guide=FALSE),colour="purple",size=5)
  
  LOQ_y = lower_Q_pred[which(abs(up_noise - lower_Q_pred) == min(abs(up_noise - lower_Q_pred)) )]
  p1 <- p1 +  geom_point(data=data.frame(x= LOQ_pred, y= y_LOQ_pred, shape='LOQ'),aes(x=x,y=y,shape = shape, guide=FALSE),colour=orange1,size=5)
  LOD_string <- paste('LOD=',round(LOD_pred, digits=1),sep='');   LOQ_string <- paste('LOQ=',round(LOQ_pred, digits=1),sep=''); 
  # p1 <- p1 +  scale_shape_manual( values=c(16,17), guide=FALSE)# ,labels=c(LOD_string, LOQ_string)) 
  
  p1 <- p1 + guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),     shape = guide_legend(order = 2)) + guides(shape=FALSE)
  #p1 <- p1 + xlim(0, LOQ_pred) + ylim(0,LOQ_y)
  
  #p1 <- p1 + annotate("text", x = 0.5*max(Cdata), y = 0.1*max(Idata), label = LOD_string, colour = "purple", size=5)
  #p1 <- p1 + annotate("text", x = 0.6*max(Cdata), y = 0.2*max(Idata), label = LOQ_string, colour = orange1, size=5)
  
  p1 <- p1 + ggtitle(paste(datain$NAME,': ', LOD_string,', ',LOQ_string,sep="")) + theme(plot.title = element_text(size = 20))
  
  print(p1);
  dev.off()
  
  
  
  
    
}