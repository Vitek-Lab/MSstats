
.bilinear <- function(x, intercept, slope, change) { 
    ifelse(x < change, intercept, slope*(x-change)+intercept) 
}

.second <- function(x,intercept,slope,slope2){
    slope2*x*x + slope*x + intercept
}

.bilinear_LOD <- function(x, intercept, slope, change) { 
    ifelse(x < change, intercept, slope*(x-change)+intercept) 
}

.bilinear_LOD_nonlin <- function(x, intercept, slope, change,exponent) { 
    ifelse(x < change, intercept, slope*(x-change)**exponent+intercept) 
}

.linear <- function(x,intercept,slope){
    slope*x + intercept
}

.invlinear <- function(y,intercept,slope){
    (y - intercept)/slope
}

.fun_error <- function(dir,peptide,site,message){
    filename = paste(dir,"/",paste(peptide,"_",site,".txt",sep=""),sep="")
    #print(filename)
    write(message,file=filename)
  
}

.ggplotColours <- function(n=6, h=c(0, 360) +15){
    if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


.theme_complete_bw <- function(base_size = 12, base_family = "") {
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(
            axis.line =         element_blank(),
            axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "black", vjust = 1, margin = margin(5,0,10,0)),
            axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1, margin = margin(0,5,0,10)),
            axis.ticks =        element_line(colour = "black"),
            axis.title.x =      element_text(size = base_size, vjust = 0.5),
            axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
            axis.ticks.length = unit(0.15, "cm"),
      
            legend.background = element_rect(colour=NA), 
            legend.key =        element_rect(fill = NA, colour = "white", size = 0.25),
            legend.key.size =   unit(1.2, "lines"),
            legend.text =       element_text(size = base_size * 0.8),
            legend.title =      element_text(size = base_size * 0.8, face = "bold", hjust = 0),
            legend.position =   "right",
      
            panel.background = element_rect(fill = "white", colour = NA), 
            panel.border =     element_rect(fill = NA, colour = "black"), 
            panel.grid.major = element_line(colour = "gray75", size = 0.1,  linetype = 'dotted'), 
            panel.grid.minor = element_line(colour = "gray75", size = 0.25,  linetype = 'dotted'), 
            panel.margin =     unit(0.25, "lines"),
      
            strip.background = element_rect(fill = NA, colour = NA), 
            strip.text.x =     element_text(colour = "black", size = base_size * 0.8),
            strip.text.y =     element_text(colour = "black", size = base_size * 0.8, angle = -90),
      
            plot.background =  element_rect(colour = NA, fill = "white"),
            plot.title =       element_text(size = base_size * 1.2,  margin = margin(0,0,10,0)),
            plot.margin =      unit(c(1, 1, 0.5, 0.5), "lines"))
}