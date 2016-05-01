
#The goal of this function is to calculate the LoD/LoQ of the data provided in the data frame.
#The function returns a new data frame containing the value of the LoD/LoQ

nonlinear_quantlim <- function(datain){
  
  
  #Need to rename variables as needed:
  
  names(datain)[names(datain) == 'CONCENTRATION'] <- 'C'
  names(datain)[names(datain) == 'INTENSITY'] <- 'I' 
  
  
  ##Coarse:
  B <- 500
  BB <- 200
  NMAX = 30
  Npoints <- 30
  Maxit = 3
 
  datain <- datain[!is.na(datain$I) & !is.na(datain$C),]  
  datain <- datain[!is.infinite(datain$I) & !is.infinite(datain$C),]  
  
  
  datain <- datain[order(datain$C),] 
  tmp_nob <- subset(datain,datain$C >0)  
  tmp_all <- datain
  tmp_blank <-subset(datain,datain$C == 0)
  
  #Calculate the value of the noise:
  #Use the zero concentration to calculate the LOD:
  
  noise = mean(tmp_blank$I)
  var_noise = var(tmp_blank$I)
  
  pb <- txtProgressBar(min = 0, max = 1, initial = 0, char = "%",
                       width = 40, title, label, style = 1, file = "")
  
  #Need to change this to account for limited number of samples (4):
  #The confidence interval for future observations of normal observations with unknown mean and variance:
  #t(alpha/2,dof = n-1)*s*sqrt(1+1/n)
  
  n_blank = length(unique(tmp_blank$I))

  fac_low = qt(1-0.25/2,n_blank - 1)*sqrt(1+1/n_blank)  
  fac = qt(1-0.05/2,n_blank - 1)*sqrt(1+1/n_blank)
  fac_high = qt(1-0.01/2,n_blank - 1)*sqrt(1+1/n_blank)
    
  #2.575
  
  low_noise = noise  - fac* sqrt(var_noise)
  
  up_noise_low = noise + fac_low * sqrt(var_noise)
  up_noise = noise   +  fac* sqrt(var_noise)
  up_noise_high = noise + fac_high*sqrt(var_noise)
  
  unique_c = sort(unique(tmp_all$C));   var_v <- rep(0.0, length(unique_c))
  weights  <- rep(0.0, length(tmp_all$C));
  weights_nob  <- rep(0.0, length(tmp_nob$C));
  
  #Calculate variance for all concentrations:
  ii = 1
  for (j in unique_c){
    data_f <- subset(tmp_all, C == j)
    var_v[ii] <- var(data_f$I)#/mean(data_f$log2Int)**2
    ii = ii +1
  }
  
  
  #Log scale discretization:
  xaxis_orig_2 <- exp(c( seq( from = log(10+0), to = log(1+max(unique_c)), by = log(1+max(unique_c))/Npoints ))) -10 #0.250 to go fast here
  xaxis_orig_2 <- unique(sort(c(xaxis_orig_2,unique_c))) 
  
   
  ###Loop here to make sure that the discretization is fine enough##
  
  x_new <- NULL
  for(it in 1:Maxit){
    

  
    if(!is.null(x_new)) xaxis_orig_2 <- c(xaxis_orig_2,x_new)
    
    y.loessm04 <- loess(y ~ x, span=0.6, data.frame(x=log(1+unique_c), y=log(1+var_v))) 
    y.predic_log <- predict(y.loessm04, data.frame(x=log(1+xaxis_orig_2)))
    y.predic_log_unique <- predict(y.loessm04, data.frame(x=log(1+unique_c)))
    
    #var_v_s_log <- exp(y.predic_log) -1
    #Can use this instead of the origianl variance for everything
    var_v_s_log <- pmax(exp(y.predic_log) -1, rep(1.0,length(y.predic_log)))
    var_v_s_log_unique <- pmax(exp(y.predic_log_unique) -1, rep(1.0,length(y.predic_log_unique)))  
    
    
    ##
    #In the following, consider that the smoothed variance is that from the log:
    
    var_v_s <- var_v_s_log
    var_v_s_unique <- var_v_s_log_unique
    ##
    
    
    if(1){# Full bootstrap
      
      outB <- matrix(NA_real_, nrow=B, ncol=length(xaxis_orig_2))
      outBB_pred <- matrix(NA_real_, nrow=B*BB, ncol=length(xaxis_orig_2))
      
      change_B <- rep(NA, B)
      for (j in 1:B){ #Number of first boostrap samples j
        
        
        setTxtProgressBar(pb, j/B, title = NULL, label = NULL)
        
        #if(j %% 10 == 0) print(paste("Boot=",j))
        lin.blank_B <- NULL
        tmpB <- tmp_all[sample(1:nrow(tmp_all), replace=TRUE),] #Pick **  observations with replacement among all that are available.
        #Blank samples are included
        
        weights = rep(0,length(tmpB$C) )
        for (kk in 1:length(tmpB$C)){
          weights[kk] = 1/var_v_s_unique[which( unique_c == tmpB$C[kk])]
        }   
        
        noise_B = mean(sample(tmp_blank$I,length(tmp_blank$I),replace=TRUE)) #Mean of resampled noise (= mean of noise)
        
        ii = 0;
        if(1)     while(ii < NMAX){#Number of ii trials for bilinear fit < NMAX
          ii = ii + 1
          {
          change = median(tmpB$C)*runif(1)*0.25
          slope=median(tmpB$I)/median(tmpB$C)*runif(1)
          
          sink(tempfile());
          #Set intercept at noise and solve for the slope and change
          fit.blank_B <- NULL
          fit.blank_B <- tryCatch({nlsLM( I ~ .bilinear_LOD(C , noise_B, slope, change),data=tmpB, trace = TRUE,start=c(slope=slope, change=change), weights = weights,
                                          control = nls.lm.control(nprint=1,ftol = sqrt(.Machine$double.eps)/2, maxiter = 50))}, error = function(e) {NULL}
          )
          
          
          sink();
          }
          
          out_loop = 0
          if(!is.null(fit.blank_B)){ #Converges but cannot have a real threshold here anyway
            if(  summary(fit.blank_B)$coefficient[2] < min(tmpB$C)  ){
              fit.blank_B <- NULL
              out_loop =1
            }
          }
          
          
          if(!is.null(fit.blank_B) && out_loop == 0){ #Boostrap again to see whether bilinear fit is real:
            
            change_BB <- rep(NA,NMAX)
            #print(paste('Change=',summary(fit.blank_B)$coefficient[2]))
            bb = 0;
            while(bb < NMAX){#Number of second boostrap samples bb < NMAX
              bb = bb + 1
              
              iii = 0
              while(iii < NMAX){#Number of iii trials for convergence of bilinear
                iii = iii + 1
                tmpBB <- tmpB[sample(1:nrow(tmpB), replace=TRUE),] 
                change = median(tmpBB$C)*runif(1)*0.25
                slope=median(tmpBB$I)/median(tmpBB$C)*runif(1)
                
                weightsB = rep(0,length(tmpB$C) )
                for (kk in 1:length(tmpBB$C)){
                  weightsB[kk] = 1/var_v_s_unique[which( unique_c == tmpBB$C[kk])]
                }  
                
                
                #Need to also bootstrap for the value of the mean:
                #Pick with replacement blank samples:
                noise_BB = noise
                
                sink(tempfile());
                
                
                fit.blank_BB <- NULL
                fit.blank_BB <- tryCatch({nlsLM( I ~ .bilinear_LOD(C , noise_BB, slope, change),data=tmpBB, trace = TRUE,start=c(slope=slope, change=change), weights = weightsB,
                                                 control = nls.lm.control(nprint=1,ftol = sqrt(.Machine$double.eps)/2, maxiter = 50))}, error = function(e) {NULL}
                )
                
                sink();
                
                if(!is.null(fit.blank_BB)){ 
                  change_BB[bb] = summary(fit.blank_BB)$coefficient[2]
                }else{
                  change_BB[bb] = NA 
                }
                
                
                
                if(!is.null(fit.blank_BB)) break
                
              } #Number of iii trials for convergence of bilinear
              
              #if(is.null(fit.blank_BB)) print(paste('No convergence for coefficient ',j,bb))
              
            }#Number of second boostrap samples bb < NMAX
            
            CI_change <- quantile(change_BB,probs=c(0.05,0.95),na.rm= TRUE) #95% Confidence interval for the value of change
            #Ensure that the 95% confidence interval is included inside the concentration range:
            if(is.na(CI_change[1]) || is.na(CI_change[2])){
              fit.blank_B  <- NULL
              out_loop =1
            }
            
            
            if(!is.na(CI_change[1]) &&  !is.na(CI_change[2])) if(CI_change[1] < min(tmp_all$C) | CI_change[2] > max(tmp_all$C)){
              fit.blank_B  <- NULL
              out_loop =1
            }else{
              #fit.blank_B  <- NULL#print('Acceptable fit')
            }
            
          } #Boostrap again to see whether bilinear fit is real
          
          if(out_loop == 1) break #Bilinear fit converges but CI not acceptable
          if(!is.null(fit.blank_B)) break #Could never find a converged bilinear fit
        } #Number of ii trials for bilinear fit < NMAX
        
        #fit.blank_B <- NULL
        if(is.null(fit.blank_B)){ #Do linear fit:
          ll = 0
          while(ll < NMAX){
            ll = ll + 1
            slope = median(tmpB$I)/median(tmpB$C)*runif(1)
            intercept = noise*runif(1)
            sink(tempfile());
            lin.blank_B <-  tryCatch({nlsLM( I ~ .linear(C , intercept, slope),data=tmpB, trace = TRUE,start=c(intercept=intercept, slope=slope), weights = weights,
                                             control = nls.lm.control(nprint=1,ftol = sqrt(.Machine$double.eps)/2, maxiter = 50))}, error = function(e) {NULL}
            )
            sink();
            if(!is.null(lin.blank_B)) break
          }
        } #Do linear fit if is.null(fit.blank_B)
        
        
        #Store the curve fits obtained via bootstrap with bilinear and linear:
        if (!is.null(fit.blank_B)){
          outB[j,]  <- .bilinear_LOD(xaxis_orig_2, noise_B, summary(fit.blank_B)$coefficient[1] , summary(fit.blank_B)$coefficient[2])           
          change_B[j] <- summary(fit.blank_B)$coefficient[2]
          
        } else{
          if(!is.null(lin.blank_B)){#If linear fit, change = 0 anyway
            outB[j,]  <- .linear(xaxis_orig_2,summary(lin.blank_B)$coefficient[1] , summary(lin.blank_B)$coefficient[2] )
            change_B[j] <- 0
          }
          else{
            outB[j,]  <- rep(NA, length(xaxis_orig_2))
            change_B[j] <- NA
          }
        }
        
        for (jj in 1:BB){#Calculate predictions
          
          
          outBB_pred[(j-1)*BB +jj,] <- outB[j,] + rnorm( length(xaxis_orig_2),0,sqrt(var_v_s_log))
          
        }
        
      } # Number of first bootstrap samples j <=B
      
      
      #Calculate the variance of the fits: 
      var_bilinear <- apply(outB, 2, var,na.rm = TRUE)
      mean_bilinear <- apply(outB, 2, mean,na.rm = TRUE)
      
      #Calculate confidence interval based on quantiles:
      # Do not use: based on normal distribution
      ##lower_CI = mean_bilinear  - 1.96* sqrt(var_bilinear)
      ##upper_CI = mean_bilinear  + 1.96* sqrt(var_bilinear)
      
      
      mean_pred <- apply(outBB_pred, 2, mean,na.rm = TRUE) 
      var_pred <- apply(outBB_pred, 2, var,na.rm = TRUE) 
      
      lower_Q = apply(outB, 2, quantile, probs=c(0.05) ,na.rm = TRUE)
      upper_Q = apply(outB, 2, quantile, probs=c(0.95) ,na.rm = TRUE)
      
      lower_Q_pred = apply(outBB_pred, 2, quantile, probs=c(0.05) ,na.rm = TRUE)
      upper_Q_pred = apply(outBB_pred, 2, quantile, probs=c(0.95) ,na.rm = TRUE)    
      
      
      
      
    }#Full bootstrap method

    
    LOD_pred_mean_low = 0; LOD_pred_high =0;
    #Calculate the LOD/LOQ from the prediction interval:
    i_before = which(diff(sign( up_noise - mean_bilinear  ))!=0) #before sign change
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise - mean_bilinear)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise - mean_bilinear)[i_after]
      
      #Linear interpollation to find where the function changes sign:
      LOD_pred_mean =  x1  - f1*(x2-x1)/(f2-f1)
      LOD_pred = LOD_pred_mean
      y_LOD_pred = up_noise
    } else{ 
      LOD_pred = 0
      y_LOD_pred = up_noise
    }
    
    #Calculate the LOD with the upper-upper and upper-lower limits of the noise (Calculated to make sure that we have a large enough resolution)
    
    i_before = which(diff(sign( up_noise_low - mean_bilinear  ))!=0) #before sign change
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_low - mean_bilinear)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_low - mean_bilinear)[i_after]
      
      #Linear interpollation to find where the function changes sign:
      LOD_pred_mean_low =  x1  - f1*(x2-x1)/(f2-f1)
      LOD_pred_low = LOD_pred_mean_low
      y_LOD_pred_low = up_noise
    } else{ 
      LOD_pred_low = 0
      y_LOD_pred_low = up_noise
    }  
    
    i_before = which(diff(sign( up_noise_high - mean_bilinear  ))!=0) #before sign change
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_high - mean_bilinear)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_high - mean_bilinear)[i_after]
      
      #Linear interpollation to find where the function changes sign:
      LOD_pred_mean_high =  x1  - f1*(x2-x1)/(f2-f1)
      LOD_pred_high = LOD_pred_mean_high
    } else{ 
      LOD_pred_mean_high = 0
    }      
    
    
    i_before = which(diff(sign( up_noise_low - mean_bilinear  ))!=0) #before sign change
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_low - mean_bilinear)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_low - mean_bilinear)[i_after]
      
      #Linear interpollation to find where the function changes sign:
      LOD_pred_mean_low =  x1  - f1*(x2-x1)/(f2-f1)
      LOD_pred_low = LOD_pred_mean_low
      y_LOD_pred_low = up_noise
    } else{ 
      LOD_pred_low = 0
      y_LOD_pred_low = up_noise
    }  
    
    i_before = which(diff(sign( up_noise_high - mean_bilinear  ))!=0) #before sign change
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_high - mean_bilinear)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_high - mean_bilinear)[i_after]
      
      #Linear interpollation to find where the function changes sign:
      LOD_pred_mean_high =  x1  - f1*(x2-x1)/(f2-f1)
      LOD_pred_high = LOD_pred_mean_high
    } else{ 
      LOD_pred_mean_high = 0
    }  
    
    



    
    
    
    
    #Do a linear fit to find the intersection:
    i_before = which(diff(sign(up_noise - lower_Q_pred))!=0)
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise - lower_Q_pred)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise - lower_Q_pred)[i_after]
      x_inter =  x1  - f1*(x2-x1)/(f2-f1)
      LOQ_pred = x_inter
      y_LOQ_pred = up_noise
    } else{
      LOQ_pred = 0
      y_LOQ_pred = up_noise
      
    }
    if(length(LOD_pred) > 1) print('multiple intersection between fit and upper bound of noise, picking first')
    LOQ_pred=LOQ_pred[1]
    
    
    i_before = which(diff(sign( up_noise_low - mean_bilinear  ))!=0) #before sign change
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_low - mean_bilinear)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_low - mean_bilinear)[i_after]
      
      #Linear interpollation to find where the function changes sign:
      LOD_pred_mean_low =  x1  - f1*(x2-x1)/(f2-f1)
      LOD_pred_low = LOD_pred_mean_low
      y_LOD_pred_low = up_noise
    } else{ 
      LOD_pred_low = 0
      y_LOD_pred_low = up_noise
    }  
    
    i_before = which(diff(sign( up_noise_high - mean_bilinear  ))!=0) #before sign change
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_high - mean_bilinear)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_high - mean_bilinear)[i_after]
      
      #Linear interpollation to find where the function changes sign:
      LOD_pred_mean_high =  x1  - f1*(x2-x1)/(f2-f1)
      LOD_pred_high = LOD_pred_mean_high
    } else{ 
      LOD_pred_mean_high = 0
    }  
    
    
    #Do a linear fit to find the intersection:
    i_before = which(diff(sign(up_noise - lower_Q_pred))!=0)
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_low - lower_Q_pred)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_low - lower_Q_pred)[i_after]
      x_inter =  x1  - f1*(x2-x1)/(f2-f1)
      LOQ_pred_low = x_inter
    } else{
      LOQ_pred_low = 0
    }
    #Do a linear fit to find the intersection:
    i_before = which(diff(sign(up_noise - lower_Q_pred))!=0)
    if(length(i_before)>0){
      i_after = i_before+1
      x1 = xaxis_orig_2[i_before]; f1 = (up_noise_high - lower_Q_pred)[i_before]
      x2 = xaxis_orig_2[i_after];  f2 = (up_noise_high - lower_Q_pred)[i_after]
      x_inter =  x1  - f1*(x2-x1)/(f2-f1)
      LOQ_pred_high = x_inter
    } else{
      LOQ_pred_high = 0
    }    
    
    
    #Want to require at least one point between the various LOD's
    
    found_low_LOD = 0; found_high_LOD = 0; x_low = -7; x_high = -7
    found_low_LOQ = 1; found_high_LOQ = 1;
    for( x in xaxis_orig_2 ){
      
      if( x < LOD_pred_mean_high &&  x > LOD_pred){
        found_high_LOD =1 ;x_high = x;
      }
      
      if( x > LOD_pred_mean_low &&  x < LOD_pred){
        found_low_LOD =1 ; x_low = x;
      }
      
      
      if( x < LOQ_pred_high &&  x > LOQ_pred){
        found_high_LOQ =1 ;x_high = x;
      }
      
      if( x > LOQ_pred_low &&  x < LOQ_pred){
        found_low_LOQ =1 ; x_low = x;
      }
    
      }
    
    
    
    x_new <- NULL
    
    if(found_low_LOD == 0){
      x_new =  c(x_new, 0.5*(LOD_pred_low + LOD_pred))
    }
    if(found_high_LOD == 0){
      x_new = c(x_new, 0.5*(LOD_pred_high + LOD_pred))
    }
    if(found_low_LOQ== 0){
      x_new =  c(x_new, 0.5*(LOQ_pred_low + LOQ_pred))
    }
    if(found_high_LOQ == 0){
      x_new = c(x_new, 0.5*(LOQ_pred_high + LOQ_pred))
    }   
    
    
    
    if(is.null(x_new)) break;
    
  }
    
    
    
    #if(is.null(x_new)) break
    
    
    return(
      data.frame(as.data.frame(list(CONCENTRATION = xaxis_orig_2, MEAN=mean_bilinear,LOW= lower_Q_pred, UP = upper_Q_pred, LOD= rep(LOD_pred, length(upper_Q_pred)),  LOQ = rep(LOQ_pred, length(upper_Q_pred)),
                NAME = rep(datain$NAME[1], length(upper_Q_pred)),
                METHOD = rep("NONLINEAR", length(upper_Q_pred))
                )
      ))
           )
}