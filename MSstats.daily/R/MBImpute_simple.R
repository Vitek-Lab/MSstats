# Model-based imputaion

.MBimpute <- function(m, my.pi=0.05, compute_pi=TRUE){
# Imputes missing values 
# Expects the data to be filtered with model-based filtering prior to imputation
# Handles no sibling peptides with missing groups (think not this version)
# 
# Input:
#   input: An m x n matrix of intensities
#          DEPRICATED:(peptides x (n+2 samples)) matrix of expression data
#          where first column contains peptide identifier and second column
#          contains protein identifier
#   treatment:  vector indicating the treatment group of each sample ie [1 1 1 1 2 2 2 2...]
#   my.pi: PI value, estimate of the proportion of peptides missign completely at random,
#          as compared to censored at lower abundance levels
#  
# Output: list of:
#   y.impute: k x n matrix of 'complete' data (original and imputed values for missing data),
#             k <= m, as some peptides may have been filtered out due to low information content
#
# p-values computed after model-based filtering and imputation


# calculate PI or used one passed in as parameter:
  if (compute_pi){
  	my.pi <- .eigen_pi(m, toplot=T)
  } # else should be taken from the user interface, can we verify it?


	prot.info<-m[,2:1]

  # run factors
  n.run <- ncol(m)-2

  # Match to protein
  all.proteins <- unique(prot.info[,2])
  all.proteins <- all.proteins[order(all.proteins)]


  y_imputed <- NULL
  name<-NULL
  
  cat("Imputing...\n")

  for (k in 1:length(all.proteins)){
    y_raw<-m[which(m$PROTEIN %in% all.proteins[k]),c(3:ncol(m))]
    
    if (nrow(y_raw) == 0) next
    
    # re-evaluate data parameters after filtering
    n.peptide <- nrow(y_raw)
    y <- as.vector(t(y_raw))
    n <- length(y)

	### minimum value in this protein
    c.guess <- min(y, na.rm=T)
    peptide <-rep(rep(1:n.peptide, each=n.run))

    # re-calculate number of observed values per treatment group for each peptide
    # yuliya: may be able to use previous computation here for efficiency... but as long
    # as this works do not care to change at this time

    # calculate pooled variance for each protein
    grp<-apply(y_raw,2, function(x) sum(!is.na(x)))
    mpos <-which.max(grp)
    
    go <- .protein_var(y_raw) # yuliya: function, see below : here used group information
    overall_var <- go$overall_var

  #  num <-0
  #  den <- 0

    # calculate pooled variance for each peptide;
    # if only 1 onservation in a dx group assign the overall variance
    # for(i in 1:n.peptide){
      # y.i <- na.omit(y[peptide==i & treatment==unique(treatment)[mpos]])
      # p2 <- var(y.i)
      # if (is.na(p2)) p2 <- 0
      # present <- length(y.i)
      # num <- num+(p2*(present-1))
      # den <- den + (present-1)
    # }
    # pep_var <- num/den
    # if (pep_var==0) pep_var <- overall_var
    
    ## !! need to check : one obsat each feature and run. Therefore, can't calculate variance
	pep_var <- overall_var
	
    peptides.missing <- rowSums(is.na(y_raw))

    f.run <- factor(rep(colnames(y_raw), n.peptide))
    f.peptide <- factor(peptide)

    # estimate rough model parameters
    # create model matrix for each protein and
    # remove any peptides with missing values
    ii <- (1:n)[is.na(y)]
    
    if (n.peptide != 1){

  	## ~f.peptide+f.treatment
  	## however, f.treatment should be f.run
  	## however, one observation per run. singular matrix
  	
      X  <- model.matrix(~f.peptide, contrasts = list(f.peptide="contr.sum") )
    } else {
    
      ## !! need to check : f.treatment?
      X <- model.matrix(~f.peptide, contrasts=list(f.peptide="contr.sum"))
    }
    if(length(ii) > 0){
      y.c <- y[-ii]
      X.c <- X[-ii,]
    } else {
      y.c <- y
      X.c <- X
    }

    # calculate initial beta values and residuals
    beta <- drop(solve(t(X.c) %*% X.c) %*% t(X.c) %*% y.c)

    # compute initial delta's
    peptides.missing[peptides.missing==0] <- 0.9 ##?? if not, Inf for delta.y?
    delta.y <- as.numeric(1/sqrt(pep_var*peptides.missing))
    dd <- delta.y[as.numeric(peptide)]


    # calculate cutoff values for each peptide
    c_hat = rep(NA, n.peptide)
    for(j in 1:n.peptide) {
      c_hat[j] = min(y_raw[j, ], na.rm = T)
    }
    c_h <- c_hat[as.numeric(peptide)]

    if(n.peptide==1){ ## need to check
      y.predict <- model.matrix(~f.peptide,contrasts = list(f.peptide="contr.sum"))%*% beta
    } else {
      y.predict <- model.matrix(~f.peptide,contrasts = list(f.peptide="contr.sum"))%*% beta
    }

    zeta <- dd*(c_h - y.predict)
    prob.cen <- pnorm(zeta, 0, 1)/(my.pi + (1-my.pi)*pnorm(zeta, 0, 1))
    choose.cen <- runif(n) < prob.cen
    set.cen <- is.na(y) &choose.cen
    set.mar <- is.na(y) &!choose.cen
    kappa <- my.pi + (1 - my.pi)*dnorm(zeta,0, 1)

    # compute information
   # I_beta <- t(X) %*% diag(as.vector(dd^2*(1 - kappa*(1 + my.Psi.dash(zeta, my.pi))))) %*% X
  #  I_GRP <- I_beta[rev(rev(1:n.peptide+1)[1:(n.u.treatment - 1)]), rev(rev(1:n.peptide+1)[1:(n.u.treatment - 1)]), drop = F]
 #   if (!is.null(dim(I_GRP))) I_GRP = det(I_GRP)

    # Imputation: Replace missing values with random numbers drawn from the estimated likelihood model
    sigma <- 1/dd
    y.impute <- t(y_raw)
    if(sum(set.cen) > 0)
      y.impute[set.cen] <- .rnorm.trunc(sum(set.cen), y.predict[set.cen], sigma[set.cen], hi=rep(c.guess, n)[set.cen])

    if(sum(set.mar) > 0)
      y.impute[set.mar] <- rnorm(n, y.predict, sigma)[set.mar]
      
    y.impute.return <- t(y.impute)

    name<-prot.info[prot.info$PROTEIN==all.proteins[k],]
	
	y.impute.return<-cbind(name,y.impute.return)
	
    #row.names(y.impute.ret) <- row.names(y_raw)
    tryCatch(y_imputed <- rbind(y_imputed, y.impute.return), error = function(e){
      print(y.impute.return)
      print(dim(y_imputed))
      print(dim(y.impute.return))      
      stop("Error in rbind() for y.impute.net")      
    })
    
  } #end for each protein
  

  cat("Done imputing.\n")
  return(list(y_imputed=y_imputed,pi=as.matrix(my.pi)))
}


############# end imputation ###################


############## function pi #############
.eigen_pi <- function(m, toplot=T){ 
# Compute PI - proportion of observations missing at random
# INPUT: m - matrix of abundances, numsmaples x numpeptides
#       toplot - T/F plot mean vs protportion missing, curve and PI 
# OUTPUT: pi - 
# 
# Shelley Herbrich, June 2010, created for EigenMS and TamuQ 
#
# (1) compute 1) ave of the present values from each petide
#             2) number of missing and present values for each peptide
  
  m = m[,-(1:2)] # yuliya: no 2 columns of ids, already on log scale
  # m = log(m)
  #remove completely missing rows
  m = m[rowSums(m, na.rm=T)!=0,]

  pepmean <- apply(m, 1, mean, na.rm=T)
  propmiss <- rowSums(is.na(m))/ncol(m)

  smooth_span <- (0.4)
  fit <- lowess(pepmean, propmiss, f=smooth_span)
  PI <- fit$y[fit$x==max(pepmean)]

  count <- 1
  while (PI<=0){
    smooth_span <- smooth_span-.1
    fit <- lowess(pepmean, propmiss, f=smooth_span)
    PI <- fit$y[fit$x==max(pepmean)]
    count <- count + 1
    if (count > 500) break
  }

  if (toplot){
  st <- paste("PI: ", PI) 
  plot(pepmean, propmiss, xlab="x", ylab="y", cex=0.5) #plot data point
  lines(fit)
  title("Lowess Regression", sub = st,
      cex.main = 2,   font.main= 3, col.main= "purple",
      cex.sub = 1, font.sub = 3, col.sub = "red")
  }
  return (pi=PI)
}





######################################################
.protein_var <- function(Y_raw){
  n.peptide <- nrow(Y_raw)
  y <- as.vector(t(Y_raw))
  n <- length(y)
  n.run <- ncol(Y_raw)
#  n.u.treatment <- length(unique(treatment))
  peptide <-rep(rep(1:n.peptide, each=n.run))

  n.present<-ifelse(!is.na(Y_raw),1,0)
  peptides.missing <- rowSums(is.na(Y_raw))

  f.run <- factor(rep(colnames(Y_raw), n.peptide))
  f.peptide <- factor(peptide)

  # estimate rough model parameters
  # create model matrix for each protein and
  # remove any peptides with missing values
  ii <- (1:n)[is.na(y)]
  if (n.peptide != 1){
  	## ~f.peptide+f.treatment
  	## however, f.treatment should be f.run
  	## however, one observation per run. singular matrix
    X  <- model.matrix(~f.peptide, contrasts = list(f.peptide="contr.sum") )
  } else {
  	## !! need to check
    X <- model.matrix(~f.peptide, contrasts=list(f.run="contr.sum"))
  }
  
  if(length(ii) > 0){
    y.c <- y[-ii] ## remove missing value
    X.c <- X[-ii,]
  } else {
    y.c <- y
    X.c <- X
  }

  # calculate initial beta values and residuals
  beta <- drop(solve(t(X.c) %*% X.c) %*% t(X.c) %*% y.c)
  Y_hat <- X.c %*% beta
  Y_temp <- Y_raw
  Y_temp <- as.numeric(t(Y_temp))
  Y_temp[!is.na(Y_temp)] <- Y_hat
  Y_temp <- matrix(Y_temp, nrow = n.peptide, byrow = T)
  Y_hat <- Y_temp
  ee <- Y_raw - Y_hat

  effects <- X.c %*% beta
  resid <- y.c - effects
  overall_var <- var(resid)
  return(list(overall_var=det(overall_var)))
}


################################################

.rnorm.trunc <- function (n, mu, sigma, lo=-Inf, hi=Inf){
# Calculates truncated noraml
  p.lo <- pnorm (lo, mu, sigma)
  p.hi <- pnorm (hi, mu, sigma)
  u <- runif (n, p.lo, p.hi)
  return (qnorm (u, mu, sigma))
}
# End rnorm.trunc





