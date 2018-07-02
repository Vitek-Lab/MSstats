#' Estimate the optimal size of training data for classification problem
#'
#' @param data output from function \code{\link{dataProcess}}
#' @param n_sample number of different sample size to simulate. Default is 5 
#' @param sample_incr number of samples per condition to increase at each step. Default is 20 
#' @param protein_desc the fraction of proteins to reduce at each step. Proteins are ranked based on their mean abundance across all the samples. Default is 0.2. If protein_desc = 0.0, protein number will not be changed. 
#' @param iter number of times to repeat simulation experiments. Default is 10
#' @description For classification problem (such as disgnosys of disease), calculate the mean predictive accuray under different size of training data for future experiments of a Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment based on simulation.
#' @details The function fits intensity-based linear model on the input prelimiary data \emph{data} and uses variance components and mean abundance to simulate new training data with different sample size and protein number. Random forest model is fitted on simulated train data and used to predict the input preliminary data \emph{data}. The above procedure is repeated \emph{iter} times. Mean predictive accuracy and variance under different size of training data are reported. 
#' @return \emph{meanPA} is the mean predictive accuracy matrix under different size of training data.
#' @return \emph{varPA} is variance of predictive accuracy under different size of training data.
#' @author Ting Huang, Meena Choi, Olga Vitek.
#'
#' Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
#' @references T. Huang et al.  TBD  2018
#' @examples # Consider the training set from a colorectal cancer study
#' # Subjects are from control group or colorectal cancer group
#' # 72 proteins were targeted with SRM
#' require(MSstatsBioData)
#' set.seed(1235)
#' data(SRM_crc_training)
#' QuantCRCSRM <- dataProcess(SRM_crc_training, normalization = FALSE)
#' # estimate the mean predictive accuray under different sizes of training data
#' # n_sample is the number of different sample size to simulate
#' # Datasets with 10 different sample size and 3 different protein numbers are simulated 
#' result.crc.srm <- designSampleSizeClassification(data=QuantCRCSRM, 
#' n_sample = 10, 
#' sample_incr = 10, 
#' protein_desc = 0.33, 
#' iter = 50)
#' result.crc.srm$meanPA # mean predictive accuracy
#' @importFrom randomForest randomForest combine
#' 
#' @export
designSampleSizeClassification <- function(data, 
                                           n_sample = 5, 
                                           sample_incr = 20, 
                                           protein_desc = 0.2,
                                           iter = 10) {
  
    ## save process output in each step
    allfiles <- list.files()
    filenaming <- "msstats"
    
    if (length(grep(filenaming, allfiles)) == 0) {
    
        finalfile <- "msstats.log"
        processout <- NULL
    
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
                        as.matrix(c(" ", " ", "MSstats - designSampleSizeClassification function", " "), ncol=1))
    
    ## check input is correct
    ## data format
    rawinput <- c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                  "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                  "Condition", "BioReplicate", "Run", "Intensity")
    ## check data
    if (length(setdiff(toupper(rawinput),toupper(colnames(data$ProcessedData)))) == 0) {
        processout <- rbind(processout,
                            "The required input - data : did not process from dataProcess function. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison.")
    }
    
    ## check n_sample
    if (!((n_sample>0) & n_sample%%1==0)) {
        processout <- rbind(processout,
                            "The required input - n_sample : was not a positive integer. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set n_sample as a positive integer.")
    }
    
    ## check sample_incr
    if (!((sample_incr>0) & sample_incr%%1==0)) {
        processout <- rbind(processout,
                            "The required input - sample_incr : was not a positive integer. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set sample_incr as a positive integer.")
    }
    
    ## check protein_desc
    if (!(is.numeric(protein_desc) & (protein_desc >= 0) & (protein_desc < 1))) {
        processout <- rbind(processout,
                            "The required input - protein_desc : was not a numeric value between 0 and 1 - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set protein_desc as a numeric value.")
    }
    
    ## check iteration
    if (!((iter>0) & iter%%1==0)) {
        processout <- rbind(processout,
                            "The required input - iter : was not a positive integer. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set iter as a positive integer.")
    }
    
    processout <- rbind(processout, paste0("n_sample = ", n_sample))
    processout <- rbind(processout, paste0("sample_incr = ", sample_incr))
    processout <- rbind(processout, paste0("protein_desc = ", protein_desc))
    processout <- rbind(processout, paste0("iter = ", iter))
    
    ## Estimate the mean abundance and variance of each protein in each phenotype group
    parameters <- .estimateVar(data)
    processout <- rbind(processout, 
                        "Calculated variance component. - okay")
    write.table(processout, file=finalfile, row.names=FALSE)
    
    ## Prepare the parameters for simulation experiment
    mu <- parameters$mu
    sigma <- parameters$sigma
    promean <- parameters$promean
    valid_x <- as.data.frame(parameters$X)
    valid_y <- as.factor(parameters$Y)

    ## Generate the vector of training sample size to simulate
    ngroup <- length(unique(data$ProcessedData$GROUP_ORIGINAL)) # Number of phenotype groups
    train_size <- seq.int(from = sample_incr, to = sample_incr * n_sample, length.out = n_sample)
    train_size <- train_size*ngroup
    message(" Size of training data to simulate: ", paste(train_size, collapse = ', '))
    
    ## Generate the vector of protein number to simulate
    nproteins <- nrow(mu)
    if(protein_desc == 0.0){ # no change of protein number
        protein_num <- nproteins
    } else{ # decrease the protein number by nproteins*protein_desc each time
        m_prot <- round(nproteins*protein_desc)
        protein_num <- seq.int(from = m_prot, to = nproteins, by = m_prot)
    }
     message(" Number of proteins to simulate: ", paste(protein_num, collapse = ', '))
    
    processout <- rbind(processout, 
                        "Prepare simulation paramters. - okay")
    write.table(processout, file=finalfile, row.names=FALSE)
    
    message(" Start to run the simulation...")
    
    PA <- list()
    for (i in 1:iter) { ## Number of iterations
        message("  Iteration: ", i)
        accur <- matrix(rep(0, times=length(train_size) * length(protein_num)), nrow=length(protein_num))
        
        ## simulate train data with different size
        for (n in seq_along(protein_num)) { 
            message("    Protein Number: ", protein_num[n])
            
            ## select proteins based on their mean abundance
            selectedPros<-order(promean,decreasing = TRUE)[1:protein_num[n]]
            mu_2 <- mu[selectedPros,]
            sigma_2 <-sigma[selectedPros,]

            ## Retired code: Simulate an independant validation data
            # valid <- .sampleSimulation(valid_size, mu_2, sigma_2) 
            # valid_x <- as.data.frame(valid$X)
            # colnames(valid_x) <- rownames(mu)
            # valid_y <- as.factor(valid$Y)
            
            for (m in seq_along(train_size)) { ## simulate samples in the training data
                ##Simulate training data 
                train <- .sampleSimulation(train_size[m], mu_2, sigma_2) 
                x <- as.data.frame(train$X)
                colnames(x) <- rownames(mu_2)
                y <- as.factor(train$Y)
                
                #Train random forest on training data
                rf <- randomForest::randomForest(x=x, y=y, mtry = )
                
                ## Retired code: parallel computing for random forest
                #registerDoSNOW()
                #mcoptions <- list(set.seed=FALSE)
                #rf <- foreach(ntree=rep(10, 10), .combine=randomForest::combine, .multicombine=TRUE,
                          #.packages='randomForest', .options.multicore=mcoptions) %dopar% {
                             # randomForest(x=x,
                                           #y=y,
                                           #ntree=ntree)}
                
                ## Calculate predictive accuracy on validation data
                rf.pred <- predict(rf, valid_x) #Predict validation data
                accuracy <- sum(diag(table(rf.pred,valid_y))) / length(rf.pred)
                accur[n,m] <- accuracy
            }
        }
        PA[[i]] <- accur
    }
    
    ## Calculate the mean accuracy and variance
    meanPA <- matrix(rep(0, times=length(train_size) * length(protein_num)), nrow=length(protein_num))
    varPA <- matrix(rep(0, times=length(train_size) * length(protein_num)), nrow=length(protein_num))
    
    for (n in seq_along(protein_num)) {
        for (m in seq_along(train_size)) {
            temp <- NULL
            for (i in 1:iter) {
                temp <- c(temp, PA[[i]][n, m])
            }
            meanPA[n, m] <- mean(temp)
            varPA[n, m] <- var(temp)
        }
    }
    rownames(meanPA) <- paste0("prot", protein_num)
    colnames(meanPA) <- paste0("tra", train_size)
    rownames(varPA) <- paste0("prot", protein_num)
    colnames(varPA) <- paste0("tra", train_size)
    
    ##
    processout <- rbind(processout, 
                        "The number of sample and proteins is estimated - okay")
    write.table(processout, file=finalfile, row.names=FALSE)
    
    return(list(meanPA = meanPA, 
                varPA = varPA))
}
    
    

#' (1) Estimate the mean abundance and variance of each protein in each phenotype group. (2) Estimate the protein abundance of the input data
#'
#' @param data The run level data reported by dataProcess function.
#' @return \emph{mu} is the mean abundance matrix of each protein in each phenotype group; 
#' @return \emph{sigma} is the sd matrix of each protein in each phenotype group;
#' @return \emph{promean} is the mean abundance of each protein across all the samples.
#' @return \emph{X} is the protein abundance matrix for \emph{data}
#' @return \emph{Y} is group information vector for \emph{data}
#' @keywords internal
.estimateVar <- function(data) {
    
    ## get the list of proteins in the data
    rqall <- data$RunlevelData
    rqall$Protein <- factor(rqall$Protein)
    proteins <- levels(rqall$Protein)
    
    ## generate a fake contrast matrix
    groups <- unique(data$ProcessedData$GROUP_ORIGINAL)
    comparison <- matrix(rep(0, length(groups)), nrow=1)
    comparison[1, 1] <- 1
    comparison[1, 2] <- -1
    row.names(comparison) <- paste0(groups[1], groups[2])
    
    ## Estimate variance components and mean abundances in each group 
    ResultOneComparison <- groupComparison(contrast.matrix=comparison, data=data)
    res <- ResultOneComparison$fittedmodel
    
    ## Store the results
    VarComponent <- data.frame(Error=NA, Subject=NA, GroupBySubject=NA)
    GroupVar <- matrix(rep(0.0, length(res) * length(groups)), ncol = length(groups))
    GroupMean <- matrix(rep(0.0, length(res) * length(groups)), ncol = length(groups))
    SampleMean <- NULL # mean across all the samples
    count <- 0
    rnames <- NULL
    for (i in 1:length(res)) {
    
        # note: when run is fixed, we can obtain the same variance of error for both case-control and time course studies.
        fit.full <- res[[i]]
        
        ## if fit.full==NA (class(fit.full)=="try-error)
        if (is.null(fit.full)) {
            ## !!!!! but if we have NULL for last protein?
            next
          
        } else {
          
            ## get variance component
            if (class(fit.full) != "lmerMod") {
                VarComponent[i, "Error"] <- summary(fit.full)$sigma^2
            } else {
                stddev  <-  c(sapply(VarCorr(fit.full), function(el) attr(el, "stddev")),attr(VarCorr(fit.full), "sc"))
                VarComponent[i, "Error"] <- stddev[names(stddev) == ""]^2
                    
                if (sum(names(stddev) %in% "SUBJECT_NESTED.(Intercept)") > 0) {
                    VarComponent[i, "Subject"] <- stddev[names(stddev) == "SUBJECT_NESTED.(Intercept)"]^2
                }
                if (sum(names(stddev) %in% "SUBJECT.(Intercept)") > 0) {
                    VarComponent[i, "Subject"] <- stddev[names(stddev) == "SUBJECT.(Intercept)"]^2
                }
                if (sum(names(stddev) %in% "SUBJECT:GROUP.(Intercept)") > 0) {
                    VarComponent[i, "GroupBySubject"] <- stddev[names(stddev) == "SUBJECT:GROUP.(Intercept)"]^2
                }
            }
            ## get mean abundance in each group
            abun <- summary(fit.full)$coefficients[, 1]
            abun[-1] <- abun[1] + abun[-1]
            if(length(abun) == length(groups)){
                count <- count + 1
                rnames <- c(rnames, proteins[i])
                GroupMean[count, ] <- abun # group mean
                GroupVar[count, ] <- rep(sqrt(sum(VarComponent[i, ], na.rm = TRUE)), times=length(abun)) # group variance
                if (class(fit.full) == "lm") {
                    SampleMean <- c(SampleMean, mean(fit.full$model$ABUNDANCE, na.rm = T))
                } else{ # class(fit.full) == "lmerMod"
                    SampleMean <- c(SampleMean, mean(fit.full@frame$ABUNDANCE, na.rm = T))
                }
            } else{
                message("Protein ", proteins[i], " is missing from at least one group and not considered in the simulation.")
            }
        }  
    } ## end-loop
    
    GroupMean <- GroupMean[1:count, ] # remove the rows for unqualitied proteins
    GroupVar <- GroupVar[1:count, ] # remove the rows for unqualitied proteins
    rownames(GroupMean) <- rnames
    rownames(GroupVar) <- rnames
    names(SampleMean) <- rnames
    
    # get the protein abundance matrix of the input preliminary data
    Quant <- quantification(data)
    rownames(Quant) <- Quant[,1]
    Quant <- Quant[,-1] # columns are samples and rows are proteins
    valid_x <- Quant # generate data matrix
    out <- strsplit(as.character(colnames(valid_x)),'_') # generate group vector
    valid_y <- do.call(rbind, out)[,1] # first column is group information
    valid_x <-apply(valid_x,1, function(x) .random.imp(x)) # impute missing values
    valid_x<- as.data.frame(valid_x)
    
    return(list(promean = SampleMean, # mean protein abundance across all the samples
                mu=GroupMean, # mean protein abundance in each group
                sigma=GroupVar, # standard deviation in each group
                X = valid_x, # protein abundance matrix
                Y = valid_y # vector with group information
                ))
}
    
#' Simulate extended datasets for sample size estimation
#'
#' @param m number of samples to simulate
#' @param mu a matrix of mean abundance in each phenotype group and each protein
#' @param sigma a matrix of variance in each phenotype group and each protein
#' @return \emph{X} A simulated matrix with required sample size
#' @return \emph{Y} Group information corresponding with \emph{X}
#' @keywords internal
.sampleSimulation <- function(m, mu, sigma) {
    
    nproteins <- nrow(mu)
    ngroup <- ncol(mu)
    ## Determine the size of each phenotype group
    samplesize <- .determineSampleSizeinEachGroup(m, ngroup)
    
    ## Simulate the data matrix
    sim_matrix <- matrix(rep(0, nproteins * m), ncol=m)
    for (i in 1:nproteins) {
        abun <- NULL
        for (j in 1:ngroup) {
            abun <- c(abun, rnorm(samplesize[j], mu[i, j], sigma[i, j]))
        }
        sim_matrix[i, ] <- abun
    }
    sim_matrix <- t(sim_matrix)
    colnames(sim_matrix) <- rownames(mu)
    #Simulate the phenotype information
    group <- rep(c(1:ngroup), times=samplesize)
    
    index <- sample(length(group), length(group))
    sim_matrix <- sim_matrix[index, ]
    group <- group[index]
    
    return(list(X=sim_matrix,
                Y=as.factor(group)))
}
    
#' Determine the size of each phenotype group
#'
#' @param m sample size
#' @param ngroup number of phenotype groups
#' @return vector of sample size in each group
#' @keywords internal
.determineSampleSizeinEachGroup <- function(m, ngroup) {
    samplesize <- vector(mode="numeric", length=ngroup)
    counter <- 0
    
    while (counter < m) {
        for (i in 1:ngroup) {
            if (counter < m) {
                counter <- counter + 1
                samplesize[i] <- samplesize[i] + 1
            }
        }
    }
    return(samplesize)
}
    
#' For each protein, impute the missing values based on the observed values
#'
#' @param data protein abundance data for one protein.
#' @return Imputed protein abundance data
#' @keywords internal
.random.imp <- function (data){
    missing <- is.na(data) # count missing values
    n.missing <- sum(missing)
    data.obs <- data[!missing] # keep the observed values
    imputed <- data
    # impute missing values by randomly selecting observed values
    imputed[missing] <- sample (data.obs, n.missing, replace=TRUE)
    return (imputed)
} 
    
#############################################
## designSampleSizeClassificationPlots
#############################################
#' Visualization for sample size calculation in classification problem
#' @param data output from function \code{\link{designSampleSizeClassification}}
#' @description To illustrate the mean classification accuracy under different protein number and sample size. The input is the result from function \code{\link{designSampleSizeClassification}}.
#' @return Plot for sample size estimation. x-axis : sample size, y-axis: mean predictive accuracy. Color: different protein number.
#' @details Data in the example is based on the results of sample size calculation in classification problem from function \code{\link{designSampleSizeClassification}}
#' @author Ting Huang, Meena Choi, Olga Vitek. 
#'
#' Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
#' @references T. Huang et al.  TBD  2018
#' @examples # Consider the training set from a colorectal cancer study
#' # Subjects are from control group or colorectal cancer group
#' # 72 proteins were targeted with SRM
#' require(MSstatsBioData)
#' set.seed(1235)
#' data(SRM_crc_training)
#' QuantCRCSRM <- dataProcess(SRM_crc_training, normalization = FALSE)
#' # estimate the mean predictive accuray under different sizes of training data
#' # n_sample is the number of different sample size to simulate
#' # Datasets with 10 different sample size and 3 different protein numbers are simulated 
#' result.crc.srm <- designSampleSizeClassification(data=QuantCRCSRM, 
#' n_sample = 10, 
#' sample_incr = 10, 
#' protein_desc = 0.33, 
#' iter = 50)
#' designSampleSizeClassificationPlots(data=result.crc.srm)
#' @importFrom reshape2 melt
#' @export
designSampleSizeClassificationPlots <- function(data) {
    
    ## ggplot needs a long format dataframe
    ## get the mean accuracy
    meandata <- as.data.frame(data$meanPA)
    meandata$Protein_number <- rownames(meandata)
    meandata <- reshape2::melt(meandata, id.vars = "Protein_number", variable.name = "Train_size", value.name = "mean")
    
    ## get the variance
    vardata <- as.data.frame(data$varPA)
    vardata$Protein_number <- rownames(vardata)
    vardata <- reshape2::melt(vardata, id.vars = "Protein_number", variable.name = "Train_size", value.name = "var")
    
    ## perform the join
    plotdata <- merge(meandata, vardata, all=TRUE)
    # get standard deviation column
    plotdata$sd <- sqrt(plotdata$var)
    # make sure train size is numeric
    plotdata$Train_size <- gsub("tra", "", plotdata$Train_size)
    plotdata$Train_size <- as.numeric(as.character(plotdata$Train_size))  
    # make sure Protein_number is ordered factor
    plotdata$Protein_number <- gsub("prot", "", plotdata$Protein_number)
    plotdata$Protein_number <- factor(plotdata$Protein_number, levels = sort(as.numeric(unique(plotdata$Protein_number))))
    
    ## make the plot
    p <- ggplot(data = plotdata, aes(x= Train_size, y= mean, group = Protein_number, colour = Protein_number)) +
        geom_point() +
        geom_line() +
        labs(title="Sample size estimation", x="Sample size", y='Mean accuracy') +
        guides(color=guide_legend(title="Protein Number"))
    
    print(p)
}



