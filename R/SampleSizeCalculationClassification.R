#' For classification problem, estimate the optimal number of training samples under different size of validation data
#'
#' param data Dataframe of the run level data reported by MSStats dataProcess() function or dataframe with columns "RUN", "Protein", "LogIntensities", "GROUP", "SUBJECT". Protein intensities must be log2 transformed.
#' param n number of steps, which is the number of different sample sizes to simulate. 
#' param step number of samples per condition to increase at each step. 
#' param noise the fraction of protein variance to change at each step. Positive values indicates increasing noise in the data while negative value means decreasing noise.
#' param iter number of times to repeat simulation experiments
#'
#' @import doMC
#' @importFrom randomForest randomForest combine
#' 
#' @export
designSampleSizeClassification <- function(data, 
                                           n = 5, step = 20, noise = 0.0, iter = 5) {
  
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
    
    ## check n
    if (!(n>0) & (is.integer(n))) {
        processout <- rbind(processout,
                            "The required input - n : was not a positive integer. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set n as a positive integer.")
    }
    
    ## check step
    if (!(step>0) & (is.integer(step))) {
        processout <- rbind(processout,
                            "The required input - step : was not a positive integer. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set step as a positive integer.")
    }
    
    ## check noise
    if (!is.numeric(noise)) {
        processout <- rbind(processout,
                            "The required input - noise : was not a numeric value - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set noise as a numeric value.")
    }
    
    ## check iteration
    if (!(iter>0) & (is.integer(iter))) {
        processout <- rbind(processout,
                            "The required input - iter : was not a positive integer. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)
        
        stop("Please set iter as a positive integer.")
    }
    
    processout <- rbind(processout, paste0("n = ", n))
    processout <- rbind(processout, paste0("step = ", step))
    processout <- rbind(processout, paste0("noise = ", noise))
    processout <- rbind(processout, paste0("iter = ", iter))
    
    ## Estimate the mean abundance and variance of each protein in each phenotype group
    parameters <- .estimateVar(data)
    
    ##
    processout <- rbind(processout, 
                        "Calculated variance component. - okay")
    write.table(processout, file=finalfile, row.names=FALSE)
    
    mu <- parameters$mu
    sigma <- parameters$sigma
    
    ngroup <- length(unique(data$ProcessedData$GROUP_ORIGINAL)) ## Number of phenotype groups
    ## Generate the vector of training sample size
    train_size <- seq.int(from = step, to = step * n, length.out = n)
    train_size <- train_size*ngroup
    
    ## Generate the vector of validation sample size
    validation_size <- seq.int(from = step, to = step*n, length.out = n)
    validation_size <- validation_size * ngroup
    
    message(" Start to run the simulation...")
    PA <- list()
    for (i in 1:iter) { ## Number of iterations
        message("  Iteration: ", i)
        rep <- list()
        new_sigma_1 <- sigma
        accur <- matrix(rep(0, times=length(train_size) * length(validation_size)), nrow=length(train_size))
        
        for (m1 in seq_along(train_size)) {
            message("    Sample size of training data: ", train_size[m1])
            train <- .sampleSimulation(train_size[m1], mu, new_sigma_1) #Simulate training data
            x <- as.data.frame(train$X)
            y <- as.factor(train$Y)
            #Train random forest on training data
            #rf <- randomForest::randomForest(x=x, y=y, ntree = 100) 
            # parallel computing for random forest
            doMC::registerDoMC()
      
            rf <- foreach(ntree=rep(10, 10), .combine=randomForest::combine, .multicombine=TRUE,
                          .packages='randomForest') %dopar% {
                              randomForest(x=x,
                                           y=y,
                                           ntree=ntree)
                              }
            new_sigma_1 <- new_sigma_1 * (noise + 1.0)
            new_sigma_2 <- sigma
            
            for (m2 in seq_along(validation_size)) { ## Different validation sample sizes
                message("      Sample size of validation data: ", validation_size[m2])
                valid <- .sampleSimulation(validation_size[m2], mu, new_sigma_2) #Simulate validation data
                ## Calculate predictive accuracy
                valid_x <- as.data.frame(valid$X)
                valid_y <- as.factor(valid$Y)
                rf.pred <- predict(rf, valid_x) #Predict validation data
                accuracy <- sum(diag(table(rf.pred,valid_y))) / validation_size[m2]
                accur[m1,m2] <- accuracy
                new_sigma_2 <- new_sigma_2 * (noise + 1.0)
            }
        }
        PA[[i]] <- accur
    }
    
    ## Calculate the mean accuracy and variance
    meanPA <- matrix(rep(0, times=length(train_size) * length(validation_size)), nrow=length(train_size))
    varPA <- matrix(rep(0, times=length(train_size) * length(validation_size)), nrow=length(train_size))
    
    for (m1 in seq_along(train_size)) {
        for (m2 in seq_along(validation_size)) {
            temp <- NULL
            for (i in 1:iter) {
                temp <- c(temp, PA[[i]][m1, m2])
            }
            meanPA[m1, m2] <- mean(temp)
            varPA[m1, m2] <- var(temp)
        }
    }
    rownames(meanPA) <- paste0("tra", train_size)
    colnames(meanPA) <- paste0("val", validation_size)
    rownames(varPA) <- paste0("tra", train_size) 
    colnames(varPA) <- paste0("val", validation_size)
    opt <- train_size[max.col(t(meanPA))]
    names(opt) <- paste("val", validation_size, sep="")
    ##
    processout <- rbind(processout, 
                        "The number of sample is calculated. - okay")
    write.table(processout, file=finalfile, row.names=FALSE)
    
    return(list(optTrainSize = opt, 
                meanPA = meanPA, 
                varPA = varPA))
}
    
    

#' Estimate the mean abundance and variance of each protein in each phenotype group
#'
#' param data The run level data reported by dataProcess function.
#' return mu is the mean abundance matrix of each protein in each phenotype group; sigma is the sd matrix of each protein in each phenotype group.
#' keywords internal
.estimateVar <- function(data) {
    
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
    Var <- matrix(rep(0.0, length(res) * length(groups)), ncol = length(groups))
    MeanAbun <- matrix(rep(0.0, length(res) * length(groups)), ncol = length(groups))
    count <- 0
    
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
            count <- count + 1
            abun <- summary(fit.full)$coefficients[, 1]
            abun[-1] <- abun[1] + abun[-1]
            MeanAbun[count, ] <- abun
            Var[count, ] <- rep(sqrt(sum(VarComponent[i, ], na.rm = TRUE)), times=length(abun))
        }  
    } ## end-loop
    
    return(list(mu=MeanAbun[1:count, ], 
                sigma=Var[1:count, ]))
}
    
    

#' Simulate extended datasets for sample size estimation
#'
#' param m number of samples to simulate
#' param mu a matrix of mean abundance in each phenotype group and each protein
#' param sigma a matrix of variance in each phenotype group and each protein
#' 
#' return A simulated matrix with required sample size
#' keywords internal
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
#' param m sample size
#' param ngroup number of phenotype groups
#'
#' return vector of sample size in each group
#' keywords internal
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
    
    
#############################################
## designSampleSizeClassificationPlots
#############################################
#' @importFrom reshape2 melt
#' @export
designSampleSizeClassificationPlots <- function(data) {
    
    ## ggplot needs a long format dataframe
    ## get the mean accuracy
    meandata <- as.data.frame(data$meanPA)
    meandata$Training_size <- rownames(meandata)
    meandata <- reshape2::melt(meandata, id.vars = "Training_size", variable.name = "Valid_size", value.name = "mean")
    
    ## get the variance
    vardata <- as.data.frame(data$varPA)
    vardata$Training_size <- rownames(vardata)
    vardata <- reshape2::melt(vardata, id.vars = "Training_size", variable.name = "Valid_size", value.name = "var")
    
    ## perform the join
    plotdata <- merge(meandata, vardata, all=TRUE)
    
    plotdata$Training_size <- gsub("tra", "", plotdata$Training_size)
    vardata$Training_size <- as.numeric(as.character(plotdata$Training_size))
    plotdata$Valid_size <- gsub("val", "", plotdata$Valid_size)
    plotdata$Valid_size <- as.character(plotdata$Valid_size)
    plotdata$Valid_size <- factor(plotdata$Valid_size, levels=sort(as.numeric(unique(plotdata$Valid_size))))
    
    ## make the plot
    plotdata$sd <- sqrt(plotdata$var)
    p <- ggplot(data = plotdata , 
                aes(x= as.numeric(as.character(Training_size)), y= mean, group = Valid_size, colour = Valid_size)) +
        geom_point() +
        geom_line() +
        labs(title="Sample size estimation", x="Size of training data", y='Mean accuracy') +
        guides(color=guide_legend(title="Size of validation data"))
    
    print(p)
}



