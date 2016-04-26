###########Function of FeatureSelection###########
#Input is the data frame 'work' after conducting normalization
#Output is the data frame 'work' after removing noisy features 
#This process removes features with interference across the MS runs
#remove_proteins_with_interference==TRUE allows the algorithm to remove the whole protein if deem interfered

.feature_selection <- function(work, 
                               remove_proteins_with_interference, 
                               max.iter=10, 
                               Improve_Margin_DIA=0.15, 
                               Step2='Subplot'){


  #Label-free experiments (SWATH)
	if(nlevels(work$LABEL) == 1 & length(unique(work$TRANSITION)) > 1){ #1.2

    #Set an arbitrary value to start the while loop
    Cut <- '60%'
    Iteration <- 0
    Cut.History <- vector()
    work.0 <- work

    while(Cut != '100%' & Iteration < max.iter){ #While.SWATH
      #We can merely work on one peptide at a time
      N.Prot <- nlevels(work$PROTEIN) 
      
      if(N.Prot > 200) {
        Max.Iter <- 5
      }
       
      #Create a data frame storing the score of interference for each fragment
      Interference.Score <- data.frame(Protein=vector(), Peptide=vector(), Feature=vector(), Interference=vector())
      Index.DIA <- 0

      #loop over each protein
      for(i in 1:N.Prot){ #F.2

        sub1 <- subset(work, PROTEIN==unique(PROTEIN)[i])
        N.Pep <- length(unique(sub1$PEPTIDE))
          
        ## show progress
        message(paste("Selection features or peptides for protein ", unique(sub1$PROTEIN), "(", i, " of ", N.Prot, ")"))

        #Compute the score of interference at each fragment of each peptide
        for(j in 1:N.Pep){ #F.3 

          sub2 <- subset(sub1, PEPTIDE==unique(PEPTIDE)[j])
          sub2$FEATURE <- factor(sub2$FEATURE, levels=unique(sub2$FEATURE))
          NF <- length(unique(sub2$FEATURE))	

          #First, remove the bottom 30% fragments with the lowest average abundance of each peptide
          #Do so only when there are at least 8 MS2 features belong this peptide
          if( NF >= 8 ){
            M1 <- tapply(sub2$ABUNDANCE, sub2$FEATURE, function(x) mean(x, na.rm=TRUE)) 
            Top70.Cut <- round(length(M1) * 0.7) 
            Keep.Top70 <- names(M1[rank(-M1) <= Top70.Cut])

            #Here we only work on the top 70% most abundant fragments in each peptide
            sub2.2 <- subset(sub2, FEATURE %in% Keep.Top70)
          } else {
            sub2.2 <- sub2
          }

          sub2.2$FEATURE <- factor(sub2.2$FEATURE, levels=unique(sub2.2$FEATURE))
          N.Feature <- nlevels(sub2.2$FEATURE) 
          sub2.2$RUN <- factor(sub2.2$RUN, levels=unique(sub2.2$RUN))

          #Use TMP to do robust run quantification
          data_tmp = dcast(RUN ~ FEATURE, data=sub2.2, value.var='ABUNDANCE', keep=TRUE)
          rownames(data_tmp) <- data_tmp$RUN
          data_tmp <- data_tmp[, -1]
          data_tmp[data_tmp == 1] <- NA

          TMP <- medpolish(data_tmp, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
          TMP.Run <- TMP$overall + TMP$row           

          #Calculate the score of interference on each fragment
          #In fact, var(data_tmp[,k]-TMP.Run)=var(TMP$res[,k])
          for(k in 1:N.Feature){ #F.4
             
            Interference <- var(TMP$res[,k], na.rm=TRUE)
            Interference.Score[(Index.DIA+1), 'Protein'] <- as.character(unique(sub2.2$PROTEIN)) 
            Interference.Score[(Index.DIA+1), 'Peptide'] <- as.character(unique(sub2.2$PEPTIDE)) 
            Interference.Score[(Index.DIA+1), 'Feature'] <- colnames(data_tmp)[k]
            Interference.Score[(Index.DIA+1), 'Interference'] <- Interference
                  
            Index.DIA <- Index.DIA+1
          } #F.4
          #End of loop for quantifying the amount of interference at each feature at this peptide

        } #F.3
        #End of loop for the peptide

      } #F.2
      #End of loop over each protein

      #Overall, the percentage of most interfered features are removed was determined by where the decrease of the interference stabalized
      Interference.Score$Protein_Feature <- paste(Interference.Score$Protein, Interference.Score$Feature, sep='_')
	    QQ <- quantile(Interference.Score$Interference, probs=seq(0, 1, 0.05))
      QQ2 <- diff(QQ); QQ3 <- QQ2/QQ[-1]
      
	    #Stop at where the improvement is less than the improve margin. Remove 60% features with most interference if this cutoff cannot be determined by this criteria
	    suppressWarnings(Pos <- max(which(QQ3 < Improve_Margin_DIA)))
	    Cut <- '60%'; if(Pos != -Inf){Cut <- names(QQ3[Pos])}; if(Pos < 12){Cut <- '60%'}
	    Cut.Noise <- QQ[Cut]	              
      Interfere <- Interference.Score[Interference.Score$Interference < Cut.Noise,]

      #Whichever peptides only has one feature left should be removed as it indicates this whole peptide has interference
      Interfere$Protein_Peptide <- paste(Interfere$Protein, Interfere$Peptide, sep='_')
      Table <- table(Interfere$Protein_Peptide)
      Drop <- names(Table[Table == 1])
      Interfere <- subset(Interfere, !(Protein_Peptide %in% Drop))

      #Now use the remaining features to compute the average amount of interference on each peptide
      #Remove the most interfered peptide. We prepare the case where not all of the peptides are proteotypic
      Keep <- unique(Interfere$Protein_Feature)
      work$Protein_Feature <- paste(work$PROTEIN, work$FEATURE, sep='_')
      work <- subset(work, Protein_Feature %in% Keep)
      Iteration <- Iteration + 1
      Cut.History[Iteration] <- Cut
      if(Iteration == 3){
        if(sum(Cut.History == rep('60%', 3)) == 3) {
          break
        }
      }
      
	    print(paste('Iteration', Iteration, 'kept', Cut, ' of the features', sep=' '))	 
	    	        
    } #While.SWATH

    #Summarize the profile for each peptide to avoid more influence of peptides with more MS2 peaks
    work$Protein_Peptide <- paste(work$PROTEIN, work$PEPTIDE, sep="_")
    N.Peptide <- length(unique(work$Protein_Peptide))
	  	
		Profile.Pep <- data.frame(Protein=vector(), Peptide=vector(), RUN=vector(), Run.Intensity=vector())
		Info <- unique(work[,c("RUN", "GROUP", "SUBJECT", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "SUBJECT_NESTED")])
		row.names(Info) <- NULL; Index.Pep <- 0

		for(k in 1:N.Peptide){#F.9
      Temp <- subset(work, Protein_Peptide == unique(Protein_Peptide)[k])

			Temp$RUN <- factor(Temp$RUN)
      Temp$FEATURE <- factor(Temp$FEATURE)
      data_TMP = dcast(RUN ~ FEATURE, data=Temp, value.var='ABUNDANCE', keep=TRUE)
  		rownames(data_TMP) <- data_TMP$RUN
  	  data_TMP <- data_TMP[,-1]

      TMP <- medpolish(data_TMP, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
      TMP.Run <- TMP$overall + TMP$row         
			N <- length(TMP.Run)
      Prot <- unique(Temp$PROTEIN)
      Pep <- unique(Temp$PEPTIDE)

			Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'Protein'] <- rep(as.character(Prot), N)
			Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'Peptide'] <- rep(as.character(Pep), N)
			Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'RUN'] <- as.numeric(names(TMP.Run))
			Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'Run.Intensity'] <- TMP.Run

			Index.Pep <- Index.Pep + N
 
    }#F.9  End of loop for peptides

		Profile.Pep <- merge(Profile.Pep, Info, by='RUN')
		Profile.Pep$Protein_Peptide <- paste(Profile.Pep$Protein, Profile.Pep$Peptide, sep='_')

		#Remove the peptides which have different profiles than the others
    if(toupper(Step2) == 'SUBPLOT'){ #F9.2
      N.Prot <- length(unique(Profile.Pep$Protein))
      Keep.Pep <- vector()
      Index <- 0
      
      for(i in 1:N.Prot){ #F.10

        Out <- data.frame(Protein=vector(), Peptide=vector(), Error=vector())
				sub <- subset(Profile.Pep, Protein == unique(Protein)[i])
				sub$RUN <- factor(sub$RUN); sub$Peptide <- factor(sub$Peptide)
        sub_tmp = dcast(RUN ~ Peptide, data=sub, value.var='Run.Intensity', keep=TRUE)
        rownames(sub_tmp) <- sub_tmp$RUN
        sub_tmp <- sub_tmp[,-1]

        tmp <- medpolish(sub_tmp, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
        tmp.run <- tmp$overall + tmp$row   
				N <- dim(sub_tmp)[2]

				#Record the variance components of each peptide
				Out[1:N, 'Protein'] <- rep(as.character(unique(sub$Protein)),N)
				Out[1:N, 'Peptide'] <- as.character(colnames(sub_tmp))

				for(j in 1:N){

					Out[j, 'Error'] <- var(sub_tmp[, j] - tmp.run, na.rm=TRUE) 
				
				}

				Out$Protein_Peptide <- paste(Out$Protein, Out$Peptide, sep='_')

				R <- Out$Error[order(Out$Error)]
				R3 <- diff(R)/R[-length(R)]

				#The cutoff is determined while the increase of error for the next peptide is more than 25%
				#Only selecting one peptide in SWATH could be dangerous and unreasonable, so at least select 2 peptide.
				Pos <- suppressWarnings(min(which(R3 > 0.20)))
				if( Pos == Inf ) {
          Pos <- 2
				}
				if( Pos == 1 ) {
          Pos <- which(R3 > 0.20)[2]				
				}
        
				Cut.Pep <- R[Pos]	
				Keep.Pep[(Index+1):(Index+Pos)] <- Out[Out$Error <= Cut.Pep, 'Protein_Peptide'] 
				Protein <- as.character(unique(sub$Protein))
				print(paste('There are', Pos, 'peptide(s) being kept in Protein', Protein, '(', i, 'out of', N.Prot, ')', sep=' '))	 

				Index <- Index+Pos

      } #F.10 (End of loop for protein)

    } 
    
		if( toupper(Step2) == 'WHOLEPLOT' ){

      N.Prot <- length(unique(Profile.Pep$Protein))
      repeated <- .checkRepeated(work)
      Out2 <- data.frame(Protein=vector(), Peptide=vector(), Error.Port=vector(), Error.Pep=vector())				
      Index.2 <- 0				

      for( k in 1:N.Prot ){ #F.11
				
        sub <- subset(Profile.Pep, Protein == unique(Protein)[k])
    
				sub$GROUP <- factor(sub$GROUP)
        sub$SUBJECT <- factor(sub$SUBJECT)
        sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
        sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)
        sub$SUBJECT_NESTED <- factor(sub$SUBJECT_NESTED)
		    sub$RUN <- factor(sub$RUN)
    
        singleSubject <- .checkSingleSubject(sub)
        TechReplicate <- .checkTechReplicate(sub) ## use for label-free model
    
				## case-control
        if (!repeated) {
          if (!TechReplicate | singleSubject) {
            fit.full <- lm(Run.Intensity ~ GROUP , data = sub)
          } else {
            fit.full <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = sub)
          }
        } else { ## time-course
          if (singleSubject) {
            fit.full <- lm(Run.Intensity ~ GROUP , data = sub)
      		} else { ## no single subject
            if (!TechReplicate) {
              fit.full <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = sub)
            } else {
              fit.full <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT), data = sub) ## SUBJECT==SUBJECT_NESTED here
            }
          }	
        } ## time-course

        N.Pep <- length(unique(sub$Peptide))
        sub$Res <- summary(fit.full)$res
        
        for(j in 1:N.Pep){ #F.12

          Out2[(Index.2+1), 'Protein'] <- as.character(unique(sub$Protein))
          Out2[(Index.2+1), 'Peptide'] <- as.character(unique(sub$Peptide)[j])
          Residual <- sub[sub$Peptide == unique(sub$Peptide)[j], 'Res']
          Out2[(Index.2+1), 'Error.Prot'] <- var(Residual, na.rm=TRUE)

          subsub <- subset(sub, Peptide == unique(Peptide)[j])

          ## case-control
          if (!repeated) {
            if (!TechReplicate | singleSubject) {
              fit.sub <- lm(Run.Intensity ~ GROUP , data = subsub)
            } else {
              fit.sub <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = subsub)
            }
          } else { ## time-course
            if (singleSubject) {
              fit.sub <- lm(Run.Intensity ~ GROUP , data = subsub)
            } else { ## no single subject
              if (!TechReplicate) {
                fit.sub <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = subsub)
              } else {
                fit.sub <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT), data = subsub) ## SUBJECT==SUBJECT_NESTED here
              }
            }	
          } ## time-course

          Out2[(Index.2+1), 'Error.Pep'] <- var(summary(fit.sub)$res, na.rm=TRUE)
          Index.2 <- Index.2+1

        } #F.12 End of loop for peptides

				#Now conduct the selection process for each protein

      }#F.11 End of loop for proteins


    } #F9.2
 	
		work.3 <- subset(work, Protein_Peptide %in% Keep.Pep)

    ##If we do not allow any of the proteins being totally removed, we keep the most abundant peptide if it is totally removed
    if( !remove_proteins_with_interference ){ #F.7
			
      Prot2 <- unique(work.3$PROTEIN)
      Prot1 <- unique(work.0$PROTEIN)
      Loss <- setdiff(Prot1, Prot2)
      work.0$Protein_Peptide <- paste(work.0$PROTEIN, work.0$PEPTIDE, sep='.') 
      Keep.Peptide <- vector()
      IK <- 0
				
      #Recover the proteins which are totally removed
      if( length(Loss) > 0 ){
        sub.Loss <- subset(work.0, PROTEIN %in% Loss)
        N2.Prot <- length(unique(sub.Loss$PROTEIN))
        
        for(j in 1:N2.Prot){

          subsub.Loss <- subset(sub.Loss, PROTEIN == unique(sub.Loss$PROTEIN)[j])
					N.Pep2 <- length(unique(subsub.Loss$PEPTIDE))

					#Get the most abundant peptide of this protein
					AVG.Pep <- tapply(subsub.Loss$ABUNDANCE, subsub.Loss$Protein_Peptide, function(x) mean(x, na.rm=TRUE))
					Top1 <- names(AVG.Pep[rank(-AVG.Pep) == 1])
					Keep.Peptide[IK+1] <- Top1
          IK <- IK+1
				}

				work.keep <- subset(work.0, Protein_Peptide %in% Keep.Peptide)
				work.4 <- rbind(work.3, work.keep)
				work.4$Protein_Feature <- NULL
        work.4$Protein_Peptide <- NULL
				work <- work.4
      } else {
				#The case of no proteins was totally removed
				work.3$Protein_Feature <- NULL
        work.3$Protein_Peptide <- NULL
				work <- work.3
			}
		} else {

      work.3$Protein_Feature <- NULL; work.3$Protein_Peptide <- NULL

      #data frame 'work.3' contain the selected features and 'work' is the data frame to be returned 
      work <- work.3

    } #F.7   
	} # 1.2 : End of label-free SWATH experiment


	#Label-free experiments (Shotgun; DDA)
	if(nlevels(work$LABEL) == 1 & nlevels(work$TRANSITION) == 1){ #1.2.2

    # In DDA (shotgun) experiment, there is no need to filter out bad fragments in each peptide
    work.2 <- work		

    #Remove the peptides with large unexplained variation within same protein roughly have same profile
		#Iteration <- 0; Cut.History <- vector(); Cut <- '60%'
    #if (ByProtein==TRUE) max.iter <- 1

    #while(Cut != '100%' & Iteration < max.iter){ #DDA_While

    # Create the data frame storing the selected features of each Protein
    #Out <- data.frame(Protein=vector(), Feature=vector(), Model.Based.Error=vector())
    #Index.FS <- 0
    N.Prot <- length(unique(work.2$PROTEIN))	
    Keep.Feature <- vector()
    index.KF <- 0

    ## Calculate the variance component estimates of each feature (peptide) for each protein
    for(i in 1:N.Prot){ #DDA_1.3

      Out <- data.frame(Protein=vector(), Feature=vector(), Model.Based.Error=vector())
      Index.FS <- 0
      DDA.1 <- subset(work.2, PROTEIN == unique(PROTEIN)[i])
      N.Feature <- length(unique(DDA.1$FEATURE))
      DDA.1$RUN <- factor(DDA.1$RUN, levels=unique(DDA.1$RUN))
      DDA.1$FEATURE <- factor(DDA.1$FEATURE, levels=unique(DDA.1$FEATURE))

      ## show progress
			message(paste("Selection features or peptides for protein ", unique(DDA.1$PROTEIN), "(", i, " of ", N.Prot, ")"))
				
			## Assess the noise based on model for each protein, assuming every peptide has same profile within each protein##
		  ## The peptides with different profiles have large variance component ##
             
			#Use TMP to do robust run quantification
      data_tmp = dcast(RUN ~ FEATURE, data=DDA.1, value.var='ABUNDANCE', keep=TRUE)
  		rownames(data_tmp) <- data_tmp$RUN
  	  data_tmp <- data_tmp[,-1]

  		TMP <- medpolish(data_tmp, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
	    TMP.Run <- TMP$overall + TMP$row           

			Error <- vector()
			#######Get the per-peptide variance and store it#######
			if( N.Feature == 1 ){
        #Keep the protein with only one peptide quantified although thr variance component cannot be assessed
				Error <- 0

				Out[(Index.FS+1), 'Protein'] <- as.character(unique(DDA.1$PROTEIN))
				Out[(Index.FS+1), 'Feature'] <- as.character(unique(DDA.1$FEATURE))
				Out[(Index.FS+1), 'Model.Based.Error'] <- Error
				#Index.FS <- Index.FS+N.Feature   
				Out$Protein_Feature <- paste(Out$Protein, Out$Feature, sep='_')
				Keep.Feature[(index.KF+1)] <- as.character(Out$Protein_Feature)				
				index.KF <- index.KF + 1
				Pos <- 1
      } else {

        for(j in 1:N.Feature){ #DDA_1.3.1

          Error[j] <- var((data_tmp[,j]-TMP.Run), na.rm=TRUE)

				} #DDA_1.3.1

				## Pull out the results to the designated data frame
				Out[(Index.FS+1):(Index.FS+N.Feature), 'Protein'] <- rep(as.character(unique(DDA.1$PROTEIN)), N.Feature)
				Out[(Index.FS+1):(Index.FS+N.Feature), 'Feature'] <- names(data_tmp)
				Out[(Index.FS+1):(Index.FS+N.Feature), 'Model.Based.Error'] <- Error
				Out$Protein_Feature <- paste(Out$Protein, Out$Feature, sep='_')

				#Index.FS <- Index.FS+N.Feature   
				R <- Out$Model.Based.Error[order(Out$Model.Based.Error)]
				R <- R[!is.na(R)]
				R2 <- diff(R); R3 <- R2/R[-length(R)]
				#Hope to keep at least 2 peptides in each protein
				K <- suppressWarnings(min(which(R3 > 0.2), na.rm=TRUE))
        Pos <- 2
				Big <- which(R3 > 0.2)
				ifelse(K == 1 & N.Feature != 2 & length(Big) > 1, Pos <- Big[2], Pos <- K)
        
				if(K == 1 & length(Big) == 1) {
          Pos <- N.Feature
				}
				
        if(N.Feature == 2 & K == Inf) {
          Pos <- 2
        }
					
        if(K == Inf) {
          Pos <- 2
        }

        Cut <- R[Pos]
        Keep.Feature[(index.KF+1):(index.KF+Pos)] <-  Out[Out$Model.Based.Error <= Cut, 'Protein_Feature']
        index.KF <- index.KF + Pos
			}

			message(paste(Pos, "peptides were selected in protein", unique(DDA.1$PROTEIN), "(", i, " of ", N.Prot, ")"))
	                                              
		} # DDA_1.3 (End of loop for proteins)

		work$Protein_Feature <- paste(work$PROTEIN, work$FEATURE, sep='_')			
		work <- subset(work, Protein_Feature %in% Keep.Feature)
		work$Protein_Feature <- NULL

	} #1.2.2  : End of the DDA experiment 

	## Label-based experiments
  if( nlevels(work$LABEL) == 2 ){ #2.4

    N.Prot <- length(unique((work$PROTEIN)))
		work.0 <- work

    #Create a data frame storing the output
    Out <- data.frame(Protein=vector(), Peptide=vector(), Feature=vector(), Label=vector(), Interference.Score=vector())   
    Index <- 0
    		
    #Loop over proteins: One protein at a time
    for(i in 1:N.Prot){ #2.4.1
         
      Sub <- subset(work, PROTEIN == unique(PROTEIN)[i])
      N.Pep <- length(unique(Sub$PEPTIDE))
       			
      ## show progress
			message(paste("Selection features or peptides for protein ", unique(Sub$PROTEIN), "(", i, " of ", N.Prot, ")"))
				
      #Loop for the peptides in this protein
      for(j in 1:N.Pep){ #2.4.2

        Sub2 <- subset(Sub, PEPTIDE == unique(PEPTIDE)[j])    
        
        #Heavy peptide
        Sub2.H <- subset(Sub2, LABEL=='H')

        Sub2.H$RUN <- factor(Sub2.H$RUN, levels=unique(Sub2.H$RUN))
        Sub2.H$FEATURE <- factor(Sub2.H$FEATURE, levels=unique(Sub2.H$FEATURE))
            
        #Apply TMP to estimate the run effect
        data_w = dcast(RUN ~ FEATURE, data=Sub2.H, value.var='ABUNDANCE', keep=TRUE)
        rownames(data_w) <- data_w$RUN
        data_w <- data_w[,-1]
        data_w[data_w<=1] <- NA

  			TMP <- medpolish(data_w, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
        TMP.Run <- TMP$overall + TMP$row

        N.Feature <- length(unique(Sub2.H$FEATURE))

        #Calculate the variance of the residuals for each feature
        for( k in 1:N.Feature ){ #2.4.3

          #If there is only one feature in this peptide, no reproducibility can be assessed
					if( N.Feature == 1 ){
              #No need to remove this peptide at this moment
							Error <- 0
          } else {    
							Res <- data_w[,k] - TMP.Run
							Error <- var(Res, na.rm=TRUE)
          }
              
					Out[(Index+1), 'Protein'] <- as.character(unique(Sub2.H$PROTEIN))
					Out[(Index+1), 'Peptide'] <- as.character(unique(Sub2.H$PEPTIDE))
					Out[(Index+1), 'Feature'] <- as.character(unique(Sub2.H$FEATURE)[k])
					Out[(Index+1), 'Label'] <- 'H'
					Out[(Index+1), 'Interference.Score'] <- Error

					Index <- Index+1
				} #2.4.3 : End of loop for calculating the interference score for each feature in this peptide at heavy

        #Light peptide
        Sub2.L <- subset(Sub2, LABEL=='L')

        Sub2.L$RUN <- factor(Sub2.L$RUN, levels=unique(Sub2.L$RUN))
        Sub2.L$FEATURE <- factor(Sub2.L$FEATURE, levels=unique(Sub2.L$FEATURE))
        Sub2.L <- subset(Sub2.L, GROUP!=0)  

        #Apply TMP to estimate the run effect
        data_w.L = dcast(RUN ~ FEATURE, data=Sub2.L, value.var='ABUNDANCE', keep=TRUE)
				rownames(data_w.L) <- data_w.L$RUN
				data_w.L <- data_w.L[,-1]
				data_w.L[data_w.L<=1]<-NA

				TMP.L <- medpolish(data_w.L, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
				TMP.Run.L <- TMP.L$overall + TMP.L$row

				N.Feature <- length(unique(Sub2.L$FEATURE))

				#Calculate the variance of the residuals for each feature on Light
				for(k in 1:N.Feature){ #2.4.3.2

          #If there is only one feature in this peptide, no reproducibility can be assess. The removals of this single features will be decided by the consensus later.
					if ( N.Feature == 1 ) {
					  Error <- 0
					} else {    
							Res <- data_w.L[,k] - TMP.Run.L
							Error.L <- var(Res, na.rm=TRUE)
          }
                   
					Out[(Index+1), 'Protein'] <- as.character(unique(Sub2.L$PROTEIN))
					Out[(Index+1), 'Peptide'] <- as.character(unique(Sub2.L$PEPTIDE))
					Out[(Index+1), 'Feature'] <- as.character(unique(Sub2.L$FEATURE)[k])
					Out[(Index+1), 'Label'] <- 'L'
					Out[(Index+1), 'Interference.Score'] <- Error.L

					Index <- Index+1
        } #2.4.3.2 : End of loop for calculating the interference score for each feature at this peptide on Light     
      } #2.4.2 : End of loop for the peptides
    } #2.4.1 :End of loop over proteins

    #Preparing for the case of shared peptide, we use the column Protein_Feature
    Out$Prot.F <- paste(Out$Protein, Out$Feature, sep='_') 
       
    #The heuristic approach now is to keep searching the cutoff until the improvement is less than 25%
    Out.L <- subset(Out, Label=='L'); Out.H <- subset(Out, Label=='H')
    Q <- quantile(Out.L$Interference.Score, probs = seq(0, 1, 0.05), na.rm = TRUE)
    Change <- diff(Q)/Q[-1]
    suppressWarnings(P <- max(which(Change<0.25), na.rm=TRUE))
    Cut <- '60%'
    if(P != -Inf){
      Cut <- names(Change[P])
    } 
    if(P < 10){
      Cut <- '50%'
    }
    L.Cut <- Q[Cut]

    #For the internal standard, remove no more than 20% of the features as they are supposed to be of good quality
    Q.H <- quantile(Out.H$Interference.Score, probs = seq(0, 1, 0.05), na.rm = TRUE)
    Change.H <- diff(Q.H)/Q.H[-1];
    suppressWarnings(P.H <- max(which(Change.H<0.25), na.rm=TRUE))
    Cut.H <- '80%'
    if( P.H != -Inf ){
      Cut.H <- names(Change.H[P.H])
    }
    if( P.H < 16 ){
      Cut.H <- '80%'
    }
    H.Cut <- Q.H[Cut.H]
            	      
    ##Remove features with interference, which is the first step of the removal.      
    Out.L$Flag.Repro <- 'OK'
    Out.H$Flag.Repro <- 'OK'

    Out.L[is.na(Out.L$Interference.Score), 'Interference.Score'] <- 10000 
    Out.H[is.na(Out.H$Interference.Score), 'Interference.Score'] <- 10000 

    Out.L[Out.L$Interference.Score >= L.Cut, 'Flag.Repro'] <- 'Noisy' 
    Out.H[Out.H$Interference.Score >= H.Cut, 'Flag.Repro'] <- 'Noisy' 

    ##Extract the remaining features before proceeding to the next step
    Keep.L <- Out.L[Out.L$Interference.Score < L.Cut, 'Prot.F']
    Keep.H <- Out.H[Out.H$Interference.Score < H.Cut, 'Prot.F']
       
    #Only keep the features when both of heavy and light are retained
    Keep <- intersect(Keep.L, Keep.H)
    work$Protein_Feature <- paste(work$PROTEIN, work$FEATURE, sep='_')

    #Keep a column of record. Perhaps it may be of interest to see what the reason is the features were removed, reproducibility or consistency.
    work$Filter.Repro <- NA
    work[work$Protein_Feature %in% Keep, 'Filter.Repro'] <- 'Keep'
    work[!(work$Protein_Feature %in% Keep), 'Filter.Repro'] <- 'Flagged'
                   
    #Remove the features due to the reproducibility
    work.keep <- subset(work, Protein_Feature %in% Keep)

    #Also remove the peptide where there is only one transition remained
		#When every transitions but one is interfered, it is more often than not the whole peptide is interfered.
    #Do not assume all the peptides are proteotypic
    #Only do so when a lot of proteins were quantified
			
    Total.Protein <- length(unique(work$PROTEIN))     
		#if(Total.Feature > 50){             
    work.keep$Protein.Peptide <- paste(work.keep$PROTEIN, work.keep$PEPTIDE, sep='.')                  
    Freq <- tapply(work.keep$FEATURE, work.keep$Protein.Peptide, function(x) length(unique(x)))
    OneFeature <- names(Freq[Freq == 1])	
    work.keep <- subset(work.keep, !(Protein.Peptide %in% OneFeature))
    #}

    N.Peptide <- length(unique(work.keep$Protein.Peptide))

    #Normalize the endogenous by equalizing the internal standard before assessing the unique profile per protein
    message(paste("Normalize the endogenous by equalizing the internal standards before assessing the unique profile per protein"))

		for(k in 1:N.Peptide){ #SRM.3

      Peptide.Name <- unique(work.keep$Protein.Peptide)[k]
      sub <- subset(work.keep, Protein.Peptide == Peptide.Name)
      Median.H <- median(sub[sub$LABEL == 'H', 'ABUNDANCE'], na.rm=TRUE)
				
      N.Run <- length(unique(work.keep$RUN))
			for(r in 1:N.Run){

			  Diff <- -median(sub[sub$RUN==r & sub$LABEL=='H', "ABUNDANCE"], na.rm=TRUE) + Median.H
        sub[sub$RUN==r, "ABUNDANCE"] <- sub[sub$RUN==r, "ABUNDANCE"] + Diff

      }

      work.keep[work.keep$Protein.Peptide==Peptide.Name, 'ABUNDANCE'] <- sub$ABUNDANCE
      message(paste(k, "Out of", N.Peptide, "peptides have been normalized."))

    }#SRM.3
		work.keep$Protein.Peptide <- NULL
			
		#Summarize the profile for each peptide to avoid more influence of peptides with more MS2 peaks
		#As endogenous has been normalized, looking at endogenous only is sufficient
		work <- subset(work.keep, LABEL=='L')
    work$Protein_Peptide <- paste(work$PROTEIN, work$PEPTIDE, sep="_")
    N.Peptide <- length(unique(work$Protein_Peptide))
	  	
    Profile.Pep <- data.frame(Protein=vector(), Peptide=vector(), RUN=vector(), Run.Intensity=vector())
    Info <- unique(work[,c("RUN", "GROUP", "SUBJECT", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "SUBJECT_NESTED")])
    row.names(Info) <- NULL
    Index.Pep <- 0

    for(k in 1:N.Peptide){#F.9
      Temp <- subset(work, Protein_Peptide == unique(Protein_Peptide)[k])

      Temp$RUN <- factor(Temp$RUN); Temp$FEATURE <- factor(Temp$FEATURE)
      data_TMP = dcast(RUN ~ FEATURE, data=Temp, value.var='ABUNDANCE', keep=TRUE)
      rownames(data_TMP) <- data_TMP$RUN
      data_TMP <- data_TMP[, -1]

      TMP <- medpolish(data_TMP, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
      TMP.Run <- TMP$overall + TMP$row         
			N <- length(TMP.Run)
      Prot <- unique(Temp$PROTEIN)
      Pep <- unique(Temp$PEPTIDE)

      Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'Protein'] <- rep(as.character(Prot), N)
      Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'Peptide'] <- rep(as.character(Pep), N)
      Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'RUN'] <- as.numeric(names(TMP.Run))
      Profile.Pep[(Index.Pep+1):(Index.Pep+N), 'Run.Intensity'] <- TMP.Run

      Index.Pep <- Index.Pep + N
 
		}#F.9  End of loop for peptides

		Profile.Pep <- merge(Profile.Pep, Info, by='RUN')
		Profile.Pep$Protein_Peptide <- paste(Profile.Pep$Protein, Profile.Pep$Peptide, sep='_')

		#Remove the peptides which have different profiles than the others
		if(toupper(Step2) == 'SUBPLOT'){ #F9.2

			N.Prot <- length(unique(Profile.Pep$Protein))
			Keep.Pep <- vector(); Index <- 0

			for(i in 1:N.Prot){ #F.10

				Out <- data.frame(Protein=vector(), Peptide=vector(), Error=vector(), Avg.Int=vector())
				sub <- subset(Profile.Pep, Protein==unique(Protein)[i])
				sub$RUN <- factor(sub$RUN); sub$Peptide <- factor(sub$Peptide)
        sub_tmp = dcast(RUN ~ Peptide, data=sub, value.var='Run.Intensity', keep=TRUE)
  			rownames(sub_tmp) <- sub_tmp$RUN
   			N <- dim(sub_tmp)[2] - 1
        sub_tmp <- sub_tmp[,-1]

        tmp <- medpolish(sub_tmp, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
        tmp.run <- tmp$overall + tmp$row   

				#Record the variance components of each peptide
				Out[1:N, 'Protein'] <- rep(as.character(unique(sub$Protein)),N)
				if(N > 1) {
          Out[1:N, 'Peptide'] <- as.character(colnames(sub_tmp))
				}
        
				if(N == 1) {
          Out[1:N, 'Peptide'] <- as.character(unique(sub$Peptide))
				}

				for(j in 1:N){
					if(N > 1) {
            Out[j, 'Error'] <- var(sub_tmp[,j]-tmp.run, na.rm=TRUE) 
					}
          
					if(N > 1) {
            Out[j, 'Avg.Int'] <- mean(sub_tmp[,j], na.rm=TRUE)
					}
          
					if(N == 1) {
            Out[j, 'Error'] <- 0; Out[j, 'Avg.Int'] <- 10
					}				
				}

				Out$Protein_Peptide <- paste(Out$Protein, Out$Peptide, sep='_')

				if (N == 1) {
        			
					Keep.Pep[(Index+1)] <- Out$Protein_Peptide
          Pos <- 1
					Index <- Index+1 

				} else if( N == 2 ) {

					#If there are only 2 peptides, they will have same variance component so take the high abundence
					Out <- Out[order(-Out$Avg.Int), ]
          Pos <- 1
					Keep.Pep[(Index+1)] <- Out$Protein_Peptide[1]
					Index <- Index+1 

				} else {

					R <- Out$Error[order(Out$Error)]
          R <- R[!is.na(R)]
					R3 <- diff(R)/R[-length(R)]

					#The cutoff is determined while the increase of error for the next peptide is more than 20%
					Pos <- suppressWarnings(min(which(R3 > 0.20)))
					if(Pos == Inf) Pos <- N

					Cut.Pep <- R[Pos]	
					Keep.Pep[(Index+1):(Index+Pos)] <- Out[Out$Error <= Cut.Pep, 'Protein_Peptide'] 
					Index <- Index+Pos
				}

				Protein <- as.character(unique(Out$Protein))
				print(paste('There are', Pos, 'peptide(s) kept in Protein', Protein, '(', i, 'out of', N.Prot, ')', sep=' '))	 

      } #F.10 (End of loop for protein)
    } 
    
    if(toupper(Step2) == 'WHOLEPLOT'){

      N.Prot <- length(unique(Profile.Pep$Protein))
      repeated <- .checkRepeated(work)
      Out2 <- data.frame(Protein=vector(), Peptide=vector(), Error.Port=vector(), Error.Pep=vector())				
      Index.2 <- 0				

      for(k in 1:N.Prot){ #F.11
				
				sub <- subset(Profile.Pep, Protein == unique(Protein)[k])
    
				sub$GROUP <- factor(sub$GROUP)
        sub$SUBJECT <- factor(sub$SUBJECT)
    		sub$GROUP_ORIGINAL <- factor(sub$GROUP_ORIGINAL)	
        sub$SUBJECT_ORIGINAL <- factor(sub$SUBJECT_ORIGINAL)
        sub$SUBJECT_NESTED <- factor(sub$SUBJECT_NESTED)
        sub$RUN <- factor(sub$RUN)
    
        singleSubject <- .checkSingleSubject(sub)
        TechReplicate <- .checkTechReplicate(sub) ## use for label-free model
    
				## case-control
        if (!repeated) {
          if (!TechReplicate | singleSubject) {
            fit.full <- lm(Run.Intensity ~ GROUP , data = sub)
          } else {
            fit.full <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = sub)
          }
        } else { ## time-course
          if (singleSubject) {
            fit.full <- lm(Run.Intensity ~ GROUP , data = sub)
          } else { ## no single subject
        	  if (!TechReplicate) {
              fit.full <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = sub)
            } else {
              fit.full <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT), data = sub) ## SUBJECT==SUBJECT_NESTED here
            }
          }	
        } ## time-course

        N.Pep <- length(unique(sub$Peptide))
				sub$Res <- summary(fit.full)$res
        
        for(j in 1:N.Pep){ #F.12

          Out2[(Index.2+1), 'Protein'] <- as.character(unique(sub$Protein))
          Out2[(Index.2+1), 'Peptide'] <- as.character(unique(sub$Peptide)[j])
					Residual <- sub[sub$Peptide==unique(sub$Peptide)[j], 'Res']
					Out2[(Index.2+1), 'Error.Prot'] <- var(Residual, na.rm=TRUE)

          subsub <- subset(sub, Peptide==unique(Peptide)[j])

					## case-control
          if (!repeated) {
            if (!TechReplicate | singleSubject) {
              fit.sub <- lm(Run.Intensity ~ GROUP , data = subsub)
            } else {
              fit.sub <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = subsub)
            }
          } else { ## time-course
            if (singleSubject) {
              fit.sub <- lm(Run.Intensity ~ GROUP , data = subsub)
            } else { ## no single subject
              if (!TechReplicate) {
                fit.sub <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) , data = subsub)
              } else {
                fit.sub <- lmer(Run.Intensity ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT), data = subsub) ## SUBJECT==SUBJECT_NESTED here
              }
            }	
          } ## time-course

          Out2[(Index.2+1), 'Error.Pep'] <- var(summary(fit.sub)$res, na.rm=TRUE)
          Index.2 <- Index.2+1

				} #F.12 End of loop for peptides

				#Now conduct the selection process for each protein

			}#F.11 End of loop for proteins
    } #F9.2
 	
    work.keep$Protein_Peptide <- paste(work.keep$PROTEIN, work.keep$PEPTIDE, sep='_')
    work.keep2 <- subset(work.keep, Protein_Peptide %in% Keep.Pep)
    work.keep2$Protein_Feature <- NULL; work.keep2$Filter.Repro <- NULL	

    ##If we do not allow any of the proteins being totally removed, we keep the most abundant peptide if it is totally removed
		if( !remove_proteins_with_interference ){ #S.7
			
      work <- work.0
      Prot2 <- unique(work.keep2$PROTEIN)
      Prot1 <- unique(work$PROTEIN)
      Loss <- setdiff(Prot1, Prot2)
      work$Protein_Peptide <- paste(work$PROTEIN, work$PEPTIDE, sep='_') 
      Keep.Peptide <- vector()
      IK <- 0

      if (length(Loss) > 0) {
        sub.Loss <- subset(work, PROTEIN %in% Loss)
        N2.Prot <- length(unique(sub.Loss$PROTEIN)) 
					
        for(j in 1:N2.Prot){

          subsub.Loss <- subset(sub.Loss, PROTEIN==unique(sub.Loss$PROTEIN)[j])
					N.Pep2 <- length(unique(subsub.Loss$PEPTIDE))
          
					#Get the most abundant peptide of this protein
					AVG.Pep <- tapply(subsub.Loss$ABUNDANCE, subsub.Loss$Protein_Peptide, function(x) mean(x, na.rm=TRUE))
					Top1 <- names(AVG.Pep[rank(-AVG.Pep)==1])
					Keep.Peptide[IK+1] <- Top1
          IK <- IK+1
				}
        work.keepPep <- subset(work, Protein_Peptide %in% Keep.Peptide)
        work.4 <- rbind(work.keep2, work.keepPep)
        work.4$Protein_Feature <- NULL
        work.keep2 <- work.4

      } else {
        #The case of no proteins was totally removed
        work.keep2$Protein_Feature <- NULL
        work.keep2$Protein.Peptide <- NULL
			}
		} else {

      work.keep2$Protein_Feature <- NULL
      work.keep2$Protein.Peptide <- NULL

    } #S.7	
			
		work <- work.keep2
			
	} #2.4 : End of case for label-based experiment (Remove interference)

	#Factorize some important variables because some of the levels have been removed
	work$PROTEIN <- factor(work$PROTEIN, levels=unique(work$PROTEIN))
	work$PEPTIDE <- factor(work$PEPTIDE, levels=unique(work$PEPTIDE))
	work$FEATURE <- factor(work$FEATURE, levels=unique(work$FEATURE))
	return(work)

} #End of function '.feature_selection'

