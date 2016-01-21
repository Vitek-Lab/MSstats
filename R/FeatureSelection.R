
###########Function of FeatureSelection###########
#Input is the data frame 'work' after conducting normalization
#Output is the data frame 'work' after removing noisy features 
#This process removes the irreproducible features across the MS runs


.feature_selection <- function(work, LeaveOneOut){
	
  	#Label-free experiments (SWATH)
	if(nlevels(work$LABEL)==1 & length(unique(work$TRANSITION))>1){ #1.2
    
       #We can merely work on one peptide at a time
       N.Prot <- nlevels(work$PROTEIN) 
       
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

             #First, remove the bottom 30% fragments of each peptide
             M1 <- tapply(sub2$ABUNDANCE, sub2$FEATURE, function(x) mean(x, na.rm=TRUE)) 
             Top70.Cut <- round(length(M1)*0.7) 
             Keep.Top70 <- names(M1[rank(-M1)<=Top70.Cut])

             #Here we only work on the top 70% most abundant fragments in each peptide
             sub2.2 <- subset(sub2, FEATURE %in% Keep.Top70)
             sub2.2$FEATURE <- factor(sub2.2$FEATURE, levels=unique(sub2.2$FEATURE))
             N.Feature <- nlevels(sub2.2$FEATURE) 
             sub2.2$RUN <- factor(sub2.2$RUN, levels=unique(sub2.2$RUN))

             #Use TMP to do robust run quantification
             data_tmp = dcast(RUN ~ FEATURE, data=sub2.2, value.var='ABUNDANCE', keep=TRUE)
  		 rownames(data_tmp) <- data_tmp$RUN
  	       data_tmp <- data_tmp[,-1]
  	       data_tmp[data_tmp==1]<-NA

  		 TMP <- medpolish(data_tmp, na.rm=TRUE, eps = 0.01, maxiter = 100, trace.iter = FALSE)
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

          #Overall, the 30% most irreproducible features are gone
          Interference.Score$Protein_Feature <- paste(Interference.Score$Protein, Interference.Score$Feature, sep='_')
          Interference.Score$Rank.Overall <- rank(Interference.Score$Interference) 
          Top70.2 <- round(dim(Interference.Score)[1]*0.7)
          
          Interfere <- Interference.Score[Interference.Score$Rank <= Top70.2,]

          #Whichever peptides only has one feature left should be removed as it indicates this peptide is irreproducibility
          Interfere$Protein_Peptide <- paste(Interfere$Protein, Interfere$Peptide, sep='_')
          Table <- table(Interfere$Protein_Peptide)
          Drop <- names(Table[Table==1])
          Interfere <- subset(Interfere, !(Protein_Peptide %in% Drop))

          #Now use the remaining features to compute the average amount of interference on each peptide
          #Retain 50% most reproducible peptide. We prepare the case where not all of the peptides are proteotypic
          Keep <- unique(Interfere$Protein_Feature)
          work$Protein_Feature <- paste(work$PROTEIN, work$FEATURE, sep='_')
          work.2 <- subset(work, Protein_Feature %in% Keep)
          work.2$Protein_Peptide <- paste(work.2$PROTEIN, work.2$PEPTIDE, sep='_')

          
          N.Peptide <- length(unique(work.2$Protein_Peptide))

          #Create a data frame storing the average amount of interference at each peptide
          Inter.Pep <- data.frame(Protein_Peptide=vector(), Avg.Inter=vector())
          Index.Pep <- 0

             #Work on one protein at a time. 
             for(h in 1:N.Peptide){  #F.5
  
             temp1 <- subset(work.2, Protein_Peptide==unique(Protein_Peptide)[h])
             N.temp1 <- length(unique(temp1$FEATURE))           

             #Use TMP to do robust run quantification
             data_pep = dcast(RUN ~ FEATURE, data=temp1, value.var='ABUNDANCE', keep=TRUE)
  		 rownames(data_pep) <- data_pep$RUN
  	       data_pep <- data_pep[,-1]
  	       data_pep[data_pep==1]<-NA

  		 TMP2 <- medpolish(data_pep, na.rm=TRUE, eps = 0.01, maxiter = 100, trace.iter = FALSE)         
             
             #Create a vector to save the interference score for each feature
             Inter <- vector()

                #Calculate the score of interference on each fragment
                #In fact, var(data_tmp[,k]-TMP.Run)=var(TMP$res[,k])
                for(k in 1:N.temp1){ #F.6
             
                  Inter[k] <- var(TMP2$res[,k], na.rm=TRUE)

                                      } #F.6
                #End of loop for quantifying the amount of interference at each feature at this peptide
 
                Inter.Pep[(Index.Pep+1), 'Protein_Peptide'] <- as.character(unique(temp1$Protein_Peptide))                 
                Inter.Pep[(Index.Pep+1), 'Avg.Inter'] <- mean(Inter, na.rm=TRUE)               

                Index.Pep <- Index.Pep+1
                                   } #F.5
             #End of loop for computing the average amount of interference at each peptide

          #Retain the best 40% peptides with most reproducible quantification
          M2 <- quantile(Inter.Pep$Avg.Inter, na.rm=TRUE, probs=0.4)
          Keep40 <- Inter.Pep[Inter.Pep$Avg.Inter <= M2, 'Protein_Peptide']
          work.3 <- subset(work.2, Protein_Peptide %in% Keep40)

          work.3$Protein_Feature <- NULL; work.3$Protein_Peptide <- NULL

          #data frame 'work.3' contain the selected features and 'work' is the data frame to be returned 
          work <- work.3

	   
		} # 1.2 : End of label-free SWATH experiment

		#Label-free experiments (Shotun; DDA)
		if(nlevels(work$LABEL)==1 & nlevels(work$TRANSITION)==1){ #1.2.2

			# In DDA (shotgun) experiment, there is no need to filter out bad fragments in each peptide
			N.Prot <- length(unique(work$PROTEIN))

			# Create the data frame storing the selected features of each Protein
			FeatureSelection.Out <- data.frame(Protein=vector(), Peptide=vector(), Model.Based.Error=vector(), Rank=vector(), Filter=vector())
			Index.FS <- 0

                  #Decide whether the experiment is time-course or case control
                  TC <- .checkRepeated(work)

			#Determine if there is only one subject and if there is any technical replicates
                  Single.Subject <- .checkSingleSubject(work)
			Technical.Rep <- .checkTechReplicate(work)

			###Rank the peptides for each protein
			for(i in 1:N.Prot){ #DDA_1.3

				DDA.1 <- subset(work, PROTEIN == unique(PROTEIN)[i])
				N.Pep <- length(unique(DDA.1$PEPTIDE))

				## show progress
				message(paste("Selection features or peptides for protein ", unique(DDA.1$PROTEIN), "(", i, " of ", N.Prot, ")"))
				
				## Create a data frame storing the assessment of the noise for each peptide
				Pep.Out <- data.frame(Peptide=rep(NA, N.Pep), Model.Based.Error=rep(NA,N.Pep), Flag=rep(NA,N.Pep))

				## Assess the noise based on model for each peptide##

				for(j in 1:N.Pep){ # DDA_1.4

					DDA.Pep <- subset(DDA.1, PEPTIDE==unique(PEPTIDE)[j])

                              if(TC & !Single.Subject){ #DDA_1.4.1
                              #Time-course experiment with more than one bio-replicates

                              	if(Technical.Rep){ #DDA_1.4.1.1
                                    #With technical replicates                              	
                                    LM.DDA <- try(lm(ABUNDANCE ~ GROUP * SUBJECT, data=DDA.Pep), silent=TRUE)

                                    } else {
                                    #Without technical replicates
                              	LM.DDA <- try(lm(ABUNDANCE ~ GROUP + SUBJECT, data=DDA.Pep), silent=TRUE)

                                    } #DDA_1.4.1.1
                              } else {
                              #Case-Control

                                    if(Technical.Rep & !Single.Subject){ #DDA_1.4.1.2
                                    #With technical replicates
                                    LM.DDA <- try(lm(ABUNDANCE ~ GROUP + SUBJECT, data=DDA.Pep), silent=TRUE)

                                    } else {
                        		#Without technical replicates or with only one bio-replicates
   						LM.DDA <- try(lm(ABUNDANCE ~ GROUP, data=DDA.Pep), silent=TRUE)

                                    } #DDA_1.4.1.2
                              } #DDA_1.4.1

					Pep.Out[j, 'Peptide'] <- as.character(unique(DDA.Pep$PEPTIDE))
					
					## If there are missing values such that only one condition is left, the model based error score is not calculable
					if(class(LM.DDA)=='try-error'){ #DDA_1.5
						Pep.Out[j, 'Model.Based.Error'] <- NA
					} else {
						Pep.Out[j, 'Model.Based.Error'] <- summary(LM.DDA)$sigma
                	} #DDA_1.5
                	
				} # DDA_1.4 (End of loop for each peptide)

				## We choose the top one-third of the peptides in DDA now. (Again, deciding the optimal number of the peptides we should choose takes too long. We are working on this.)
				## We automatically excluded the precursors whose model-based error is not able to be quantified. This would be due to too many missing values such that only one condition has the peaks
				N.Select <- max(round(dim(Pep.Out[!(is.na(Pep.Out$Model.Based.Error)),])[1]/3+0.01), 1)
				Pep.Out$Rank <- rank(Pep.Out$Model.Based.Error)
				Pep.Out[Pep.Out$Rank <= N.Select, 'Flag'] <- 'Selected'
				Pep.Out[!(Pep.Out$Rank <= N.Select), 'Flag'] <- 'Noisy'
				Pep.Out[is.na(Pep.Out$Model.Based.Error), 'Flag'] <- 'Noisy'


				## Pull out the results to the designated data frame, and flag the noisy peptide
				FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Pep), 'Protein'] <- as.character(unique(DDA.1$PROTEIN))
				FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Pep), 'Peptide'] <- Pep.Out$Peptide
				FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Pep), 'Model.Based.Error'] <- Pep.Out$Model.Based.Error
				FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Pep), 'Rank'] <- Pep.Out$Rank
				FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Pep), 'Filter'] <- Pep.Out$Flag
				Index.FS <- Index.FS+N.Pep   
				                                              
			} # DDA_1.3 (End of loop for proteins)

			## If not all of the peptides are proteotypic, some proteins may have the same peptides
			work$Protein_Peptide <- paste(work$PROTEIN, work$PEPTIDE, sep='_')            
			FeatureSelection.Out$Protein_Peptide <- paste(FeatureSelection.Out$Protein, FeatureSelection.Out$Peptide, sep='_')
          
			## Create the label for selected, or removed, peptides
			work$Filter <- NA; 
			Selection <- unique(FeatureSelection.Out[FeatureSelection.Out$Filter=='Selected','Protein_Peptide'])
			work[work$Protein_Peptide %in% Selection, 'Filter'] <- 'Selected'
			work[!(work$Protein_Peptide %in% Selection), 'Filter'] <- 'Flagged'
			work$Protein_Peptide <- NULL
			work <- subset(work, Filter=='Selected')
		
		} #1.2.2  : End of the DDA experiment 

		## Label-based experiments
    	if(nlevels(work$LABEL)==2){ #2.4

    		N.Prot <- length(unique((work$PROTEIN)))

    		#Create a data frame storing the output
    		Out <- data.frame(Protein=vector(), Peptide=vector(), Feature=vector(), Label=vector(), Interference.Score=vector())   
    		Index <- 0
    		
       		#Loop over proteins: One protein at a time
       		for(i in 1:N.Prot){ #2.4.1
         
       			Sub <- subset(work, PROTEIN==unique(PROTEIN)[i])
       			N.Pep <- length(unique(Sub$PEPTIDE))
       			
       			## show progress
				message(paste("Selection features or peptides for protein ", unique(Sub$PROTEIN), "(", i, " of ", N.Prot, ")"))
				
           		#Loop for the peptides in this protein
           		for(j in 1:N.Pep){ #2.4.2

            		Sub2 <- subset(Sub, PEPTIDE==unique(PEPTIDE)[j])    
        
            		#Heavy peptide
            		Sub2.H <- subset(Sub2, LABEL=='H')

            		Sub2.H$RUN <- factor(Sub2.H$RUN, levels=unique(Sub2.H$RUN))
            		Sub2.H$FEATURE <- factor(Sub2.H$FEATURE, levels=unique(Sub2.H$FEATURE))
            
            		#Apply TMP to estimate the run effect
            		data_w = dcast(RUN ~ FEATURE, data=Sub2.H, value.var='ABUNDANCE', keep=TRUE)
  					rownames(data_w) <- data_w$RUN
  	      			data_w <- data_w[, -1]
  	      			data_w[data_w <=1 ] <- NA

  					TMP <- medpolish(data_w, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
	      			TMP.Run <- TMP$overall + TMP$row

            		N.Feature <- length(unique(Sub2.H$FEATURE))

            		#Calculate the variance of the residuals for each feature
                 	for(k in 1:N.Feature){ #2.4.3

                    	#If there is only one feature in this peptide, no reproducibility can be assessed
						if(N.Feature==1){
                            #No need to remove this peptide at this moment
							Error <- 0
						}else{    
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
					data_w.L <- data_w.L[, -1]
					data_w.L[data_w.L <= 1] <- NA

					TMP.L <- medpolish(data_w.L, na.rm=TRUE, eps = 0.01, maxiter = 1000, trace.iter = FALSE)
					TMP.Run.L <- TMP.L$overall + TMP.L$row

					N.Feature <- length(unique(Sub2.L$FEATURE))

					#Calculate the variance of the residuals for each feature on Light
					for(k in 1:N.Feature){ #2.4.3.2

						#If there is only one feature in this peptide, no reproducibility can be assess. The removals of this single features will be decided by the consensus later.
						if(N.Feature==1){
							Error <- 0
						}else{    
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
       
       		#The heuristic approach now is to keep searching the cutoff until the percentage of the improvement is less than 30%, but do not remove more than half of the features and remove at least a quarter of them
            Out.L <- subset(Out, Label=='L'); Out.H <- subset(Out, Label=='H')
            Q <- quantile(Out.L$Interference.Score, probs = seq(0, 1, 0.05), na.rm = TRUE)
            Change <- diff(Q)/Q[-1]; Change.half <- Change[-(1:11)]
            P <- which(Change.half > 0.3); K <- min(P, na.rm=TRUE); P95 <- which(names(Change.half)=='95%')

            name.P <- ifelse(K>P95, '90%', names(Change.half[K-1]))
            name.P <- ifelse(K!=1, name.P, '55%')

            L.Cut <- Q[name.P]

            #For the internal standard, remove no more than 20% of the features as they are supposed to be of good quality
            Q.H <- quantile(Out.H$Interference.Score, probs = seq(0, 1, 0.05), na.rm = TRUE)
            Change.H <- diff(Q.H)/Q.H[-1]; Change.H.half <- Change.H[-(1:11)]
            P.H <- which(Change.H.half>0.25); K2 <- min(P.H, na.rm=TRUE); P85 <- which(names(Change.H.half)=='85%')

			if(K2 == Inf){
				name.P.H <- '90%'
            } else {
				name.P.H <- ifelse(K2<P85, '80%', names(Change.H.half[K2-1])); 
      	 		H.Cut <- Q.H[name.P.H]
            }
      
            ##Remove the irreproducible features, which is the first step of the removal.      
       		Out.L$Flag.Repro <- 'OK'; Out.H$Flag.Repro <- 'OK'
       		Out.L[is.na(Out.L$Interference.Score), 'Interference.Score'] <- 10000 
       		Out.H[is.na(Out.H$Interference.Score), 'Interference.Score'] <- 10000 

       		Out.L[Out.L$Interference.Score >= L.Cut, 'Flag.Repro'] <- 'Noisy' 
       		Out.H[Out.H$Interference.Score >= H.Cut, 'Flag.Repro'] <- 'Noisy' 

            ##Extract the remaining features befoer proceeding to the next step
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
 				OneFeature <- names(Freq[Freq==1])	
                work.keep <- subset(work.keep, !(Protein.Peptide %in% OneFeature))
                  #}
                work.keep$Protein.Feature <-NULL; work.keep$Protein.Peptide <- NULL 
	
		if(LeaveOneOut){#LOO		  
			#First it is safer to factorize each categorical variable, for which we use, since some of the levels are gone.
            #This step is precautionary for some potential, but rare, errors.
            work.keep$PROTEIN <- factor(work.keep$PROTEIN); work.keep$FEATURE <- factor(work.keep$FEATURE)
            work.keep$PEPTIDE <- factor(work.keep$PEPTIDE)
            work.keep$GROUP_ORIGINAL <- factor(work.keep$GROUP_ORIGINAL)
            work.keep$SUBJECT_ORIGINAL <- factor(work.keep$SUBJECT_ORIGINAL)
            work.keep$RUN <- factor(work.keep$RUN); work.keep$GROUP <- factor(work.keep$GROUP)

            N.Protein <- length(unique(work.keep$PROTEIN))
            #The normalized intensity after correcting the peptide-by-run interaction
            work.keep$ABUNDANCE.N <- NA

			#Store the results of the leave-one-out, and also store the values of the initial error as E0
            LOO <- data.frame(Protein=vector(), Feature=vector(), Remove_Order=vector(), Improve=vector(), Error=vector())
            Index <- 0; E0 <- vector()

            #Analyze one protein at a time
            for(c in 1:N.Protein){ #2.5.0

            	temp <- subset(work.keep, PROTEIN==unique(PROTEIN)[c])
                N.peptide <- length(unique(temp$PEPTIDE))
				temp$FEATURE <- factor(temp$FEATURE)
                temp$PEPTIDE <- factor(temp$PEPTIDE)

                #Normalize the internal standard for each peptide to adjust the peptide-specific deviation
                for(j in 1:N.peptide){ #2.5.1

					Temp2 <- subset(temp, PEPTIDE==unique(PEPTIDE)[j])
					N.Run <- length(unique(Temp2$RUN))

					#Compute the overall median on internal standard
					Median <- median(Temp2[Temp2$GROUP==0, 'ABUNDANCE'], na.rm=TRUE)

				    #Equalize the median of the internal standard at each run, and apply the same shifts on its counterpart at the endogenous
					for(k in 1:N.Run){ #2.5.2

						temp[temp$PEPTIDE==unique(temp$PEPTIDE)[j] & temp$LABEL=='L' & temp$RUN==k, 'ABUNDANCE.N'] <- temp[temp$PEPTIDE==unique(temp$PEPTIDE)[j] & temp$LABEL=='L' & temp$RUN==k, 'ABUNDANCE']- median(temp[temp$PEPTIDE==unique(temp$PEPTIDE)[j] & temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE'], na.rm=TRUE) + Median
						temp[temp$PEPTIDE==unique(temp$PEPTIDE)[j] & temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE.N'] <- temp[temp$PEPTIDE==unique(temp$PEPTIDE)[j] & temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE']- median(temp[temp$PEPTIDE==unique(temp$PEPTIDE)[j] & temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE'], na.rm=TRUE) + Median

					} #2.5.2

				} #2.5.1 : Finish the loop over the peptides on the adjustment (normalization)

                #Start the leave-one-out algorithm
                list.f <- unique(temp$FEATURE)
                N.Feature <- length(list.f)
                  
				#After the peptides with internal standards are 'normalized', there is only need for modeling the light (without the internal standard)
				temp.L <- subset(temp, GROUP!=0)
				temp.L$RUN <- factor(temp.L$RUN); temp.L$FEATURE <- factor(temp.L$FEATURE)

                #Estimate of the random error when every feature presents
                data.TMP <- dcast(RUN ~ FEATURE, data=temp.L, value.var='ABUNDANCE.N', keep=TRUE)
				data.TMP <- data.TMP[,-1]
                data.TMP[data.TMP <=1] <- NA
                Out.TMP <- medpolish(data.TMP, maxiter=100, trace.iter=FALSE, na.rm=TRUE)
                Error.0 <- sd(Out.TMP$residuals, na.rm=TRUE); E0[c] <- Error.0

				#leave-one-out can only be meaningfully conducted when there are more than 2 features
                if(N.Feature >2){ #2.5.A
                	W.Pep <- data.frame(Protein=vector(), Feature=vector(), Remove_Order=vector(), Improve=vector(), Error=vector())

					#Record the initial error estimates
                    W.Pep[1, 'Protein'] <- as.character(unique(temp.L$PROTEIN))
                    W.Pep[1, 'Feature'] <- 'All Features'; W.Pep[1, 'Remove_Order'] <- NA
			      	W.Pep[1, 'Error'] <- round(Error.0, digits=4); W.Pep[1, 'Improve'] <- NA
                        
					#Determine the worst feature one at a time
                    #The number of iterations is at most (N.Feature-2) as the leave-one-out cannot be conducted when 2 features left
					for(iter in 1:(N.Feature-2)){ #2.5.2
						#The feature contributed the most when it is removed
						Error.I <- vector()
                        Error.Decline <- vector()
                        K <- N.Feature-iter+1
                              
						for(I in 1:K){ #2.5.3

							Leave <- list.f[I]
                            sub3 <- subset(temp.L, FEATURE!=Leave)
                            M_sub3 = dcast(RUN ~ FEATURE, data=sub3, value.var='ABUNDANCE.N', keep=TRUE)
                            M_sub3 <- M_sub3[,-1]; M_sub3[M_sub3 <= 1] <- NA
                            Out.sub3 <- medpolish(M_sub3, maxiter=1000, trace.iter=FALSE, na.rm=TRUE)
                            Error.I[I] <- sd(Out.sub3$res, na.rm=TRUE)
							Error.Decline[I] <- (Error.0-Error.I[I])/Error.0

                        } #2.5.3

						#Index of the worst feature
                        Location <- which(Error.Decline==max(Error.Decline, na.rm=TRUE))
                        Worst.Feature <- list.f[Location]
						Error.Worst <- Error.I[Location]
                        Increment <- Error.Decline[Location]

                        W.Pep[(iter+1), 'Protein'] <- as.character(unique(temp.L$PROTEIN))
                        W.Pep[(iter+1), 'Feature'] <- as.character(Worst.Feature)
                        W.Pep[(iter+1), 'Remove_Order'] <- iter
						W.Pep[(iter+1), 'Error'] <- round(Error.Worst, digits=4)
                        W.Pep[(iter+1), 'Improve'] <- round(Increment, digits=3)
						#At the end of this iteration, the worst feature is removed from the data frame
                        Error.0 <- Error.Worst
                        temp.L <- subset(temp.L, FEATURE != Worst.Feature)
 						list.f <- unique(temp.L$FEATURE)                             
 
                    }#2.5.2 (End of iterations of leave-one-out in this protein)
                  
					n.row <- dim(W.Pep)[1] 
                  	LOO[(Index+1):(Index+n.row),] <- W.Pep
                  	Index <- Index + n.row                 
 
				} else { #The case of two-feature proteins
                  	LOO[(Index+1),] <- c(as.character(unique(temp.L$PROTEIN)), 'No more than 2 features', NA, NA, Error.0)
                  	Index <- Index+1
                } #2.5.A

				print(paste(c, 'out of total', N.Protein, 'proteins (', as.character(unique(W.Pep$Protein)), ')are done for leave-one-out.', sep=' '))

			} #2.5.0
            #End of loop for proteins       
         
                  
			#Remove the features with obvious interference
            Error0 <- quantile(E0, probs=seq(0,1,0.05), na.rm=TRUE)
            Improve <- quantile(as.numeric(LOO$Improve), probs=seq(0,1,0.05), na.rm=TRUE)
            E40 <- Error0['40%']; E20 <- Error0['20%']; E60 <- Error0['60%']; E30 <- Error0['30%']
			IM60 <- Improve['60%']; IM80 <- Improve['80%']; IM85 <- Improve['85%']
            Remove.Protein <- vector(); Remove.Feature <- vector(); Keep.Protein <- vector(); Consistent.Protein <- vector()
            Index.RP <- 0; Index.RF <- 0; Index.KP <- 0; Index.CP <- 0
			LOO$P.F <- paste(LOO$Protein, LOO$Feature, sep='.')
			LOO$Peptide.1 <-  sub("(.*?)_(.*?)_.*", "\\1", LOO$Feature)
			LOO$Peptide.2 <-  sub("(.*?)_(.*?)_.*", "\\2", LOO$Feature)
			LOO$Peptide <-  paste(LOO$Peptide.1, LOO$Peptide.2, sep='_')

			for(k in 1:N.Protein){ #2.6

				sub <- subset(LOO, Protein==unique(LOO$Protein)[k])
                Name <- as.character(unique(sub$Protein))
                DIM <- dim(sub)[1]; N.Pep <- length(unique(sub))
                N.Pep <- length(unique(sub$Peptide))-1
 
				#Case of 2-feature protein
				if(DIM == 1){ #2.7
			
     				if(sub$Error >= E40){ #2.7.1
						Remove.Protein[Index.RP+1] <- Name
						Index.RP <- Index.RP+1
					} #2.7.1
	                #Case of 3-feature protein
				} else if(DIM == 2) { 

					if(sub$Error[2] >=E40){ #2.7.2
						Remove.Protein[Index.RP+1] <- Name
      					Index.RP <- Index.RP+1

					} else if(sub$Error[1] < E20){
						Keep.Protein[Index.KP+1] <- Name
						Index.KP <- Index.KP+1

					} else if(sub$Error[1] > E40 & sub$Error[2] <= E40){
						Remove.Feature[Index.RF+1] <- as.character(sub$P.F[2])
  					} #2.7.2

				#Case of more than 3 features	(more than 1 peptide)
				} else {

					if(sub$Error[1] <= E30){ #2.7.3
						Keep.Protein[Index.KP+1] <- Name
						Index.KP <- Index.KP+1
					
					} else if(sub$Error[1] > E30 & sub$Error[1] <= E60){

						#Large improvement at the end of any peptide is an indication of different but consistent profiles from different peptides
						Pep <- sub[-1,]; Key <- vector(); Index.K <- 0
						
						for(p in 1:N.Pep){

							Pep.1 <- subset(Pep, Peptide==unique(Pep$Peptide)[p])
							
							if(dim(Pep.1)[1] > 1){
								Key[Index.K+1] <- Pep.1$Improve[length(Pep.1$Improve)]
							} else {
								Key[Index.K+1] <- 0
							}
							Index.K <- Index.K+1

						}

						if(max(Key, na.rm=TRUE) >= IM80){
							Consistent.Protein[Index.CP+1] <- Name
							Index.CP <- Index.CP+1
						} else {
							suppressWarnings(I <- min(which(Pep$Error < E30), na.rm=TRUE))

							if(I==Inf){ #Remove the whole protein if the error never be small
								Remove.Protein[Index.RP+1] <- Name
								Index.RP <- Index.RP+1
							} else {
                                Gone <- Pep[1:(I-1), 'P.F']; G <- length(Gone)
                                Remove.Feature[(Index.RF+1):(Index.RF+G)] <- as.character(Gone) 
								Index.RF <- Index.RF + G
							}
						}
							
					} else if(sub$Error[1] > E60){

						#Large improvement at the end of any peptide is an indication of different but consistent profiles from different peptides
						Pep <- sub[-1,]; Key <- vector(); Index.K <- 0
						
						for(p in 1:N.Pep){

							Pep.1 <- subset(Pep, Peptide==unique(Pep$Peptide)[p])
							if(dim(Pep.1)[1] > 1){
								Key[Index.K+1] <- Pep.1$Improve[length(Pep.1$Improve)]
							} else {
								Key[Index.K+1] <- 0
							}
							Index.K <- Index.K+1

						}

						if(max(Key, na.rm=TRUE) >= IM85){
							Consistent.Protein[Index.CP+1] <- Name
							Index.CP <- Index.CP+1 

						} else {

							suppressWarnings(I40 <- min(which(Pep$Error <= E40), na.rm=TRUE))
							if(I40==Inf){
								Remove.Protein[Index.RP+1] <- Name
								Index.RP <- Index.RP+1
							} else{
                                Gone2 <- Pep[1:(I-1), 'P.F']; G2 <- length(Gone2)
                                Remove.Feature[(Index.RF+1):(Index.RF+G2)] <- as.character(Gone2) 
								Index.RF <- Index.RF + G2
							}
						}

					} #2.7.3

				} #2.7
                      
			} #2.6 End of loop of the proteins

			Remove.P <- setdiff(setdiff(Remove.Protein, Keep.Protein), Consistent.Protein)
		            		
       		work.keep.LOO <- subset(work.keep, !(PROTEIN %in% Remove.P))
			work.keep.LOO$Protein.Feature <- paste(work.keep.LOO$PROTEIN, work.keep.LOO$FEATURE, sep='.')
			work.keep.LOO2 <- subset(work.keep.LOO, !(Protein.Feature %in% Remove.Feature))
			work.keep.LOO2$Protein_Feature <- NULL; work.keep.LOO2$ABUNDANCE.N <- NULL; work.keep.LOO2$Protein.Feature <- NULL
       		work <- work.keep.LOO2
		} #LOO

		work <- work.keep
	} #2.4 : End of case for label-based experiment (Remove interference)

	#Factorize some important variables because some of the levels have been removed
	work$PROTEIN <- factor(work$PROTEIN, levels=unique(work$PROTEIN))
	work$PEPTIDE <- factor(work$PEPTIDE, levels=unique(work$PEPTIDE))
	work$FEATURE <- factor(work$FEATURE, levels=unique(work$FEATURE))
	return(work)

} #End of function '.feature_selection'


