
###########Function of FeatureSelection###########
#Input is the data frame 'work' after conducting normalization
#Output is the data frame 'work' after removing noisy features 
#featureSelectionGoal is either 'TPR.TNR' or 'RemoveInterference'
#'TPR.TNR' selects the features more likely to enhance the true positive rate


.feature_selection <- function(work, featureSelectionGoal='TPR.TNR'){

	if(featureSelectionGoal=='TPR.TNR'){ #1.1

  		#Label-free experiments (SWATH)
		if(nlevels(work$LABEL)==1 & length(unique(work$TRANSITION))>1){ #1.2

			N.Prot <- nlevels(work$PROTEIN)

			#Create the data frame storing the selected features of each Protein
			FeatureSelection.Out <- data.frame(Protein=vector(), Peptide=vector(), Feature=vector())
			Index.FS <- 0
  
			#Conduct the feature selection for each protein
			for(i in 1:N.Prot){ #1.3

				tem <- subset(work, PROTEIN==unique(PROTEIN)[i])
				N.Pep <- length(unique(tem$PEPTIDE)) 
				
				#Create a data frame storing the assessment of the noise for each peptide
				Pep.Out <- data.frame(Peptide=rep(NA, N.Pep), Model.Based.Error=rep(NA,N.Pep))

				#Also, create the data frames storing the retained features in each peptide
				Retained.Feature <- data.frame(Peptide=vector(), Feature=vector()); Index.f <- 0

        		#In each peptide, select the max(Top 3, Top 40%) of the features based on the within-group variation
        		#Then within-group variation is calculated, after blocking the features, for each peptide
        		for(j in 1:N.Pep){ #1.4
           
           			tem2 <- subset(tem, PEPTIDE==unique(PEPTIDE)[j])
           			tem2$FEATURE <- factor(tem2$FEATURE, levels=unique(tem2$FEATURE))

					#Determine how many features to pick in this peptide
           			N.Feature <- length(unique(tem2$FEATURE)); Top40 <- round(N.Feature*0.4)         
          			N.Select <- max(3, Top40)
          			
           			#Create a data frame storing the variance(error) for each feature 
           			F.Out <- data.frame(Feature=rep(NA, N.Feature), Error=rep(NA, N.Feature))

         			#The peptide with less than 3 fragments will not be considered.
         			if(N.Feature>=N.Select){ #1.4.1
         				
            			#The features are ranked by the within-group variation        
            			for(k in 1:N.Feature){ #1.5
   
            				tem.f <- subset(tem2, FEATURE==unique(FEATURE)[k])
            				LM.f <- lm(ABUNDANCE ~ GROUP, data=tem.f)
            				F.Out[k, 'Feature'] <- as.character(unique(tem.f$FEATURE))
            				F.Out[k, 'Error'] <- summary(LM.f)$sigma
            
                        } #1.5 
            
           				F.Out$Rank <- rank(F.Out$Error)

            			#We work on the max(Top3, Top40) features of this peptide
           				Keep.Feature <- F.Out[F.Out$Rank <= N.Select,'Feature']
            			tem2.2 <- subset(tem2, FEATURE %in% Keep.Feature)

            			#Compute the with-group variation, after blocking for features, of this peptide
            			tem2.2$FEATURE <- factor(tem2.2$FEATURE, levels=unique(tem2.2$FEATURE))
            			LM.Pep <- lm(ABUNDANCE ~ FEATURE + GROUP, data=tem2.2)
            
            			Pep.Out[j, 'Peptide'] <- as.character(unique(tem2.2$PEPTIDE))
            			Pep.Out[j, 'Model.Based.Error'] <- summary(LM.Pep)$sigma

            			#Store which fragments were retained in this peptide
            			Retained.Feature[(Index.f+1):(Index.f+N.Select), 'Peptide'] <- as.character(unique(tem2$PEPTIDE))
            			Retained.Feature[(Index.f+1):(Index.f+N.Select), 'Feature'] <- Keep.Feature

            			#Update the current index at 'Retained.Feature'
            			Index.f <- Index.f+N.Select
 
         			} else {

         				#The case of peptides with less than 3 features
            			if(N.Feature > 1){ #1.4.2
            
            				LM.Pep2 <- lm(ABUNDANCE ~ FEATURE + GROUP, data=tem2)
            				Pep.Out[j, 'Peptide'] <- as.character(unique(tem2$PEPTIDE))
         
         					#Penalize the peptides having less than 3 features
            				Pep.Out[j, 'Model.Based.Error'] <- (1.5)*summary(LM.Pep2)$sigma

            				#Store which fragments were retained in this peptide
           					Retained.Feature[(Index.f+1):(Index.f+N.Feature), 'Peptide'] <- as.character(unique(tem2$PEPTIDE))
            				Retained.Feature[(Index.f+1):(Index.f+N.Feature), 'Feature'] <- as.character(unique(tem2$FEATURE))

            				#Update the current index at 'Retained.Feature'
            				Index.f <- Index.f+N.Feature

            			} else {

            				LM.Pep2 <- lm(ABUNDANCE ~ GROUP, data=tem2)
            				Pep.Out[j, 'Peptide'] <- as.character(unique(tem2$PEPTIDE))
            				Pep.Out[j, 'Model.Based.Error'] <- 2*summary(LM.Pep2)$sigma

            				#Store which fragments were retained in this peptide
            				Retained.Feature[(Index.f+1):(Index.f+N.Feature), 'Peptide'] <- as.character(unique(tem2$PEPTIDE))
            				Retained.Feature[(Index.f+1):(Index.f+N.Feature), 'Feature'] <- as.character(unique(tem2$FEATURE))

            				#Update the current index at 'Retained.Feature'
            				Index.f <- Index.f+N.Feature
            				
            			} #1.4.2

        				warning(paste('The Peptide', unique(tem2$PEPTIDE), 'is likely to be removed since it has less than 3 fragments. This is our rule imposed in the SWATH experiment.', sep=' '))

					} #1.4.1

				} #1.4


          		## Select the best one-fourth of the peptide now. (Deciding the optimal number of the peptides chosen takes too much time, but we are working on reducing its computation time.)
          		N.Select.Pep <- max(round(dim(Pep.Out)[1]/4), 1)
          		Pep.Out$Rank <- rank(Pep.Out$Model.Based.Error)   
          		Selected.Pep <- Pep.Out[Pep.Out$Rank <= N.Select.Pep, 'Peptide']
          		
          		## Pull out the selected features in this peptide        
          		Selected.Feature <- subset(Retained.Feature, Peptide %in% Selected.Pep)     
          
          		## Put the selected features (peptides) of this proteins into the desinated data frame
          		N.F <- dim(Selected.Feature)[1]
          		FeatureSelection.Out[(Index.FS+1):(Index.FS+N.F), 'Protein'] <- as.character(unique(tem$PROTEIN))
          		FeatureSelection.Out[(Index.FS+1):(Index.FS+N.F), 'Peptide'] <- Selected.Feature$Peptide
          		FeatureSelection.Out[(Index.FS+1):(Index.FS+N.F), 'Feature'] <- Selected.Feature$Feature

          		Index.FS <- Index.FS+N.F 
			} #1.3 (end of loop for proteins)

			## If not all of the peptides are proteotypic, some proteins may have the same peptides
			work$Protein_Feature <- paste(work$PROTEIN, work$FEATURE, sep='_')            
			FeatureSelection.Out$Protein_Feature <- paste(FeatureSelection.Out$Protein, FeatureSelection.Out$Feature, sep='_')
          
			## Create the label for selected, or removed, features
			work$Filter <- NA; 
			work[work$Protein_Feature %in% unique(FeatureSelection.Out$Protein_Feature), 'Filter'] <- 'Selected'
			work[!(work$Protein_Feature %in% unique(FeatureSelection.Out$Protein_Feature)), 'Filter'] <- 'Flagged'
			work$Protein_Feature <- NULL
		  
			###If we want to remove the bad features directly
			work <- subset(work, Filter=='Selected')
           
		} # 1.2 : End of label-free SWATH experiment

		#Label-free experiments (Shotun; DDA)
		if(nlevels(work$LABEL)==1 & nlevels(work$TRANSITION)==1){ #1.2.2

			# In DDA (shotgun) experiment, there is no need to filter out bad fragments in each peptide
			N.Prot <- nlevels(work$PROTEIN)

			# Create the data frame storing the selected features of each Protein
			FeatureSelection.Out <- data.frame(Protein=vector(), Peptide=vector())
			Index.FS <- 0

			###Rank the peptides for each protein
			for(i in 1:N.Prot){ #DDA_1.3

				DDA.1 <- subset(work, PROTEIN == unique(PROTEIN)[i])
				N.Pep <- length(unique(DDA.1$PEPTIDE))

				#$ Create a data frame storing the assessment of the noise for each peptide
				Pep.Out <- data.frame(Peptide=rep(NA, N.Pep), Model.Based.Error=rep(NA,N.Pep))

				## Assess the within-group variation for each peptide##

				for(j in 1:N.Pep){ # DDA_1.4

					DDA.Pep <- subset(DDA.1, PEPTIDE==unique(PEPTIDE)[j])
					LM.DDA <- lm(ABUNDANCE ~ GROUP, data=DDA.Pep)
					Pep.Out[j, 'Peptide'] <- as.character(unique(DDA.Pep$PEPTIDE))
					Pep.Out[j, 'Model.Based.Error'] <- summary(LM.DDA)$sigma
                  
				} # DDA_1.4 (End of loop for each peptide)

				###We choose the top one-third of the peptides in DDA now. (Again, deciding the optimal number of the peptides we should choose takes too long. We are working on this.)
				N.Select <- max(round(length(unique(DDA.1$PEPTIDE))/3), 1)
				Pep.Out$Rank <- rank(Pep.Out$Model.Based.Error)
				Pep.Select <- Pep.Out[Pep.Out$Rank <= N.Select, 'Peptide']

				###Pull out the selected peptides to the designated data frame
				FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Select), 'Protein'] <- as.character(unique(DDA.1$PROTEIN))
				FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Select), 'Peptide'] <- Pep.Select
				Index.FS <- Index.FS+N.Select   
				                                              
			} # DDA_1.3 (End of loop for proteins)

			## If not all of the peptides are proteotypic, some proteins may have the same peptides
			work$Protein_Peptide <- paste(work$PROTEIN, work$PEPTIDE, sep='_')            
			FeatureSelection.Out$Protein_Peptide <- paste(FeatureSelection.Out$Protein, FeatureSelection.Out$Peptide, sep='_')
          
			###Create the label for selected, or removed, peptides
			work$Filter <- NA; 
			work[work$Protein_Peptide %in% unique(FeatureSelection.Out$Protein_Peptide), 'Filter'] <- 'Selected'
			work[!(work$Protein_Peptide %in% unique(FeatureSelection.Out$Protein_Peptide)), 'Filter'] <- 'Flagged'
			work$Protein_Peptide <- NULL
			work <- subset(work, Filter=='Selected')
		
		} #1.2.2  : End of the DDA experiment (TPR.TNR)


		#Label-based experiments (TPN.TNR)
		if(nlevels(work$LABEL)==2){#L.1

			N.Prot <- nlevels(work$PROTEIN)

			#Create a data frame to store the error score of each protein and feature
			Out.L <- data.frame(Protein=vector(), Feature=vector(), Error.score=vector(), Label=vector(), Rank=vector(), Flag=vector())  
			Out.H <- data.frame(Protein=vector(), Feature=vector(), Error.score=vector(), Label=vector(), Rank=vector(), Flag=vector())  
			Index <- 0

			#Loop over the proteins
			for(i in 1:N.Prot){ #L.2

				sub <- subset(work, PROTEIN==unique(PROTEIN)[i])
				N.Pep <- length(unique(sub$PEPTIDE)) 
				sub$ABUNDANCE.N <- NA     

				#Normalize the heavy label for each peptide as peptide-by-run interaction is commonly seen
				for(j in 1:N.Pep){ #L.3

					temp <- subset(sub, PEPTIDE==unique(PEPTIDE)[j])
					N.Run <- length(unique(temp$RUN))

					#Compute the overall median on heavy
					Median <- median(temp[temp$LABEL=='H', 'ABUNDANCE'], na.rm=FALSE)
         
             
					#Normalize the light on each run 
					for(k in 1:N.Run){ #L.4

						temp[temp$LABEL=='L' & temp$RUN==k, 'ABUNDANCE'] <- temp[temp$LABEL=='L' & temp$RUN==k, 'ABUNDANCE']- median(temp[temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE'], na.rm=TRUE) + Median
						temp[temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE'] <- temp[temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE']- median(temp[temp$LABEL=='H' & temp$RUN==k, 'ABUNDANCE'], na.rm=TRUE) + Median

					} #L.4

					sub[sub$PEPTIDE==unique(sub$PEPTIDE)[j], 'ABUNDANCE.N'] <- temp$ABUNDANCE

				} #L.3 : Finish the loop over the peptides

            	#After normalizing the heavy, the batch effect should have been corrected.
             	#Assess the within-group variation for each feature
            
                N.Feature <- length(unique(sub$FEATURE))
                sub$FEATURE <- factor(sub$FEATURE, levels=unique(sub$FEATURE))

                #Create a data frame to store the error score of each feature
                F.Out.L <- data.frame(Feature=vector(), Error.score=vector(), Label=vector())  
                F.Out.H <- data.frame(Feature=vector(), Error.score=vector(), Label=vector())  

                for(m in 1:N.Feature){  #L.5

                	Feature <- subset(sub, FEATURE==unique(FEATURE)[m])

                	#Work on heavy and light separately as the model is different
                	Feature.H <- subset(Feature, LABEL=='H')
                	Feature.L <- subset(Feature, LABEL=='L')

                	#If there are less than 2 abundances in one group, the linear model has error message.
                	LM.L <- try(lm(ABUNDANCE.N ~ GROUP_ORIGINAL, data=Feature.L),silent=TRUE)
                	
                  	if(class(LM.L)=="try-error") {  #L.5.1
						Sigma.L <- NA
					} else {
						Sigma.L <- summary(LM.L)$sigma   
					}  #L.5.1

                	LM.H <- try(suppressWarnings(lm(ABUNDANCE.N ~ 1, data=Feature.H)),silent=TRUE)
                	
                  	if(class(LM.H)=="try-error") {  #L.5.2
						Sigma.H <- NA
					} else {
						Sigma.H <- suppressWarnings(summary(LM.H))$sigma   
					}  #L.5.2
                         
					F.Out.L[m, 'Feature'] <- as.character(unique(Feature$FEATURE))
                	F.Out.L[m, 'Error.score'] <- Sigma.L; F.Out.L[m, 'Label'] <- 'L'

                	F.Out.H[m, 'Feature'] <- as.character(unique(Feature$FEATURE))
                	F.Out.H[m, 'Error.score'] <- Sigma.H; F.Out.H[m, 'Label'] <- 'H' 

				}   #L.5 : End of loop for computing the error scores at each feature

                #Rank the features based on the error scores. Rank 1 means the most noisy.
                F.Out.L$Rank <- rank(-F.Out.L$Error.score)
                F.Out.H$Rank <- rank(-F.Out.H$Error.score)

                #Flag the worst one-third features on heavy and light, respectively 
                #The reproducibility-optimized cutoff requires much running time so we do not use it for now.
                Cut.L <- round(dim(F.Out.L)[1]/3); Cut.H <- round(dim(F.Out.H)[1]/3) 

                F.Out.L$Flag <- "OK"; F.Out.H$Flag <- "OK"
                F.Out.L[F.Out.L$Rank<=Cut.L, 'Flag'] <- "Noisy"
                F.Out.H[F.Out.H$Rank<=Cut.H, 'Flag'] <- "Noisy"

                #The features with too many missing values so that the error scores are NA are also flagged
                F.Out.L[is.na(F.Out.L$Error.score), 'Flag'] <- "Noisy"
                F.Out.H[is.na(F.Out.H$Error.score), 'Flag'] <- "Noisy"

                #Put the results in a data frame which contains the results for all proteins
                Out.L[(Index+1):(Index+N.Feature), 'Protein'] <- as.character(unique(sub$PROTEIN))
                Out.L[(Index+1):(Index+N.Feature), 2:6] <- F.Out.L

                Out.H[(Index+1):(Index+N.Feature), 'Protein'] <- as.character(unique(sub$PROTEIN))
                Out.H[(Index+1):(Index+N.Feature), 2:6] <- F.Out.H

                #Update the index
                Index <- Index+N.Feature

			} #L.2 : End of loop over the proteins

                 
            #For shared peptide it creates same peptide sequence, so use protein-feature instead
            Out.L$Prot_F <- paste(Out.L$Protein, Out.L$Feature, sep='_')
            Out.H$Prot_F <- paste(Out.H$Protein, Out.H$Feature, sep='_')

            #Only the features both selected at Heavy and Light are kept
            L.list <- Out.L[Out.L$Flag=='OK', 'Prot_F']
            H.list <- Out.H[Out.H$Flag=='OK', 'Prot_F']
            Keep <- intersect(L.list, H.list)

            work$Prot_F <- paste(work$PROTEIN, work$FEATURE, sep='_')

            work$Filter <- NA
            work[work$Prot_F %in% Keep, 'Filter'] <- 'Selected'
            work[!(work$Prot_F %in% Keep), 'Filter'] <- 'Flagged'

            work$Prot_F <- NULL
            work <- subset(work, Filter=='Selected')

		}#L.1 : End of the label-based experiment (TPR.TNR)

	}#1.1 : End of option 'TPR.TNR'


	#Choose the subset of features for removing the interference
	if(featureSelectionGoal=='RemoveInterference'){ #2.1

    	#Label-free experiment
    	if(nlevels(work$LABEL)==1){ #F.1

        	stop("This is a Label-free experiment. Please use \"TPR.TNR\" as the input for the featureSelectionGoal")
                              
		} #F.1

    	#Label-based experiments
    	if(nlevels(work$LABEL)==2){ #2.4

    		N.Prot <- nlevels(work$PROTEIN)

    		#Create a data frame storing the output
    		Out <- data.frame(Protein=vector(), Peptide=vector(), Feature=vector(), Label=vector(), Interference.Score=vector())   
    		Index <- 0
    		
       		#Loop over proteins: One protein at a time
       		for(i in 1:N.Prot){ #2.4.1
         
       			Sub <- subset(work, PROTEIN==unique(PROTEIN)[i])
       			N.Pep <- length(unique(Sub$PEPTIDE))
       			
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
  	      			data_w <- data_w[,-1]
  	      			data_w[data_w==1]<-NA

  					TMP <- medpolish(data_w, na.rm=TRUE, eps = 0.01, maxiter = 100, trace.iter = FALSE)
	      			TMP.Run <- TMP$overall + TMP$row

            		N.Feature <- length(unique(Sub2.H$FEATURE))

            		#Calculate the variance of the residuals for each feature
                 	for(k in 1:N.Feature){ #2.4.3

                    	#If there is only one feature in this peptide, no reproducibility can be assess
						if(N.Feature==1){
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
					} #2.4.3 : End of loop for calculating the interference score for each feature


            		#Light peptide
            		Sub2.L <- subset(Sub2, LABEL=='L')

            		Sub2.L$RUN <- factor(Sub2.L$RUN, levels=unique(Sub2.L$RUN))
            		Sub2.L$FEATURE <- factor(Sub2.L$FEATURE, levels=unique(Sub2.L$FEATURE))
            
            		#Apply TMP to estimate the run effect
            		data_w.L = dcast(RUN ~ FEATURE, data=Sub2.L, value.var='ABUNDANCE', keep=TRUE)
					rownames(data_w.L) <- data_w.L$RUN
					data_w.L <- data_w.L[,-1]
					data_w.L[data_w.L==1]<-NA

					TMP.L <- medpolish(data_w.L, na.rm=TRUE, eps = 0.01, maxiter = 100, trace.iter = FALSE)
					TMP.Run.L <- TMP.L$overall + TMP.L$row

					N.Feature <- length(unique(Sub2.L$FEATURE))

					#Calculate the variance of the residuals for each feature on Light
					for(k in 1:N.Feature){ #2.4.3.2

						#If there is only one feature in this peptide, no reproducibility can be assess
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
					} #2.4.3.2 : End of loop for calculating the interference score for each feature on Light     

				} #2.4.2 : End of loop for the peptides
       
			} #2.4.1 :End of loop over proteins

        	
      		#Preparing for the case of shared peptide, we use the column Protein_Feature
       		Out$Prot.F <- paste(Out$Protein, Out$Feature, sep='_') 
       
       		#Flag the worst one-third of the light transition and one-fourth of the heavy
       		#Again, how to optimally decide the cutoff is an ongoing research. It is a trade-off between optimum and expensive calculation

       		Out.L <- subset(Out, Label=='L'); Out.H <- subset(Out, Label=='H')
       		L.Cut <- quantile(Out.L$Interference.Score, probs = (2/3), na.rm=TRUE)
       		H.Cut <- quantile(Out.H$Interference.Score, probs = (3/4), na.rm=TRUE)
       
       		Out.L$Flag <- 'OK'; Out.H$Flag <- 'OK'
       		Out.L[is.na(Out.L$Interference.Score), 'Interference.Score'] <- 10000 
       		Out.H[is.na(Out.H$Interference.Score), 'Interference.Score'] <- 10000 

       		Out.L[Out.L$Interference.Score >= L.Cut, 'Flag'] <- 'Noisy' 
       		Out.H[Out.H$Interference.Score >= H.Cut, 'Flag'] <- 'Noisy' 

       		Keep.L <- Out.L[Out.L$Interference.Score < L.Cut, 'Prot.F']
       		Keep.H <- Out.H[Out.H$Interference.Score < H.Cut, 'Prot.F']
       
       		#Only keep the features when both of heavy and light are retained
       		Keep <- intersect(Keep.L, Keep.H)
       		work$Protein_Feature <- paste(work$PROTEIN, work$FEATURE, sep='_')

       		work$Filter <- NA
       		work[work$Protein_Feature %in% Keep, 'Filter'] <- 'Selected'
       		work[!(work$Protein_Feature %in% Keep), 'Filter'] <- 'Flagged'
       		work$Protein_Feature <- NULL
       		work <- subset(work, Filter=='Selected')

		} #2.4 : End of case for label-based experiment (Remove interference)

	} #2.1 : End of the option 'RemoveInterference'

	return(work)

} #End of function 'FeatureSelection'

#featureSelectionGoal is either 'TPR.TNR' or 'RemoveInterference'

