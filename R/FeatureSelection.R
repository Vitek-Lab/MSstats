
###########Function of FeatureSelection###########
#Input is the data frame 'work' after conducting normalization
#Output is the data frame 'work' after removing noisy features 
#featureSubset is either 'highQuality_Significance' or 'highQuality_Reproducibility'
#'highQuality_Significance' selects the features more likely to enhance the true positive rate
#'highQuality_Reproducibility' selects features with reproducible quantification. The with peptide like PTM which has different pattern but reproducible will likely be kept in this options. 



.feature_selection <- function(work, featureSubset){

	if(toupper(featureSubset)=='HIGHQUALITY_SIGNIFICANCE'){ #1.1

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
				
				## show progress
				message(paste("Selection features or peptides for protein ", unique(tem$PROTEIN), "(", i, " of ", N.Prot, ")"))
				
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
			N.Prot <- length(unique(work$PROTEIN))

			# Create the data frame storing the selected features of each Protein
			FeatureSelection.Out <- data.frame(Protein=vector(), Peptide=vector(), Model.Based.Error=vector(), Rank=vector(), Filter=vector())
			Index.FS <- 0

			###Rank the peptides for each protein
			for(i in 1:N.Prot){ #DDA_1.3

				DDA.1 <- subset(work, PROTEIN == unique(PROTEIN)[i])
				N.Pep <- length(unique(DDA.1$PEPTIDE))

				## show progress
				message(paste("Selection features or peptides for protein ", unique(DDA.1$PROTEIN), "(", i, " of ", N.Prot, ")"))
				
				## Create a data frame storing the assessment of the noise for each peptide
				Pep.Out <- data.frame(Peptide=rep(NA, N.Pep), Model.Based.Error=rep(NA,N.Pep), Flag=rep(NA,N.Pep))

				## Assess the within-group variation for each peptide##

				for(j in 1:N.Pep){ # DDA_1.4

					DDA.Pep <- subset(DDA.1, PEPTIDE==unique(PEPTIDE)[j])
					LM.DDA <- try(lm(ABUNDANCE ~ GROUP, data=DDA.Pep), silent=TRUE)
					Pep.Out[j, 'Peptide'] <- as.character(unique(DDA.Pep$PEPTIDE))
					
					## If there are missing values such that only one condition is left, the model based error score is not calculable
					if(class(LM.DDA)=='try-error'){ #DDA_1.5
						Pep.Out[j, 'Model.Based.Error'] <- NA
					} else {
						Pep.Out[j, 'Model.Based.Error'] <- summary(LM.DDA)$sigma
                	} #DDA_1.5
                	
				} # DDA_1.4 (End of loop for each peptide)

				## We choose the top one-third of the peptides in DDA now. (Again, deciding the optimal number of the peptides we should choose takes too long. We are working on this.)
				## We automatically excluded the precursors whose model-based error is not able to be quantified. This would be there are a lot of missing peaks at this precursor
				N.Select <- max(round(dim(Pep.Out[!(is.na(Pep.Out$Model.Based.Error)),])[1]/3), 1)
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
		
		} #1.2.2  : End of the DDA experiment (Significance)


		## Label-based experiments (TPN.TNR)
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
				
				## show progress
				message(paste("Selection features or peptides for protein ", unique(sub$PROTEIN), "(", i, " of ", N.Prot, ")"))
				
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

            	## After normalizing the heavy, the batch effect should have been corrected.
             	## Assess the within-group variation for each feature
            
                N.Feature <- length(unique(sub$FEATURE))
                sub$FEATURE <- factor(sub$FEATURE, levels=unique(sub$FEATURE))

                ## Create a data frame to store the error score of each feature
                F.Out.L <- data.frame(Feature=vector(), Error.score=vector(), Label=vector())  
                F.Out.H <- data.frame(Feature=vector(), Error.score=vector(), Label=vector())  

                for(m in 1:N.Feature){  #L.5

                	Feature <- subset(sub, FEATURE==unique(FEATURE)[m])

                	## Work on heavy and light separately as the model is different
                	Feature.H <- subset(Feature, LABEL=='H')
                	Feature.L <- subset(Feature, LABEL=='L')

                	## If there are less than 2 abundances in one group, the linear model has error message.
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

                ## Rank the features based on the error scores. Rank 1 means the most noisy.
                F.Out.L$Rank <- rank(-F.Out.L$Error.score)
                F.Out.H$Rank <- rank(-F.Out.H$Error.score)

                ## Flag the worst one-third features on heavy and light, respectively 
                ## The reproducibility-optimized cutoff requires much running time so we do not use it for now.
                Cut.L <- round(dim(F.Out.L)[1]/3); Cut.H <- round(dim(F.Out.H)[1]/3) 

                F.Out.L$Flag <- "OK"; F.Out.H$Flag <- "OK"
                F.Out.L[F.Out.L$Rank<=Cut.L, 'Flag'] <- "Noisy"
                F.Out.H[F.Out.H$Rank<=Cut.H, 'Flag'] <- "Noisy"

                ## The features with too many missing values so that the error scores are NA are also flagged
                F.Out.L[is.na(F.Out.L$Error.score), 'Flag'] <- "Noisy"
                F.Out.H[is.na(F.Out.H$Error.score), 'Flag'] <- "Noisy"

                ## Put the results in a data frame which contains the results for all proteins
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

		}#L.1 : End of the label-based experiment (Significance)

	}#1.1 : End of option 'TPR.TNR'


	#Choose the subset of features for removing the interference
	if(toupper(featureSubset)=='HIGHQUALITY_REPRODUCIBILITY'){ #2.1

    	#Label-free experiment
    	if(nlevels(work$LABEL)==1 & nlevels(work$TRANSITION)==1){ #F.1
      #DDA 
        	stop("This is a Label-free DDA experiment. Please use \"TPR.TNR\" as the input for the featureSelectionGoal. DDA experiment only quantifies until the peptide level so the reproducibility of each peptide cannot be assessed.")
                              
		} else if(nlevels(work$LABEL)==1 & nlevels(work$TRANSITION)>1) {
      #DIA  

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

          

                                                                }#F.1
        #End of label-free experiment on removing irreproducible features

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

	} #2.1 : End of the option 'HighQuality_Reproducibility'

	return(work)

} #End of function '.feature_selection'

#featureSubset is either 'HighQuality_Significance' or 'HighQuality_Reproducibility'

