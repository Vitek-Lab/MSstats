
###########Function of FeatureSelection###########

FeatureSelection <- function(work, featureSelectionGoal='TPR.TNR'){

if(featureSelectionGoal='TPR.TNR'){ #1.1

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
 
         } else{

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

            } else{

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


          ###Select the best one-fourth of the peptide now. (Deciding the optimal number of the peptides chosen takes too much time, but we are working on reducing its computation time.)
          N.Select.Pep <- max(round(dim(Pep.Out)[1]/4), 1)
          Pep.Out$Rank <- rank(Pep.Out$Model.Based.Error)   
          Selected.Pep <- Pep.Out[Pep.Out$Rank <= N.Select.Pep, 'Peptide']
          ###Pull out the selected features in this peptide        
          Selected.Feature <- subset(Retained.Feature, Peptide %in% Selected.Pep)     
          
          ###Put the selected features (peptides) of this proteins into the desinated data frame
          N.F <- dim(Selected.Feature)[1]
          FeatureSelection.Out[(Index.FS+1):(Index.FS+N.F), 'Protein'] <- as.character(unique(tem$PROTEIN))
          FeatureSelection.Out[(Index.FS+1):(Index.FS+N.F), 'Peptide'] <- Selected.Feature$Peptide
          FeatureSelection.Out[(Index.FS+1):(Index.FS+N.F), 'Feature'] <- Selected.Feature$Feature

          Index.FS <- Index.FS+N.F 
                        } #1.3 (end of loop for proteins)

          ###If not all of the peptides are proteotypic, some proteins may have the same peptides
          work$Protein_Feature <- paste(work$PROTEIN, work$FEATURE, sep='_')            
          FeatureSelection.Out$Protein_Feature <- paste(FeatureSelection.Out$Protein, FeatureSelection.Out$Feature, sep='_')
          
          ###Create the label for selected, or removed, features
          work$Filter <- NA; 
          work[work$Protein_Feature %in% unique(FeatureSelection.Out$Protein_Feature), 'Filter'] <- 'Selected'
          work[!(work$Protein_Feature %in% unique(FeatureSelection.Out$Protein_Feature)), 'Filter'] <- 'Flagged'
          work$Protein_Feature <- NULL
		  
		  ###If we want to remove the bad features directly
		  work <- subset(work, Filter=='Selected')
           
                              } #1.2 (end of label-free SWATH experiment)



   #Label-free experiments (Shotun; DDA)
   if(nlevels(work$LABEL)==1 & nlevels(work$TRANSITION)==1){ #1.2.2

   #In DDA (shotgun) experiment, there is no need to filter out bad fragments in each peptide
   N.Prot <- nlevels(work$PROTEIN)

   #Create the data frame storing the selected features of each Protein
   FeatureSelection.Out <- data.frame(Protein=vector(), Peptide=vector())
   Index.FS <- 0

       ###Rank the peptides for each protein
       for(i in 1:N.Prot){ #DDA_1.3

       DDA.1 <- subset(work, PROTEIN == unique(PROTEIN)[i])
       N.Pep <- length(unique(DDA.1$PEPTIDE))

       #Create a data frame storing the assessment of the noise for each peptide
       Pep.Out <- data.frame(Peptide=rep(NA, N.Pep), Model.Based.Error=rep(NA,N.Pep))

           
           ##Assess the within-group variation for each peptide##

           for(j in 1:N.Pep){ #DDA_1.4

           DDA.Pep <- subset(DDA.1, PEPTIDE==unique(PEPTIDE)[j])
           LM.DDA <- lm(ABUNDANCE ~ GROUP, data=DDA.Pep)
           Pep.Out[j, 'Peptide'] <- as.character(unique(DDA.Pep$PEPTIDE))
           Pep.Out[j, 'Model.Based.Error'] <- summary(LM.DDA)$sigma
                  
                             } #DDA_1.4 (End of loop for each peptide)

           ###We choose the top one-third of the peptides in DDA now. (Again, deciding the optimal number of the peptides we should choose takes too long. We are working on this.)
           N.Select <- max(round(length(unique(DDA.1$PEPTIDE))/3), 1)
           Pep.Out$Rank <- rank(Pep.Out$Model.Based.Error)
           Pep.Select <- Pep.Out[Pep.Out$Rank <= N.Select, 'Peptide']

           ###Pull out the selected peptides to the designated data frame
           FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Select), 'Protein'] <- as.character(unique(DDA.1$PROTEIN))
           FeatureSelection.Out[(Index.FS+1):(Index.FS+N.Select), 'Peptide'] <- Pep.Select
           Index.FS <- Index.FS+N.Select                                                 
                          } #DDA_1.3 (End of loop for proteins)


          ###If not all of the peptides are proteotypic, some proteins may have the same peptides
          work$Protein_Peptide <- paste(work$PROTEIN, work$PEPTIDE, sep='_')            
          FeatureSelection.Out$Protein_Peptide <- paste(FeatureSelection.Out$Protein, FeatureSelection.Out$Peptide, sep='_')
          
          ###Create the label for selected, or removed, peptides
          work$Filter <- NA; 
          work[work$Protein_Peptide %in% unique(FeatureSelection.Out$Protein_Peptide), 'Filter'] <- 'Selected'
          work[!(work$Protein_Peptide %in% unique(FeatureSelection.Out$Protein_Peptide)), 'Filter'] <- 'Flagged'
          work$Protein_Peptide <- NULL
          work <- subset(work, Filter=='Selected')
                                                                 } #1.2.2 (End of the DDA experiment)








   #Label-based experiments
   if(nlevels(work$LABEL)==2){#1.2.3

        #There is an issue for Quality control samples to deal with. Will update this part soon.

                              }#1.2.3

                                } #1.1




#Choose the subset of features for removing the interference
if(featureSelectionGoal='Interference'){ #2.1


###Update this option later



                                         } #2.1


                                   } #End of function 'FeatureSelection'