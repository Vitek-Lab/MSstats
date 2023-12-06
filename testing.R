library(tidyverse)
library(MSstats)
library(rstan)
library(protDP)

setwd("/Users/kohler.d/Library/CloudStorage/OneDrive-NortheasternUniversity/Northeastern/Research/MS_data/Bulk/DDA_Choi2017")

## Load data
evidence = read.csv("MaxQ/Choi2017_DDA_MaxQuant_evidence.txt", sep="\t")
pg = read.csv("MaxQ/Choi2017_DDA_MaxQuant_proteinGroups.txt", sep="\t")
annotation = read.csv("MaxQ/Choi2017_DDA_MaxQuant_annotation.csv")

msstats_input_data = MaxQtoMSstatsFormat(evidence, annotation, pg,
                                         removeFewMeasurements = FALSE, 
                                         use_log_file = FALSE)
msstats_input_data = as.data.frame(msstats_input_data)
msstats_input_data = msstats_input_data %>% filter(!grepl(";", ProteinName))

# prot = (msstats_input_data %>% distinct(ProteinName) %>% sample_n(50) %>% c())[[1]]
prot = c("P00359")
sample = msstats_input_data %>% filter(ProteinName %in% prot)


summarized_results = dataProcess(sample, normalization = FALSE,
                                 summaryMethod="bayesian",
                                 bayes_method="MCMC", chains=4, cores=4, 
                                 n_iterations=10000, group_size=51, 
                                 use_log_file = FALSE)

dataProcessPlots(summarized_results, type="ProfilePlot")

profile_plot(summarized_results$MSstats, "P00359",
             include_summary=TRUE,
             summary_error_bars=TRUE, color_features_grey=FALSE)

summarized_results$MSstats$ProteinLevelData %>% filter(Protein == "P00360")

g_size = c(10, 25, 50, 100, 250, 500, 1000)

results = list()

times = c()
for (g in g_size){
    t0 = proc.time()
    summarized_results = dataProcess(sample, summaryMethod="bayesian",
                                     bayes_method="MCMC", chains=4, cores=4, 
                                     n_iterations=2000, group_size=g, 
                                     use_log_file = FALSE)
    t1 = proc.time()
    t = t1 - t0
    times = c(times, t[[3]])
    
    results[[as.character(g)]] = summarized_results
}

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")
colnames(comparison)=c("Condition1", "Condition2", "Condition3", "Condition4")

test.MSstats <- groupComparison(contrast.matrix=comparison, data=summarized_results)


