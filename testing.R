library(tidyverse)
library(MSstats)
library(rstan)
library(protDP)
library(bayesplot)

setwd("D://OneDrive - Northeastern University/Northeastern/Research/MS_data/Bulk/DDA_Choi2017")

## Load data
evidence = read.csv("MaxQ/Choi2017_DDA_MaxQuant_evidence.txt", sep="\t")
pg = read.csv("MaxQ/Choi2017_DDA_MaxQuant_proteinGroups.txt", sep="\t")
annotation = read.csv("MaxQ/Choi2017_DDA_MaxQuant_annotation.csv")

msstats_input_data = MaxQtoMSstatsFormat(evidence, annotation, pg,
                                         removeFewMeasurements = FALSE, 
                                         use_log_file = FALSE)
msstats_input_data = as.data.frame(msstats_input_data)
msstats_input_data = msstats_input_data %>% filter(!grepl(";", ProteinName))

## Find some proteins without missing
no_missing = msstats_input_data %>% group_by(ProteinName) %>% 
    summarise(missing = sum(is.na(Intensity))) %>% filter(missing==0) %>% 
    distinct(ProteinName)
msstats_input_data %>% filter(ProteinName %in% no_missing[[1]]) %>% 
    group_by(ProteinName) %>% summarize(peps = n_distinct(PeptideSequence)) %>% 
    filter(peps > 2)

# prot = (msstats_input_data %>% distinct(ProteinName) %>% sample_n(50) %>% c())[[1]]
prot = c("P00447")
sample = msstats_input_data %>% filter(ProteinName %in% prot)

summarized_results = dataProcess(sample, normalization = FALSE,
                                 summaryMethod="bayesian",
                                 bayes_method="MCMC", chains=4, cores=4, 
                                 n_iterations=10000, group_size=51, 
                                 use_log_file = FALSE)

dataProcessPlots(summarized_results, type="ProfilePlot")

profile_plot(summarized_results$MSstats, "P00447",
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


## determine priors
library(dglm)

## Run effect
run_info = msstats_input_data %>% group_by(ProteinName, Run) %>% 
    summarize(mean_run = mean(log2(Intensity), na.rm=TRUE),
              std_run = sd(log2(Intensity), na.rm=TRUE))
run_info %>% ggplot() + geom_histogram(aes(mean_run), bins=100)
run_info %>% ggplot() + geom_histogram(aes(std_run), bins=100)

## Mean prior
mean(run_info$mean_run, na.rm = TRUE)
sd(run_info$mean_run, na.rm = TRUE)

## Sd prior
fit <- dglm(run_info$std_run~1, family=Gamma(link="log"), 
            mustart=mean(run_info$std_run, na.rm=TRUE))
summary_fit = summary(fit)
mu <- exp(summary_fit$coefficients[[1]])
shape <- exp(abs(summary_fit$dispersion.summary$coefficients[[1]]))
scale <- mu/shape
c(shape, scale)

## Feature effect
msstats_input_data$Feature = paste(msstats_input_data$PeptideSequence,
                                    msstats_input_data$PrecursorCharge,
                                    msstats_input_data$FragmentIon,
                                    msstats_input_data$ProductCharge, sep="_")
feature_info = msstats_input_data %>% group_by(ProteinName, Feature) %>% 
    summarize(mean_feat = mean(log2(Intensity), na.rm=TRUE),
              std_feat = sd(log2(Intensity), na.rm=TRUE))
feature_info %>% ggplot() + geom_histogram(aes(mean_feat), bins=100)
feature_info %>% ggplot() + geom_histogram(aes(std_feat), bins=100)

## Mean prior
mean(feature_info$mean_feat, na.rm = TRUE)
sd(feature_info$mean_feat, na.rm = TRUE)

## Sd prior
fit <- dglm(feature_info$std_feat~1, family=Gamma(link="log"), 
            mustart=mean(feature_info$std_feat, na.rm=TRUE))
summary_fit = summary(fit)
mu <- exp(summary_fit$coefficients[[1]])
shape <- exp(abs(summary_fit$dispersion.summary$coefficients[[1]]))
scale <- mu/shape
c(shape, scale)

## Scale
msstats_input_data %>% group_by(Run) %>% 
    summarize(std = sd(log2(Intensity), na.rm=TRUE)) %>% 
    ggplot() + geom_histogram(aes(std), bins=5)


## TODO:
## Same analysis for features
## Same analysis for overall std (either singular/feature/run based)
## Test on swapping out feature vs run vs single var dims
