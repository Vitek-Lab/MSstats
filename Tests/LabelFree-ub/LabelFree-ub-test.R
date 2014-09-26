#! /usr/bin/Rscript --vanilla

rm(list = ls(all = TRUE)) 

library(ggplot2)

source('MSstats.daily/R/mainfunctions.R')
source('MSstats.daily/R/methods.R')
source('Code/MQtoMSstatsFormat.R')

## assume a maxquant evidence file with at least the following columns
## Raw file|Intensity|Proteins|Modifications|Sequence|Modified sequence|Charge|Protein group IDs|id|Retention time|Reverse|Contaminant
evidence_sample = read.delim("Tests/LabelFree-ub/data/LabelFree-ub-evidence-sample.txt", header=T, stringsAsFactors=F)

## filter out ub-sites only by looking for peptides with di-gly (gl) modified lysines 
evidence_sample = evidence_sample[grepl('(gl)',evidence_sample$Modified.sequence),]

## assume an annotation file with the following columns
## Raw.file|IsotopeLabelType|Condition|BioReplicate|Run
annotation = read.delim("Tests/LabelFree-ub/data/LabelFree-ub-evidence-keys.txt", header=T, stringsAsFactors=F)

##  converting maxquant format and annotation file to MSstats format
mssdata  = MaxQtoMSstatsFormat(evidence = evidence_sample, annotation =  annotation, proteinGroups=NULL, proteinID="Proteins", useUniquePeptide=T, summary=max, fewMeasurements="remove", experiment="DDA")

## profile plots before normalization
mssquant = dataProcess(mssdata, normalization=F, fillIncompleteRows=T)
dataProcessPlots(data=mssquant,type="ProfilePlot",featureName="Peptide", address="Tests/LabelFree-ub/data/LabelFree-ub-unadjusted-")

## normalization and other choices for data processing
mssquant = dataProcess(mssdata, normalization='globalStandards', nameStandards=c('P04406'), fillIncompleteRows=T)
dataProcessPlots(data=mssquant,type="ProfilePlot",featureName="Peptide", address="Tests/LabelFree-ub/data/LabelFree-ub-normalized-")

## reading the file with comparisons between conditions
comparisons = read.delim(file = 'Tests/LabelFree-ub/data/LabelFree-ub-evidence-contrasts.txt')
results = groupComparison(data = mssquant, contrast.matrix = comparisons, labeled = F, scopeOfBioReplication = 'restricted', scopeOfTechReplication = 'restricted', interference = F, equalFeatureVar = F, missing.action = 'impute')
