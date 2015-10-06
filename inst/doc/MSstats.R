## ---- eval=F-------------------------------------------------------------
#  install.packages(pkgs = 'MSstats_3.0.9.tar.gz', repos = NULL, type = 'source')

## ---- eval=F-------------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("MSstats")

## ---- eval=T, warning=F--------------------------------------------------
library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
?MSstats

## ---- eval=F-------------------------------------------------------------
#  setwd('/Users/Meena/Dropbox/MSstats_GitHub_document')

## ---- eval=T-------------------------------------------------------------
DDA2009.superhirn <- DDARawData

head(DDA2009.superhirn)

## ---- eval=T-------------------------------------------------------------
# default option
DDA2009.TMP <- dataProcess(raw = DDA2009.superhirn,  fillIncompleteRows = TRUE,
                           normalization = 'equalizeMedians',
                           summaryMethod = 'TMP',
                           censoredInt = "NA", cutoffCensored = "minFeatureNRun",
                           MBimpute = TRUE)

## ---- eval=T, echo=T-----------------------------------------------------
names(DDA2009.TMP)

# the data after reformatting and normalization
head(DDA2009.TMP$ProcessedData) 
# run-level summarized data
head(DDA2009.TMP$RunlevelData)
# Since this is not model-based, no model summary (here DDAskyline.quant$ModelQC=NULL).
# Only with 'summaryMethod="linear"'
head(DDA2009.TMP$ModelQC) 
# here 'TMP'
head(DDA2009.TMP$SummaryMethod) 
# predict values by AFT with 'MBimpute=TRUE'. 
# These values are matching with rownames of DDA2009.TMP$ProcessedData
head(DDA2009.TMP$PredictBySurvival) 

## ---- eval=T-------------------------------------------------------------
# no action for missing values.
DDA2009.TMP.random <- dataProcess(raw = DDA2009.superhirn, fillIncompleteRows = TRUE,
                                  normalization = 'equalizeMedians',
                                  summaryMethod = 'TMP',
                                  censoredInt=NULL)

## ---- eval=F, message=F, warning=F---------------------------------------
#  # linear mixed model (lm or lmer) with run and feature
#  DDA2009.linear <- dataProcess(raw = DDA2009.superhirn,
#                                summaryMethod = "linear", censoredInt = NULL)
#  
#  # accerated failure model with left-censored. NA intensities are assumed as censored
#  DDA2009.linear.censored <- dataProcess(raw = DDA2009.superhirn,
#                                         summaryMethod = "linear", censoredInt = "NA")

## ---- eval=F, message=F, warning=F---------------------------------------
#  # use type="QCplot" with all proteins and change the upper limit of y-axis=35
#  dataProcessPlots(data = DDA2009.TMP, type="QCplot", ylimUp=35)

## ---- eval=F-------------------------------------------------------------
#  # QCplot for only 'yeast' protein
#  dataProcessPlots(data = DDA2009.TMP, type="QCplot", ylimUp=35,
#                   which.Protein="yeast", address="yeast_eqmedians_")

## ---- eval=F-------------------------------------------------------------
#  dataProcessPlots(data = DDA2009.TMP, type="Profileplot",  ylimUp=35,
#                   featureName="NA", width=7, height=7, address="DDA2009_TMP_")

## ---- eval=F-------------------------------------------------------------
#  dataProcessPlots(data = DDA2009.TMP.random, type="Profileplot",  ylimUp=35,
#                   featureName="NA", width=7, height=7,
#                   originalPlot=FALSE, summaryPlot=TRUE, address="DDA2009_TMP_random_")

## ------------------------------------------------------------------------
?groupComparison

## ---- eval=TRUE----------------------------------------------------------
levels(DDA2009.TMP$ProcessedData$GROUP_ORIGINAL)

## ---- eval=TRUE----------------------------------------------------------
comparison1<-matrix(c(-1,1,0,0,0,0),nrow=1)
comparison2<-matrix(c(0,-1,1,0,0,0),nrow=1)
comparison3<-matrix(c(0,0,-1,1,0,0),nrow=1)
comparison4<-matrix(c(0,0,0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,0,0,0,-1,1),nrow=1)
comparison6<-matrix(c(1,0,0,0,0,-1),nrow=1)

comparison<-rbind(comparison1,comparison2,comparison3,comparison4,comparison5,comparison6)
row.names(comparison)<-c("C2-C1","C3-C2","C4-C3","C5-C4","C6-C5","C1-C6")

## ---- eval=T, message=F, warning=F---------------------------------------
DDA2009.comparisons <- groupComparison(contrast.matrix = comparison, data = DDA2009.TMP)

## ---- eval=TRUE----------------------------------------------------------
# output from groupComparison function has three data frames
names(DDA2009.comparisons)

## ---- eval=TRUE----------------------------------------------------------
# name of columns in result data.frame
names(DDA2009.comparisons$ComparisonResult) 

## ---- eval=TRUE----------------------------------------------------------
# get only significant proteins and comparisons among all comparisons
SignificantProteins <- with(DDA2009.comparisons, ComparisonResult[ComparisonResult$adj.pvalue < 0.05, ])
nrow(SignificantProteins)

## ---- eval=TRUE----------------------------------------------------------
?groupComparisonPlots

## ---- eval=F-------------------------------------------------------------
#  groupComparisonPlots(data = DDA2009.comparisons$ComparisonResult, type = 'VolcanoPlot')

## ---- eval=F-------------------------------------------------------------
#  groupComparisonPlots(data = DDA2009.comparisons$ComparisonResult, type = 'Heatmap')

## ---- eval=T-------------------------------------------------------------
raw<-DDARawData.Skyline

## ---- eval=T-------------------------------------------------------------
colnames(raw)[9] <- 'Run' 
colnames(raw)[10] <- 'Intensity' 
head(raw)

## ---- eval=T-------------------------------------------------------------
sum(is.na(raw$Intensity))
sum(raw$Intensity==0)

## ---- eval=T-------------------------------------------------------------
library(reshape2)

raw$pepprecursor<-paste(raw$PeptideSequence, raw$PrecursorCharge, sep="_")

data_w = dcast( Run ~ pepprecursor, data=raw, value.var='Intensity', fun.aggregate=sum, fill=NULL) 

newdata = melt(data_w, id.vars=c('Run'))
colnames(newdata)[colnames(newdata) %in% c("variable","value")]<-c('pepprecursor','Intensity')
  		
uniinfo<-unique(raw[,c("ProteinName","PeptideSequence","PrecursorCharge","pepprecursor")])		
newraw<-merge(newdata,uniinfo, by="pepprecursor")

uniinfo<-unique(raw[,c("Run","BioReplicate","Condition")])		
newraw<-merge(newraw,uniinfo, by="Run")

newraw$BioReplicate<-1 # it should be change based on your experiment.
newraw$FragmentIon<-"sum"
newraw$ProductCharge<-NA
newraw$IsotopeLabelType<-"L"

raw<-newraw
# now 'raw' is ready to use MSstats

## ---- eval=F-------------------------------------------------------------
#  ## 1. First, get protein ID information
#  proteinGroups <- read.table("DDA2014_proteinGroups.txt", sep = "\t", header = TRUE)

## ---- eval=F-------------------------------------------------------------
#  ## 2. Read in annotation including condition and biological replicates: annotation.csv
#  annot <- read.csv("DDA2014_annotation.csv", header = TRUE)

## ---- eval=F-------------------------------------------------------------
#  ## 3. Read in MaxQuant file: evidence.txt
#  infile <- read.table("evidence.txt", sep = "\t", header = TRUE)

## ---- eval=F, echo=FALSE-------------------------------------------------
#  load("DDA2014.evidence.RData")
#  infile <- DDA2014.evidence

## ---- eval=F, warning=F--------------------------------------------------
#  ## 4. Reformat for MSstats required input
#  ##    check options for converting format
#  ?MaxQtoMSstatsFormat

## ---- eval=F-------------------------------------------------------------
#  msstats.raw <- MaxQtoMSstatsFormat(evidence=infile, annotation=annot, proteinGroups=proteinGroups)
#  
#  ## now 'msstats.raw' is ready for MSstats

