\name{designSampleSize}
\alias{designSampleSize}
\title{Planning future experimental designs of Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiments in sample size calculation}
\description{Calculate sample size for future experiments of a Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment based on intensity-based linear model. Two options of the calculation: (1) number of biological replicates per condition, (2) power.}
\usage{
designSampleSize(data=data,desiredFC=desiredFC,FDR=0.05,numSample=TRUE,power=0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{'fittedmodel' in testing output from function groupComparison.}
  \item{desiredFC}{the range of a desired fold change which includes the lower and upper values of the desired fold change.}
  \item{FDR}{a pre-specified false discovery ratio (FDR) to control the overall false positive. Default is 0.05}
  \item{numSample}{minimal number of biological replicates per condition. TRUE represents you require to calculate the sample size for this category, else you should input the exact number of biological replicates.}
  \item{power}{a pre-specified statistical power which defined as the probability of detecting a true fold change. TRUE represent you require to calculate the power for this category, else you should input the average of power you expect. Default is 0.9}
}

\details{The function fits the model and uses variance components to calculate sample size. The underlying model fitting with intensity-based linear model with technical MS run replication. Estimated sample size is rounded to 0 decimal.
}

\value{A list of the sample size calculation results including Variable desiredFC, numSample, numPep, numTran, FDR,  and power.}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}

\references{
Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\section{Warning }{It can only obtain either one of the categories of the sample size calculation (numSample, numPep, numTran, power) at the same time. 
}

\examples{
# Consider quantitative data (i.e. QuantData) from yeast study.
# A time course study with ten time points of interests and three biological replicates.

QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

## based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
comparison<-rbind(comparison1,comparison2, comparison3)
row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")

testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=QuantData)


## Calculate sample size for future experiments:

#(1) Minimal number of biological replicates per condition

designSampleSize(data=testResultMultiComparisons$fittedmodel,numSample=TRUE,
desiredFC=c(1.25,1.75),FDR=0.05,power=0.8)


#(2) Power calculation

designSampleSize(data=testResultMultiComparisons$fittedmodel,numSample=2,
desiredFC=c(1.25,1.75),FDR=0.05,power=TRUE)
}

