\name{designSampleSize}
\alias{designSampleSize}
\title{Planning future experimental designs of label-free and label-based LC-MS, SRM, DIA experiments in sample size calculation}
\description{Calculate sample size for future experiments of a label-free and label-based LC-MS, SRM and DIA experiment based on intensity-based linear model. Four options of the calculation: (1) number of biological replicates per condition, (2) number of peptides per protein, and (3) number of transitions per peptide, (4) power.}
\usage{
designSampleSize(data=data,numSample=numSample,numPep=numPep,numTran=numTran,desiredFC=desiredFC,FDR=FDR,power=power,scopeOfBioReplication="restricted",interference=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{name of the data set.}
  \item{numSample}{minimal number of biological replicates per condition. TRUE represents you require to calculate the sample size for this category, else you should input the exact number of biological replicates.}
  \item{numPep}{minimal number of peptides per protein. TRUE represents you require to calculate the sample size for this category, else you should input the exact number of peptides.}
  \item{numTran}{minimal number of transitions per peptide. TRUE represents you require to calculate the sample size for this category, else you should input the exact number of transitions.}
  \item{desiredFC}{the range of a desired fold change which includes the lower and upper values of the desired fold change.}
  \item{FDR}{a pre-specified false discovery ratio (FDR) to control the overall false positive.}
  \item{power}{a pre-specified statistical power which defined as the probability of detecting a true fold change. TRUE represent you require to calculate the power for this category, else you should input the average of power you expect.}
  \item{scopeOfBioReplication}{choice of scope of biological replication. "restricted"(default) represents restricted scope of biological replication to the selected individuals. "expanded" represents expanded scope of biological replication to the whole population.}
   \item{interference}{choice of interference data. TRUE(default) means data contain interference transitions and need additional model interaction to address the interference. FALSE means data contain no interference transitions and no need additional model interaction to address the interference.}
}

\details{The function fits the model and uses variance components to calculate sample size. The underlying model fitting with intensity-based linear model with technical MS run replication. Data in the example is \code{\link{QuantData}}. Estimated sample size is rounded to 0 decimal.
}

\value{A list of the sample size calculation results including Variable desiredFC, numSample, numPep, numTran, FDR,  and power.}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\section{Warning }{It can only obtain either one of the categories of the sample size calculation (numSample, numPep, numTran, power) at the same time. 
}

\examples{
# Consider quantitative data (i.e. QuantData) from yeast study.
# A time course study with ten time points of interests and three biological replicates.

head(QuantData)

## Calculate sample size for future experiments:

#(1) Minimal number of biological replicates per condition

designSampleSize(data=QuantData,numSample=TRUE,numPep=3,numTran=4,
desiredFC=c(1.25,1.75),FDR=0.05,power=0.8)

#(2) Minimal number of peptides per protein

designSampleSize(data=QuantData,numSample=2,numPep=TRUE,numTran=4,
desiredFC=c(1.25,1.75),FDR=0.05,power=0.8)	

#(3) Minimal number of transitions per peptide

designSampleSize(data=QuantData,numSample=2,numPep=3,numTran=TRUE,
desiredFC=c(1.25,1.75),FDR=0.05,power=0.8)	

#(4) Power calculation

designSampleSize(data=QuantData,numSample=2,numPep=3,numTran=4,
desiredFC=c(1.25,1.75),FDR=0.05,power=TRUE)
}

