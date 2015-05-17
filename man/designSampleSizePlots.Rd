\name{designSampleSizePlots}
\alias{designSampleSizePlots}
\title{Visualization for sample size calculation}
\description{To illustrate the relationship of desired fold change and the calculated minimal number sample size which are (1) number of biological replicates per condition, (2) number of peptides per protein, (3) number of transitions per peptide, and (4) power. The input is the result from function (\code{\link{designSampleSize}}.}
\usage{
designSampleSizePlots(data=data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{output from function designSampleSize.}
%%  \item{address}{the name of folder that will store the plots. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory. An output pdf file is automatically created with the default name of "SampleSizePlot.pdf". The command address can help to specify where to store the file as well as how to modify the beginning of the file name.}
}


\details{Data in the example is based on the results of sample size calculation from function \code{\link{designSampleSize}}.
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\examples{
# Based on the results of sample size calculation from function designSampleSize, we generate a series of sample size plots for number of biological replicates, or peptides, or transitions or power plot.

QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

## based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
comparison<-rbind(comparison1,comparison2, comparison3)
row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")

testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=QuantData)


## plot the calculated sample sizes for future experiments:

#(1) Minimal number of biological replicates per condition

result.sample<-designSampleSize(data=testResultMultiComparisons$fittedmodel,numSample=TRUE,desiredFC=c(1.25,1.75),FDR=0.05,power=0.8)
designSampleSizePlots(data=result.sample)

#(2) Power

result.power<-designSampleSize(data=testResultMultiComparisons$fittedmodel,numSample=2,desiredFC=c(1.25,1.75),FDR=0.05,power=TRUE)
designSampleSizePlots(data=result.power)


}

