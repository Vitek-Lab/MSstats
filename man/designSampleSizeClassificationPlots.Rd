\name{designSampleSizeClassificationPlots}
\alias{designSampleSizeClassificationPlots}
\title{Visualization for sample size calculation in classification problem}
\description{To illustrate the mean classification accuracy under different size of training data and validation data. The input is the result from function (\code{\link{designSampleSizeClassification}}.}
\usage{
designSampleSizeClassificationPlots(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{output from function designSampleSizeClassification.}
}


\details{Data in the example is based on the results of sample size calculation in classification problem from function \code{\link{designSampleSizeClassification}}.
}

\author{Ting Huang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}
\references{
Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\examples{
# Consider quantitative data (i.e. QuantData) from yeast study.
# A time course study with ten time points of interests and three biological replicates.
QuantSRM<-dataProcess(SRMRawData)
# estimate the mean predictive accuray under different sizes of training data and validation data
result.srm<-designSampleSizeClassification(data=QuantSRM, n=5, step=10)
designSampleSizeClassificationPlots(data=result.srm)

# Consider another training set from a colorectal cancer study
# Subjects are from control group or colorectal cancer group
# 72 proteins were targeted with SRM
require(MSstatsBioData)
data(SRM_crc_training)
QuantCRCSRM<-dataProcess(SRM_crc_training)
# estimate the mean predictive accuray under different sizes of training data and validation data
# n is the number of biological replicates per condition
crc.srm<-designSampleSizeClassification(data=QuantCRCSRM, n=10, step=10)
designSampleSizeClassificationPlots(data=crc.srm)
}

