\name{designSampleSizeClassification}
\alias{designSampleSizeClassification}
\title{Calculate the optimal size of training data for classification problem}
\description{For classification problem (such as disgnosys of disease), calculate sample size of training data under different size of validation data for future experiments of a Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment based on simulation. }
\usage{
designSampleSizeClassification(data, n = 5, step = 20, noise = 0.0, iter = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{output from function \code{\link{dataProcess}}}
  \item{n}{number of steps, which is the number of different sample sizes to simulate. Default is 5}
  \item{step}{number of samples per condition to increase at each step. Default is 20}
  \item{noise}{the fraction of protein variance to change at each step. Positive value indicates increasing noise in the data while negative value means decreasing noise. Default is 0.0}
  \item{iter}{number of times to repeat simulation experiments. Default is 5}
}

\details{The function fits intensity-based linear model and uses variance components and mean abundance to simulate new data with different sample size. Training data and validation data are simulated independently. Random forest model is fitted on train data and used to predict the validation data. The above procedure is repeated \emph{iter} times. Mean predictive accuracy and variance under different size of training data and validation data are reported. For each size of validation data, the number of training samples with best predictive accuracy is also returned.
}

\value{A list of the sample size calculation results for calssification problem including Variable optTrainSize, meanPA, and varPA. \emph{optTrainSize} is the number of training samples with best predictive accuracy under each size of validation data; \emph{meanPA} is the mean predictive accuracy matrix under different size of training data and validation data; \emph{varPA} is variance of predictive accuracy under different size of training data and validation data.}

\author{Ting Huang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}

\references{
T. Huang et al.  TBD  2018
}


\examples{
# Consider quantitative data (i.e. QuantData) from yeast study.
# A time course study with ten time points of interests and three biological replicates.
QuantSRM <- dataProcess(SRMRawData)

# estimate the mean predictive accuray under different sizes of training data and validation data
# n is the number of biological replicates per condition
result.srm <- designSampleSizeClassification(data=QuantSRM, n=5, step=10)
result.srm$optTrainSize

# Consider another training set from a colorectal cancer study
# Subjects are from control group or colorectal cancer group
# 72 proteins were targeted with SRM
library(MSstatsBioData)
data(SRM_crc_training)
QuantCRCSRM <- dataProcess(SRM_crc_training)

# estimate the mean predictive accuray under different sizes of training data and validation data
# n is the number of biological replicates per condition
crc.srm <- designSampleSizeClassification(data=QuantCRCSRM, n=10, step=10)
crc.srm$optTrainSize
}

