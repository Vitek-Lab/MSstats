\name{designSampleSizeClassificationPlots}
\alias{designSampleSizeClassificationPlots}
\title{Visualization for sample size calculation in classification problem}
\description{To illustrate the mean classification accuracy under different size of training data and validation data. The input is the result from function (\code{\link{designSampleSizeClassification}}.}
\usage{
designSampleSizeClassificationPlots(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{output from function \code{\link{designSampleSizeClassification}}}
}

\details{Data in the example is based on the results of sample size calculation in classification problem from function \code{\link{designSampleSizeClassification}}.
}

\value{Plot for sample size estimation. x-axis : size of training set, y-axis: mean accuracy.}

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
result.srm <- designSampleSizeClassification(data=QuantSRM, n=5, step=10)
designSampleSizeClassificationPlots(data=result.srm)

# Consider another training set from a colorectal cancer study
# Subjects are from control group or colorectal cancer group
# 72 proteins were targeted with SRM
require(MSstatsBioData)
data(SRM_crc_training)
QuantCRCSRM <- dataProcess(SRM_crc_training)

# estimate the mean predictive accuray under different sizes of training data and validation data
# n is the number of biological replicates per condition
crc.srm <- designSampleSizeClassification(data=QuantCRCSRM, n=10, step=10)
designSampleSizeClassificationPlots(data=crc.srm)
}

