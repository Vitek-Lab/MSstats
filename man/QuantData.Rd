\name{QuantData}
\alias{QuantData}
\title{Quantitative data after data pre-processing and quality control of MS runs}
\description{This is an example of quantitative data after data pre-processing and quality control of MS runs (i.e., transformation, normalization, and filtering of the original intensities measurements) and can be obtained directly from \code{\link{dataProcess}}. }

\usage{QuantData}

\format{data.frame}

\details{The quantitative data after data pre-processing and quality control of MS runs not only contain the same variable from the raw data, but with addition variables for statistical model fitting and group comparison. For examples, Variable ABUNDANCE represents the final measurement, which could be normalized or not depending on the options you specified in \code{\link{dataProcess}}. Default option in \code{\link{dataProcess}} is with log2 transformation and normalization.

}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. Protein significance analysis in selected reaction monitoring (SRM) measurements.  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.}

\examples{
QuantData<-dataProcess(RawData)

head(QuantData)
}

