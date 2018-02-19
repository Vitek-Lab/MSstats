\name{transformMSstatsToMSnSet}
\alias{transformMSstatsToMSnSet}
\title{Transformation input format for MSstats to MSnSet class}
\description{Convert the required input format for MSstats into general format (MSnSet class in MSnbase package) on the proteomics.}

\usage{
transformMSstatsToMSnSet(data)	
}

\arguments{
  \item{data}{name of the raw (input) data set with required column for MSstats.}
}

\details{
	
\itemize{
\item{raw : See \code{\link{SRMRawData}} for the required data structure of raw (input) data.}
\item{output : After transformation, assayData includes value of Intensity. phenoData has variables of IsotopeLabelType, Condition, BioReplicate,Run. featureData has variables of ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge.}
}
}

\value{data.frame as described in details.}

\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}


\examples{
	
library(MSnbase)
quant.msnset<-transformMSstatsToMSnSet(SRMRawData)

}

