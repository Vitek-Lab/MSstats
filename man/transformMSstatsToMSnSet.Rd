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



\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.

Gatto, L. and Lilly, K.S. (2012). MSnbase-an R Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. \emph{Bioinformatics}, 28, 288-289.
}

\examples{
	
library(MSnbase)
quant.msnset<-transformMSstatsToMSnSet(SRMRawData)

}

