\name{transformMSnSetToMSstats}
\alias{transformMSnSetToMSstats}
\title{Transforms a MSnSet class dataset into a required input for MSstats}
\description{Convert MSnSet class into the required input format for MSstats}

\usage{
transformMSnSetToMSstats(ProteinName,PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge ,IsotopeLabelType, Bioreplicate,Run, Condition, data)	
}

\arguments{
  \item{data}{name of dataset with MSnSet class}
  \item{ProteinName}{name of column in the MSnSet that contains protein information. If not assigned, "ProteinAccession" column will be used.}
  \item{PeptideSequence}{name of column in the MSnSet that contains information of peptide sequence. If not assigned, "PeptideSequence" column will be used.}
  \item{PrecursorCharge}{name of column in the MSnSet that contains information of peptide charge. If not assigned, "charge" will be used.}
  \item{FragmentIon}{name of column in the MSnSet that contains information of transition. If not assigned, value of "NA" will be used.}
  \item{ProductCharge}{name of column in the MSnSet that contains information of transition charge. If not assigned, value of "NA" will be used.}
  \item{IsotopeLabelType}{name of the column in \code{phenoData} component of MSnSet that contains labeling information. If not assigned, "mz" column will be used.}
  \item{Bioreplicate}{name of the column in \code{phenoData} component of MSnSet that contains unique ids of biological replicates of the corresponding samples. If not assigned, rownames of pData(data) will be used.}
  \item{Run}{name of the column in MSnSet that contains information of experimental MS runs. If not assigned, "file" column will be used.}
  \item{Condition}{names of the columns in \code{phenoData} that correspond to the group variables of interest.  If more than one variable is listed, a concatentated variable is created based on the variables.}
}

\details{

raw : See \code{\link[MSnbase]{MSnSet}} for the general format on the proteomics. Condition must be specified. Intensity should not be specified, as this information is extracted automatically from the assayData component of the MSnSet.

}


\section{Warning}{The types of experiment that MSstats can analyze are LC-MS, SRM, DIA(SWATH) with label-free or labeled synthetic peptides. MSstats does not support for metabolic labeling or iTRAQ experiments.
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.

Gatto, L. and Lilly, K.S. (2012). MSnbase-an R Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. \emph{Bioinformatics}, 28, 288-289.
}

\examples{
	
library("MSnbase")
data(itraqdata)
class(itraqdata)

msnset <- quantify(itraqdata[10:15], method = "trap", reporters = iTRAQ4, verbose = FALSE)
msnset

pData(msnset)$group<-c("control","disease","control","disease")

transformMSnSetToMSstats(data=msnset,Condition="group")
}

