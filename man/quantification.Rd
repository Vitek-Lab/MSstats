\name{quantification}
\alias{quantification}
\title{Protein sample quantification or group quantification}
\description{Model-based quantification for each condition or for each biological samples per protein in a targeted Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment. Quantification takes the processed data set by \code{\link{dataProcess}} as input and automatically generate the quantification results (data.frame) with long or matrix format.
}
\usage{
quantification(data, type="Sample", format="matrix")
}

\arguments{
  \item{data}{name of the (processed) data set.}
  \item{type}{choice of quantification. "Sample" or "Group" for protein sample quantification or group quantification.}
  \item{format}{choice of returned format. "long" for long format which has the columns named Protein, Condition, LonIntensities (and BioReplicate if it is subject quantification), NumFeature for number of transitions for a protein, and NumPeaks for number of observed peak intensities for a protein. "matrix" for data matrix format which has the rows for Protein and the columns, which are Groups(or Conditions) for group quantification or the combinations of BioReplicate and Condition (labeled by "BioReplicate"_"Condition") for sample quantification. Default is "matrix" }
}

\details{

\itemize{
\item{Sample quantification : individual biological sample quantification for each protein. The label of each biological sample is a combination of the corresponding group and the sample ID. If there are no technical replicates or experimental replicates per sample, sample quantification is the same as run summarization from dataProcess. If there are technical replicates or experimental replicates, sample quantification is median among run quantification corresponding MS runs.}

\item{Group quantification : quantification for individual group or individual condition per protein. It is median among sample quantification.}

\item{The quantification for endogenous samples is based on run summarization from subplot model, with TMP robust estimation.}

%%\item{results : whether format="long" or "matrix", both files, "Group or SampleQuantification_longformat.csv" and "Group or SampleQuantification_dataMatrix.csv" will be stored in the assigned folder}

The input of this function is the quantitative data from function (\code{\link{dataProcess}}).
}
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}
\references{
Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\examples{
#Consider quantitative data (i.e. QuantData) from a yeast study with ten time points of interests, three biological replicates, and no technical replicates which is a time-course experiment. 
#Sample quantification shows model-based estimation of protein abundance in each biological replicate within each time point.
#Group quantification shows model-based estimation of protein abundance in each time point.

QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

# Sample quantification

sampleQuant<-quantification(QuantData)
head(sampleQuant)

# Group quantification

groupQuant<-quantification(QuantData, type="Group")
head(groupQuant)
}

