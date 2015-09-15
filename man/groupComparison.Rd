\name{groupComparison}
\alias{groupComparison}
\title{Finding differentially abundant proteins across conditions in targeted Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment}
\description{Tests for significant changes in protein abundance across conditions based on a family of linear mixed-effects  models in targeted Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS)
 experiment. It is applicable to multiple types of sample preparation, including label-free workflows, workflows that use stable isotope labeled reference proteins and peptides, and workflows that use fractionation. Experimental design of case-control study (patients are not repeatedly measured) or time course study (patients are repeatedly measured) is automatically determined based on proper statistical model.
}
\usage{
groupComparison(contrast.matrix=contrast.matrix, data=data)	
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contrast.matrix}{comparison between conditions of interests.}
  \item{data}{name of the (output of dataProcess function) data set.}
}

\details{

\itemize{
\item{contrast.matrix : comparison of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically. Command levels(QuantData$ProcessedData$GROUP_ORIGINAL) can illustrate the actual order of the levels of conditions.}

}
The underlying model fitting functions are \code{\link[stats]{lm}} and \code{\link[lme4]{lmer}} for the fixed effects model and mixed effects model, respectively.  

The input of this function is the quantitative data from function (\code{\link{dataProcess}}).

}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{

Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\section{Warning }{When a feature is missing completely in a condition or a MS run, a warning message is sent to the console notifying the user of the missing feature. Additional filtering or imputing process is required before model fitting.  
}

\examples{
#Consider quantitative data (i.e. QuantData) from yeast study with ten time points of interests, three biological replicates, and no technical replicates. 
#It is a time-course experiment and we attempt to compare differential abundance between time 1 and 7 in a set of targeted proteins. 
#In this label-based SRM experiment, MSstats uses the fitted model with expanded scope of Biological replication.  

QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

levels(QuantData$ProcessedData$GROUP_ORIGINAL)
comparison<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
row.names(comparison)<-"T7-T1"

# Tests for differentially abundant proteins with models:
# label-based SRM experiment with expanded scope of biological replication.

testResultOneComparison<-groupComparison(contrast.matrix=comparison, data=QuantData)

# table for result
testResultOneComparison$ComparisonResult


}

