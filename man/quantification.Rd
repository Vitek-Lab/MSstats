\name{quantification}
\alias{quantification}
\title{Protein sample quantification or group quantification}
\description{Model-based quantification for each condition or for each biological samples per protein in a targeted Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment. Quantification takes the processed data set by \code{\link{dataProcess}} as input and automatically generate the quantification results (data.frame) with long or matrix format.
}
\usage{
quantification(data,type="Sample",format="matrix", scopeOfTechReplication="restricted", scopeOfBioReplication="restricted", interference=TRUE, missing.action="nointeraction",equalFeatureVar=TRUE)
}

\arguments{
  \item{data}{name of the (processed) data set.}
  \item{type}{choice of quantification. "Sample" or "Group" for protein sample quantification or group quantification.}
  \item{format}{choice of returned format. "long" for long format which has the columns named Protein, Condition, LonIntensities (and BioReplicate if it is subject quantification), NumFeature for number of transitions for a protein, and NumPeaks for number of observed peak intensities for a protein. "matrix" for data matrix format which has the rows for Protein and the columns, which are Groups(or Conditions) for group quantification or the combinations of BioReplicate and Condition (labeled by "BioReplicate"_"Condition") for sample quantification. Default is "matrix" }
    \item{scopeOfBioReplication}{choice of scope of biological replication. "restricted" (default) represents restricted scope of biological replication by specifying subject term as fixed effect in the model. "expanded" represents expanded scope of biological replication by specifying subject term as random effect in the model.}
  \item{scopeOfTechReplication}{choice of scope of technical MS run replication. "restricted"(default) represents restricted scope of technical MS run replication by specifying run term as fixed effect in the model. "expanded" represents expanded scope of technical MS run replication by specifying run term as random effect in the model.}
  \item{interference}{choice of interference data. TRUE(default) means data contain interference transitions and need additional model interaction to address the interference. FALSE means data contain no interference transitions and no need additional model interaction to address the interference.}
  \item{missing.action}{specifies the action to take in presence of extreme missing values; must be one of 'nointeraction','impute', or 'remove'. Default is 'nointeraction'.}
  \item{equalFeatureVar}{logical variable for whether the model should account for heterogeneous variation among intensities from different features. Default is TRUE, which assume equal variance among intensities from features. FALSE means that we cannot assume equal variance among intensities from features, then we will account for heterogeneous variation from different features.}
%%  \item{address}{the name of folder that will store the results. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory.An output csv file is automatically created with the default name of "XXX.csv". The command address can help to specify where to store the file as well as how to modify the beginning of the file name.}
}

\details{

\itemize{
\item{Sample quantification : model-based individual biological sample quantification for each protein. The label of each biological sample is a combination of the corresponding group and the sample ID. The same model with groupComparison will be used. However, if there is only one transition in a certain protein, the estimate of variation is NA. Therefore, the result may be unreliable.}

\item{Group quantification : model-based quantification for individual group or individual condition per protein. The same model with groupComparison will be used. The quantification for reference is the average among all reference intensities.}

\item{The quantification for endogenous samples is based on the log-intensities of model-based averaging of all endogenous transitions within a specific sample. The quantification for reference sample is based on  the log-intensities of the model-based averaging among all reference transitions. The quantification of log-ratios of specific endogenous sample over reference sample can be obtained by the quantification of that endogenous sample minus the quantification of the reference sample.}

%%\item{results : whether format="long" or "matrix", both files, "Group or SampleQuantification_longformat.csv" and "Group or SampleQuantification_dataMatrix.csv" will be stored in the assigned folder}

The input of this function is the quantitative data from function (\code{\link{dataProcess}}).
}
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
#Consider quantitative data (i.e. QuantData) from a yeast study with ten time points of interests, three biological replicates, and no technical replicates which is a time-course experiment. 
#Sample quantification shows model-based estimation of protein abundance in each biological replicate within each time point.
#Group quantification shows model-based estimation of protein abundance in each time point.

QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

# Sample quantification

sampleQuant<-quantification(QuantData$ProcessedData)
head(sampleQuant)

# Group quantification

groupQuant<-quantification(QuantData$ProcessedData, type="Group")
head(groupQuant)
}

