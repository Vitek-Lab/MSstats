\name{groupComparison}
\alias{groupComparison}
\title{Finding differentially abundant proteins across conditions in LC-MS, SRM and DIA experiment}
\description{Tests for significant changes in protein abundance across conditions based on a family of linear mixed-effects  models in LC-MS, SRM, DIA experiment. Experimental design of case-control study (patients are not repeatedly measured) or time course study (patients are repeatedly measured) is automatically determined based on proper statistical model. Other choices of model specification include (1) labeling technique: label-based or label-free;  (2) scope of inference: restricted scope (Variable RUN is fixed effect) or expanded scope (Variable RUN is random effect) of technical MS run replication;  (3) interference:  including or excluding additional model interaction to account interference; (4) unequal variance between features: whether the model consider hetergenous variation among intensities between features.
}
\usage{
groupComparison(contrast.matrix=contrast.matrix, data=data, labeled=TRUE, scopeOfBioReplication="restricted", scopeOfTechReplication="expanded", interference=TRUE,featureVar=FALSE,missing.action="nointeraction")	
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contrast.matrix}{comparison between conditions of interests.}
  \item{data}{name of the (processed) data set.}
  \item{labeled}{choice of labeling technique. TRUE(default) represents the label-based study. FALSE represents label-free study.}
  \item{scopeOfBioReplication}{choice of scope of biological replication. "restricted" represents restricted scope of biological replication by specifying subject term as fixed effect in the model. "expanded" (default) represents expanded scope of biological replication by specifying subject term as random effect in the model.}
  \item{scopeOfTechReplication}{choice of scope of technical MS run replication. "restricted" represents restricted scope of technical MS run replication by specifying run term as fixed effect in the model. "expanded" (default) represents expanded scope of technical MS run replication by specifying run term as random effect in the model.}
  \item{interference}{choice of interference data. TRUE(default) means data contain interference transitions and need additional model interaction to address the interference. FALSE means data contain no interference transitions and no need additional model interaction to address the interference.}
  \item{featureVar}{logical variable for whether the model should account for heterogeneous variation among intensities from different features. Default is FALSE, which assume equal variance among intensities from features.}
  \item{missing.action}{specifies the action to take in presence of extreme missing values; must be one of 'nointeraction','impute', or 'remove'. Default is 'nointeraction'.}

}

\details{

\itemize{
\item{contrast.matrix : comparison of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically. Command levels(QuantData$GROUP_ORIGINAL) can illustrate the actual order of the levels of conditions.}
\item{labeled : choices of labeling technique. In label-based study (labeled=TRUE), scopeOfTechReplication, scopeOfBioReplication and interference work as described above. In label-free study (labeled=FALSE), no need to specify scopeOfTechReplication because biological replicates and technical MS runs are confounding. interference works as described above.}
\item{interference : the model can be specified with interaction model terms that reflect interferences in the quantified transitions. }
\item{featureVar : If the unequal variation of error for different peptide features is detected, then a possible solution is to account for the unequal error variation by means of a procedure called iteratively re-weighted least squares. featureVar=TRUE performs an iterative fitting procedure, in which features are weighted inversely proportionaly to the variation in their intensities, so that feature with large variation are given less importance in the estimation of parameters in the model.}
\item{missing.action : When peak intensities from all replicates in a condition are missing for at least one feature, there are three possible actions; (1) remove interaction (missing.action="nointeraction"), which means to assume feature demonstrate no interference across runs, (2) impute with the average minimum intensity across run (missing.action="impute"), or (3) remove the features from the dataset (missing.action="remove")}
}
The underlying model fitting functions are \code{\link[stats]{lm}} and \code{\link[lme4]{lmer}} for the fixed effects model and mixed effects model, respectively.  

The input of this function is the quantitative data from function (\code{\link{dataProcess}}). The example data is \code{\link{QuantData}}.

}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\section{Warning }{When a feature is missing completely in a condition or a MS run, a warning message is sent to the console notifying the user of the missing feature. Additional filtering or imputing process is required before model fitting.  
}

\examples{
#Consider quantitative data (i.e. QuantData) from yeast study with ten time points of interests, three biological replicates, and no technical replicates. 
#It is a time-course experiment and we attempt to compare differential abundance between time 1 and 7 in a set of targeted proteins. 
#In this label-based SRM experiment, we recommend the fitted model with expanded scope of technical replication and restricted scope of biological replication (i.e. labeled=TRUE, scopeOfTechReplication="expanded", scopeOfBioReplication="restricted").   

head(QuantData)

levels(QuantData$GROUP_ORIGINAL)
comparison<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
row.names(comparison)<-"T7-T1"

# Tests for differentially abundant proteins with models:

#(1) label-based SRM experiment with restricted scope of biological replication and expanded scope of technical MS run replication (default) with or without interference
testResultOneComparison<-groupComparison(contrast.matrix=comparison, data=QuantData)
testResultOneComparison$ComparisonResult

testresult2<-groupComparison(contrast.matrix=comparison, data=QuantData, interference=FALSE)
testresult2$ComparisonResult

#(2) label-based SRM experiment with restricted scope of technical MS run replication and restricted scope of biological replication
testresult3<-groupComparison(contrast.matrix=comparison, data=QuantData, scopeOfTechReplication="restricted")
testresult3$ComparisonResult

#(3) label-based SRM experiment with expanded scope of technical MS run replication and expanded scope of biological replication 
testresult4<-groupComparison(contrast.matrix=comparison, data=QuantData, scopeOfBioReplication="expanded")
testresult4$ComparisonResult

#(4) label-free SRM experiment with expanded scope of biological replication and interference
testresult5<-groupComparison(contrast.matrix=comparison, data=QuantData, labeled=FALSE)
testresult5$ComparisonResult

}

