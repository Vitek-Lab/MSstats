\name{testResultMultiComparisons}
\alias{testResultMultiComparisons}
\title{Significance testing result of a set of comparisons between two groups}
\description{Consider quantitative data (i.e. \code{\link{QuantData}}) from yeast study with ten time points of interests, three biological replicates, and no technical replicates. It is a time-course experiment. In this label-based SRM experiment, we recommend the fitted model with expanded scope of technical replication (i.e. labeled=TRUE, scopeOfTechReplication="expanded" as default). This is a testing result of multiple comparisons of \code{\link{QuantData}} based on the intensity-based linear model with expanded scope. The comparison is time 1 vs time 3 (T3-T1), time 1 vs time 7 (T7-T1), time 1 vs time 9 (T9-T1).}

\usage{testResultMultiComparisons}

\details{The testing result contains variable of Protein, Comparison(Label), log2 fold change(logFC), standard error (SE), T values (Tvalue), degree of freedom (DF), raw p-values (pvalue), adjusted p-values based on Benjamini and Hochberg method to collect multiple testing issue and further control false discovery rate (adj.pvalue). There are multiple lines for the same proteins.}


\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. Protein significance analysis in selected reaction monitoring (SRM) measurements.  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.}

\examples{

comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
comparison<-rbind(comparison1,comparison2, comparison3)
row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")

testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=QuantData)

testResultMultiComparisons$ComparisonResult
}

