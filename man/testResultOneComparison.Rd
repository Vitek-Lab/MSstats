\name{testResultOneComparison}
\alias{testResultOneComparison}
\title{Significance testing result of one comparison between two groups}
\description{Consider quantitative data (i.e. \code{\link{QuantData}}) from yeast study with ten time points of interests, three biological replicates, and no technical replicates. It is a time-course experiment, and we attempt to compare differential abundance between time 1 and 7 in a set of targeted proteins. In this label-based SRM experiment, we recommend the fitted model with expanded scope of technical replication (i.e. labeled=TRUE, scopeOfTechReplication="expanded" as default).}

\usage{testResultOneComparison}

\details{The testing result contains variable of Protein, Comparison(Label), log2 fold change(logFC), standard error (SE), T values (Tvalue), degree of freedom (DF), raw p-values (pvalue), adjusted p-values based on Benjamini and Hochberg method to collect multiple testing issue and further control false discovery rate (adj.pvalue).}


\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. Protein significance analysis in selected reaction monitoring (SRM) measurements.  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.}

\examples{
comparison<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
row.names(comparison)<-"T7-T1"

testResultOneComparison<-groupComparison(contrast.matrix=comparison,data=QuantData)

testResultOneComparison$ComparisonResult

}

