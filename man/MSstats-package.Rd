\name{MSstats-package}
\alias{MSstats-package}
\alias{MSstats}
\docType{package}
\title{Tools for protein significance analysis in DDA,SRM and DIA for label-free or label-based proteomic experiments}
\description{A set of tools for protein significance analysis in SRM, DDA and DIA experiments.}
\details{
\tabular{ll}{
Package: \tab MSstats\cr
Version: \tab 2.1.3\cr
Date: \tab 2014-01-05\cr
License: \tab Artistic-2.0 \cr
LazyLoad: \tab yes\cr
}
The package includes four main sections: I. explanatory data analysis (data pre-processing and quality control of MS runs), II. model-based analysis (finding differentially abundant proteins), III. statistical design of future experiments (sample size calculations), and IV. protein quantification (estimation of protein abundance). Section I  contains functions for (1) data pre-processing and quality control of MS runs (see \code{\link{dataProcess}}) and (2) visualizing for explanatory data analysis (see \code{\link{dataProcessPlots}}). Section II contains functions for (1) finding differentially abundant proteins (see \code{\link{groupComparison}}) and (2) visualizing for the testing results (see \code{\link{groupComparisonPlots}}) and for the model-based quality control (see \code{\link{modelBasedQCPlots}}). Section III contains functions for (1) calculating sample size (see \code{\link{designSampleSize}}) and (2) visualizing for the sample size calculations (see \code{\link{designSampleSizePlots}}). Section IV contains functions for (1) per-protein group quantification and patient quantification (see \code{\link{quantification}})

Examples of data or results in MSstats are (1) example of required input data format \code{\link{RawData}}; (2) example of raw data after data pre-processing \code{\link{QuantData}}; (3) results of significance testing of a single comparison \code{\link{testResultOneComparison}}; (4) results of significance testing of multiple comparisons \code{\link{testResultMultiComparisons}}.
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\keyword{ MSstats }


