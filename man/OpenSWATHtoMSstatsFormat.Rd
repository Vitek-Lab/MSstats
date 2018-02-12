\name{OpenSWATHtoMSstatsFormat}
\alias{OpenSWATHtoMSstatsFormat}
\title{Generate MSstats required input format for OpenSWATH output}
\description{Preprocess MSstats input report from OpenSWATH and convert into the required input format for MSstats.}

\usage{
OpenSWATHtoMSstatsFormat(input,
    annotation = NULL,
    filter_with_mscore = TRUE,
    mscore_cutoff = 0.01,
    useUniquePeptide = TRUE,
    fewMeasurements="remove",
    removeProtein_with1Feature = FALSE,
    summaryforMultipleRows=max)	
}

\arguments{
  \item{input}{name of MSstats input report from OpenSWATH, which includes feature-level data.}
  \item{annotation}{name of 'annotation.txt' data which includes Condition, BioReplicate, Run. Run should be the same as filename.}
  \item{filter_with_mscore}{TRUE(default) will filter out the features that have greater than mscore_cutoff in m_score column. Those features will be removed.}
  \item{mscore_cutoff}{Cutoff for m_score. default is 0.01.}
  \item{useUniquePeptide}{TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.}
  \item{fewMeasurements}{'remove' (default) will remove the features that have 1 or 2 measurements across runs.}
  \item{removeProtein_with1Feature}{TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.}
}


\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}

\examples{
# Example will be ready in next version.
}