\name{DIAUmpiretoMSstatsFormat}
\alias{DIAUmpiretoMSstatsFormat}
\title{Generate MSstats required input format for DIA-Umpire output}
\description{Convert DIA-Umpire output into the required input format for MSstats.}

\usage{
DIAUmpiretoMSstatsFormat(raw.frag, raw.pep, raw.pro,
    annotation,
    useSelectedFrag = TRUE,
    useSelectedPep = TRUE,
    fewMeasurements="remove",
    removeProtein_with1Feature = FALSE,
    summaryforMultipleRows=max)	
}

\arguments{
  \item{raw.frag}{name of FragSummary_date.xls data, which includes feature-level data.}
  \item{raw.pep}{name of PeptideSummary_date.xls data, which includes selected fragments information.}
  \item{raw.pro}{name of ProteinSummary_date.xls data, which includes selected peptides information.}
  \item{annotation}{name of annotation data which includes Condition, BioReplicate, Run information.}
  \item{useSelectedFrag}{TRUE will use the selected fragment for each peptide. 'Selected_fragments' column is required.}
  \item{useSelectedPep}{TRUE will use the selected peptide for each protein. 'Selected_peptides' column is required.}
  \item{fewMeasurements}{'remove'(default) will remove the features that have 1 or 2 measurements across runs.}
  \item{removeProtein_with1Feature}{TRUE will remove the proteins which have only 1 feature, which is the combination of peptide, precursor charge, fragment and charge. FALSE is default.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.}
}

\value{data.frame with the required format of MSstats.}

\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}

\examples{
# Manual will be updated.
# Output of DIAUmpiretoMSstatsFormat function 
# should have the same 10 columns as an example dataset.

head(DDARawData)
}