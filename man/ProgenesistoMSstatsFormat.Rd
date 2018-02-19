\name{ProgenesistoMSstatsFormat}
\alias{ProgenesistoMSstatsFormat}
\title{Generate MSstats required input format for Progenesis output}
\description{Convert Progenesis output into the required input format for MSstats.}

\usage{
ProgenesistoMSstatsFormat(input,
      annotation,
      useUniquePeptide=TRUE,
      summaryforMultipleRows=max,
      fewMeasurements="remove",
      removeOxidationMpeptides=FALSE,
      removeProtein_with1Peptide=FALSE)	
}

\arguments{
  \item{input}{name of Progenesis output, which is wide-format. 'Accession', 'Sequence', 'Modification', 'Charge' and one column for each run are required.}
  \item{annotation}{name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, Run information. It will be matched with the column name of input for MS runs.}
  \item{useUniquePeptide}{TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
}
  \item{fewMeasurements}{'remove'(default) will remove the features that have 1 or 2 measurements across runs.}
  \item{removeOxidationMpeptides}{TRUE will remove the modified peptides including 'Oxidation (M)' sequence. FALSE is default.}
  \item{removeProtein_with1Peptide}{TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.}

}

\value{data.frame with the required format of MSstats.}

\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}

\examples{
# Please check section 4.4. 
# Suggested workflow with Progenesis output for DDA in MSstats user manual.
# Output of ProgenesistoMSstatsFormat function 
# should have the same 10 columns as an example dataset.

head(DDARawData)
}