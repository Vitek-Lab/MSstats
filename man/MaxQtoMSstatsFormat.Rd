\name{MaxQtoMSstatsFormat}
\alias{MaxQtoMSstatsFormat}
\title{Generate MSstats required input format for MaxQuant output}
\description{Convert MaxQuant output into the required input format for MSstats.}

\usage{
MaxQtoMSstatsFormat(evidence, 
      annotation,
      proteinGroups,
      proteinID="Proteins", 
      useUniquePeptide=TRUE, 
      summaryforMultipleRows=max, 
      fewMeasurements="remove", 
      removeMpeptides=FALSE,
      removeOxidationMpeptides=FALSE,
      removeProtein_with1Peptide=FALSE)	
}

\arguments{
  \item{evidence}{name of 'evidence.txt' data, which includes feature-level data.}
  \item{annotation}{name of 'annotation.txt' data which includes Raw.file, Condition, BioReplicate, Run, IsotopeLabelType information.}
  \item{proteinGroups}{name of 'proteinGroups.txt' data. It needs to matching protein group ID. If proteinGroups=NULL, use 'Proteins' column in 'evidence.txt'.}
  \item{proteinID}{'Proteins'(default) or 'Leading.razor.protein' for Protein ID.}
  \item{useUniquePeptide}{TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of all.
}
  \item{fewMeasurements}{'remove'(default) will remove the features that have 1 or 2 measurements across runs. 'keep' will keep all features.}
  \item{removeMpeptides}{TRUE will remove the peptides including 'M' sequence. FALSE is default.}
  \item{removeOxidationMpeptides}{TRUE will remove the peptides including 'oxidation (M)' in modification. FALSE is default.}
  \item{removeProtein_with1Peptide}{TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.}

}

\value{data.frame with the required format of MSstats.}

\section{Warning}{MSstats does not support for metabolic labeling or iTRAQ experiments.
}

\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}

\examples{
# Please check section 4.3. Suggested workflow with MaxQuant output for DDA in MSstats user manual.
# Output of MaxQtoMSstatsFormat function should have the same 10 columns as an example dataset.

head(DDARawData)
}

