\name{MaxQtoMSstatsFormat}
\alias{MaxQtoMSstatsFormat}
\title{Generate MSstats required input format for MaxQuant output}
\description{Convert MaxQuant output into the required input format for MSstats.}

\usage{
MaxQtoMSstatsFormat(evidence, annotation,proteinGroups,
 proteinID="Proteins", useUniquePeptide=TRUE, 
 summaryforMultipleRows=max, 
 fewMeasurements="remove", removeMpeptides=TRUE)	
}

\arguments{
  \item{evidence}{name of 'evidence.txt' data, which includes feature-level data.}
  \item{annotation}{name of 'annotation.txt' data which includes Raw.file, Condition, BioReplicate, Run, IsotopeLabelType information.}
  \item{proteinGroups}{name of 'proteinGroups.txt' data. It needs to matching protein group ID. If proteinGroups=NULL, use 'Proteins' column in 'evidence.txt'.}
  \item{proteinID}{'Proteins'(default) or 'proteinGroups' in 'proteinGroup.txt' for Protein ID.}
  \item{useUniquePeptide}{TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain fun, use highest or sum of all.
}
  \item{fewMeasurements}{'remove'(default) will remove the features that have 1 or 2 measurements across runs. It can affect 'dataProcess' function with unequal variance between features option.}
  \item{removeMpeptides}{TRUE(default) will remove the peptides including 'M' sequence.}
}


\section{Warning}{MSstats does not support for metabolic labeling or iTRAQ experiments.
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
