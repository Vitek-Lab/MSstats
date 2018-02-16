\name{SkylinetoMSstatsFormat}
\alias{SkylinetoMSstatsFormat}
\title{Generate MSstats required input format for Skyline output}
\description{Preprocess MSstats input report from Skyline and convert into the required input format for MSstats.}

\usage{
SkylinetoMSstatsFormat(input,
      annotation = NULL,
      removeiRT = TRUE,
      filter_with_Qvalue = TRUE,
      qvalue_cutoff = 0.01,
      useUniquePeptide = TRUE,
      fewMeasurements="remove",
      removeOxidationMpeptides = FALSE,
      removeProtein_with1Feature = FALSE)	
}

\arguments{
  \item{input}{name of MSstats input report from Skyline, which includes feature-level data.}
  \item{annotation}{name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Skyline, use annotation=NULL (default). It will use the annotation information from input.}
  \item{removeiRT}{TRUE(default) will remove the proteins or peptides which are labeld 'iRT' in 'StandardType' column. FALSE will keep them.}
  \item{filter_with_Qvalue}{TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in DetectionQValue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.}
  \item{qvalue_cutoff}{Cutoff for DetectionQValue. default is 0.01.}
  \item{useUniquePeptide}{TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.}
  \item{fewMeasurements}{'remove' (default) will remove the features that have 1 or 2 measurements across runs.}
  \item{removeOxidationMpeptides}{TRUE will remove the peptides including 'oxidation (M)' in modification. FALSE is default.}
  \item{removeProtein_with1Feature}{TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.}
}

\value{data.frame with the required format of MSstats.}

\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}

\examples{
# Please check section 4.2. Suggested workflow with Skyline output for DDA in MSstats user manual.
# Output of SkylinetoMSstatsFormat function should have the same 10 columns as an example dataset.

head(DDARawData)
}