\name{SpectronauttoMSstatsFormat}
\alias{SpectronauttoMSstatsFormat}
\title{Generate MSstats required input format for Spectronaut output}
\description{Convert Spectronaut output into the required input format for MSstats.}

\usage{
SpectronauttoMSstatsFormat(input,
      intensity = 'PeakArea',
      filter_with_Qvalue = TRUE,
      qvalue_cutoff = 0.01,
      useUniquePeptide = TRUE,
      fewMeasurements="remove",
      removeProtein_with1Feature = FALSE,
      summaryforMultipleRows=max)	
}

\arguments{
  \item{input}{name of Spectronaut output, which is long-format. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity, F.ExcludedFromQuantification are required. Rows with F.ExcludedFromQuantification=True will be removed.}
  \item{intensity}{'PeakArea'(default) uses not normalized peak area. 'NormalizedPeakArea' uses peak area normalized by Spectronaut.}
  \item{filter_with_Qvalue}{TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.}
  \item{qvalue_cutoff}{Cutoff for EG.Qvalue. default is 0.01.}
  \item{useUniquePeptide}{TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.}
  \item{fewMeasurements}{'remove'(default) will remove the features that have 1 or 2 measurements across runs.}
  \item{removeProtein_with1Feature}{TRUE will remove the proteins which have only 1 feature, which is the combination of peptide, precursor charge, fragment and charge. FALSE is default.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.}
}


\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}

\examples{
# Please check section 5.2. 
# Suggested workflow with Spectronaut output for DIA in MSstats user manual.
# Output of SpectronauttoMSstatsFormat function 
# should have the same 10 columns as an example dataset.

head(DDARawData)
}