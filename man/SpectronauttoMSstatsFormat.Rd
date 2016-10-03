\name{SpectronauttoMSstatsFormat}
\alias{SpectronauttoMSstatsFormat}
\title{Generate MSstats required input format for Spectronaut output}
\description{Convert Spectronaut output into the required input format for MSstats.}

\usage{
SpectronauttoMSstatsFormat(input,
      summaryforMultipleRows=max)	
}

\arguments{
  \item{input}{name of Spectronaut output, which is long-format. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity, F.ExcludedFromQuantification are required. Rows with F.ExcludedFromQuantification=True will be removed.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.}
}


\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}
