\name{SpikeInDataNonLinear}
\alias{SpikeInDataNonLinear}
\title{Example dataset from an MRM spike-in experiment with a nonlinear behavior}
\description{
This dataset is part of the CPTAC 7, study 3 (Addona et al., 2009). It corresponds to the spike-in data for peptide ESDTSYVSLK at site 19. This particular data was chosen because of the concentration threshold that is present at low concentrations that warrant the use of a nonlinear method. The data is composed of 4 replicates at 10 different concentrations (including a blank sample with concentration 0).
}

\usage{SpikeInDataNonLinear}

\format{data.frame}

\details{The intensity reported is the sum of the intensity of all the different fragments of the peptide. Only the peptide being spiked (light peptide) is contained in the example data set. The intensity was normalized using the corresponding heavy peptide in log space such that intensity of the heavy remains constant for all concentrations and all replicates. The intensity was rescaled following the method described in Addona et al., 2009. The concentration and Intensity are both in units of fmol/uL.
}

\value{data.frame as described in details.}

\author{Cyril Galitzine, Olga Vitek. 

Maintainer: Cyril Galitzine (\email{cyrildgg@gmail.com}), Meena Choi (\email{mnchoi67@gmail.com})
}

\references{
T.A. Addonna et al. "Multi-site assessment of the precision and reproducibility of multiple reaction monitoring-based measurements of proteins in plasma."  \emph{Nat Biotechnol.} 2009 Jul;27(7):633-41
}

\examples{
head(SpikeInDataNonLinear)
}

