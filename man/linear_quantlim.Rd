\name{linear_quantlim}
\alias{linear_quantlim}
\title{Calculation of the LOB and LOD with a  linear fit}
\description{
This function calculates the value of the LOB (limit of blank) and LOD (limit of detection) from the (Concentration, Intensity) spiked in data.
The function also returns the values of the linear curve fit that allows it to be plotted. At least 2 blank samples (characterized by Intensity = 0) are required by this function which are used to calculate the background noise. The LOB is defined as the concentration at which the value of the linear fit is equal to the 95\% upper bound of the noise. The LOD is the concentration at which the latter is equal to the 90\% lower bound of the prediction interval (5\% quantile) of the linear fit. 
A weighted linear fit is used with weights for every unique concentration proportional to the inverse of variance between replicates.

}

\usage{
linear_quantlim(datain, alpha = 0.05, Npoints = 100, Nbootstrap = 500)	
}

\arguments{
  \item{datain}{Data frame that contains the input data. The input data frame has to contain the following columns : CONCENTRATION, INTENSITY (both of which are measurements from the spiked in experiment) and NAME which designates the name of the assay (e.g. the name of the peptide or protein)}
  \item{alpha}{Probability level to estimate the LOB/LOD}
  \item{Npoints}{Number of points to use to discretize the concentration line between 0 and the maximum spiked concentration}
  \item{Nbootstrap}{Number of bootstrap samples to use to calculate the prediction interval of the fit. This number has to be increased for very low alpha values or whenever very accurate assay characterization is required.}
}

\details{
\itemize{
\item{datain : Each line of the data frame contains one measurement from the spiked-in experiment. Multiple different INTENSITY values for the same CONCENTRATION are assumed to correspond to different replicates. Blank Samples are characterized by CONCENTRATION = 0.}
}
}

\value{
\itemize{
	\item{Data frame that contains the output of the function. It contains the following columns: i) CONCENTRATION: Concentration values at which the value of the fit is calculated ii) MEAN: The value of the curve fit iii) LOW: The value of the lower bound of the 95\% prediction interval iv) UP: The value of the upper bound of the 95\% prediction interval v) LOB: The value of the LOB (one column with identical values) vi) LOD: The value of the LOD (one column with identical values) vii) SLOPE: Value of the slope of the linear curve fit where only the spikes above LOD are considered viii) INTERCEPT: Value of the intercept of the linear curve fit where only the spikes above LOD are considered ix) NAME: The name of the assay (identical to that provided in the input) x) METHOD which is always set to LINEAR when this function is used. Each line of the data frame corresponds to a unique concentration value at which the value of the fit and prediction interval are evaluated. More unique concentrations values than in the input data frame are used to increase the accuracy of the LOB/D calculations.}
}}


\section{Warning}{The LOB and LOD can only be calculated when more than 2 blank samples are included.
The data should ideally be plotted using the companion function plot_quantlim to ensure that a linear fit is suited to the data.	
}

\author{Cyril Galitzine, Olga Vitek. 

Maintainer: Cyril Galitzine (\email{cyrildgg@gmail.com}), Meena Choi (\email{mnchoi67@gmail.com})
}

\references{
C. Galitzine et al. "Nonlinear regression improves accuracy of characterization of multiplexed mass spectrometric assays" \emph{Mol Cell Proteomics}, doi:10.1074/mcp.RA117.000322, 2018.
}

\examples{

# Consider data from a spiked-in contained in an example dataset
head(SpikeInDataLinear)

\dontrun{
# Call function
linear_quantlim_out <- linear_quantlim(SpikeInDataLinear)
}

}