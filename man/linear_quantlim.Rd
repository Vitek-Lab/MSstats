\name{linear_quantlim}
\alias{linear_quantlim}
\title{Calculation of the LOD and LOQ with a  linear fit}
\description{
This function calculates the value of the LOD (limit of detection) and LOQ (limit of quantification) from the (Concentration, Intensity) spiked in data.
The function also returns the values of the linear curve fit that allows it to be plotted. At least 2 blank samples (characterized by Intensity = 0) are required by this function which are used to calculate the background noise. The LOD is defined as the concentration at which the value of the linear fit is equal to the 95\% upper bound of the noise. The LOQ is the concentration at which the latter is equal to the 95\% lower bound of the prediction interval of the linear fit. 
A weighted linear fit is used with weights for every unique concentration proportional to the inverse of variance between replicates.

}

\usage{
linear_quantlim(datain)	
}

\arguments{
  \item{datain}{Data frame that contains the input data. The input data frame has to contain the following columns : CONCENTRATION, INTENSITY (both of which are measurements from the spiked in experiment) and NAME which designates the name of the assay (e.g. the name of the peptide or protein)}
}

\details{
\itemize{
\item{datain : Each line of the data frame contains one measurement from the spiked-in experiment. Multiple different INTENSITY values for the same CONCENTRATION are assumed to correspond to different replicates. Blank Samples are characterized by CONCENTRATION = 0.}

\item{output : Data frame that contains the output of the function. It contains the following columns: i) CONCENTRATION: Concentration values at which the value of the fit is calculated ii) MEAN: The value of the curve fit iii) LOW: The value of the lower bound of the 95\% prediction interval iv) UP: The value of the upper bound of the 95\% prediction interval v) LOD: The value of the LOD (one column with identical values) vi) LOQ: The value of the LOQ (one column with identical values) vii) NAME: The name of the assay (identical to that provided in the input) viii) METHOD which is always set to LINEAR when this function is used.Each line of the data frame corresponds to a unique concentration value at which the value of the fit and prediction interval are evaluated. More unique concentrations values than in the input data frame are used to increase the accuracy of the LOD/Q calculations.}
}
}

\section{Warning}{The LOD and LOQ can only be calculated when more than 2 blank samples are included.
The data should ideally be plotted using the companion function plot_quantlim to ensure that a linear fit is suited to the data.	
}

\author{Cyril Galitzine, Olga Vitek. 

Maintainer: Cyril Galitzine (\email{cyrildgg@gmail.com}), Meena Choi (\email{mnchoi67@gmail.com})
}
\references{

C. Galitzine et al.  TBD  2016

}

\examples{
# Consider data from a spiked-in contained in an example dataset
head(SpikeInDataLinear)

#Call function
linear_quantlim_out <- linear_quantlim(SpikeInDataLinear)

#Get values of LOD/LOQ
linear_quantlim_out$LOD[1]
linear_quantlim_out$LOQ[1]



}