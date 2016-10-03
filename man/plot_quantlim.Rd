\name{plot_quantlim}
\alias{plot_quantlim}
\title{Plot of the curve used to calculate LOB and LOD}
\description{
This function allows to plot the curve fit that is used to calculate the LOB and LOD with functions nonlinear_quantlim() and linear_quantlim(). The function outputs for each calibration curve, two pdf files each containg one plot. On the first, designated by *_overall.pdf, the entire concentration range is plotted. On the second plot, designated by *_zoom.pdf,, the concentration range between 0 and xlim_plot (if specified in the argument of the function) is plotted.
When no xlim_plot value is specified, the region close to LOB and LOD is automatically plotted.
}

\usage{
plot_quantlim(spikeindata, quantlim_out, alpha, dir_output, xlim_plot)
}

\arguments{
  \item{spikeindata}{Data frame that contains the experimental spiked in data. This data frame should be identical to that used as input by function functions nonlinear_quantlim() or linear_quantlim().  The data frame has to contain the following columns : CONCENTRATION, INTENSITY (both of which are measurements from the spiked in experiment) and NAME which designates the name of the assay (e.g. the name of the peptide or protein)}
  \item{quantlim_out}{Data frame that was output by functions nonlinear_quantlim() or linear_quantlim(). It has to contain at least the following columns: i) CONCENTRATION: Concentration values at which the value of the fit is calculated ii) MEAN: The value of the curve fit iii) LOW: The value of the lower bound of the 95\% prediction interval iv) UP: The value of the upper bound of the 95\% prediction interval v) LOB: The value of the LOB (one column with identical values) vi) LOD: The value of the LOD (one column with identical values) vii) NAME: The name of the assay (identical to that provided in the input) viii) METHOD which is LINEAR or NONLINEAR}
  \item{alpha}{Probability level to estimate the LOB/LOD}
\item{dir_output}{String containg the path of the directly where the pdf files of the plots should be output.}
\item{xlim_plot}{Optional argument containing the maximum xaxis value of the zoom plot. When no value is specified, a suitable value close to LOD is automatically chosen.}
}

\details{
\itemize{
\item{spikeindata : Each line of the data frame contains one measurement from the spiked-in experiment. Multiple different INTENSITY values for the same CONCENTRATION are assumed to correspond to different replicates. Blank Samples are characterized by CONCENTRATION = 0.}

\item{output : Each line of the data frame corresponds to a unique concentration value at which the value of the fit and prediction interval are evaluated. More unique concentrations values than in the input data frame are used to increase the accuracy of the LOB/D calculations.}
}
}

\section{Warning}{
This plotting function should ideally be used every time  nonlinear_quantlim() or linear_quantlim() are called to visually ensure that the fits and data are accurate.
}

\author{Cyril Galitzine, Olga Vitek. 

Maintainer: Cyril Galitzine (\email{cyrildgg@gmail.com}), Meena Choi (\email{mnchoi67@gmail.com})
}
\references{

C. Galitzine et al.  TBD  2016

}

\examples{
# Consider data from a spiked-in contained in an example dataset. This dataset contains 
# a significant threshold at low concentrations that is not well captured by a linear fit.

head(SpikeInDataNonLinear)

#Call function
nonlinear_quantlim_out <- nonlinear_quantlim(SpikeInDataNonLinear, alpha = 0.05)

#Get values of LOB/LOD
nonlinear_quantlim_out$LOB[1]
nonlinear_quantlim_out$LOD[1]

plot_quantlim(spikeindata = SpikeInDataLinear, quantlim_out  = nonlinear_quantlim_out,
dir_output =  getwd(), alpha = 0.05)

}