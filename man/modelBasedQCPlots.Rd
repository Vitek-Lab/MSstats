\name{modelBasedQCPlots}
\alias{modelBasedQCPlots}
\title{Visualization for model-based quality control in fitting model}
\description{To check the assumption of linear model(summaryMethod="linear" and censoredInt=NULL), modelBasedQCPlots takes the results after fitting models from function (\code{\link{dataProcess}}) as input and automatically generate two types of figures in pdf files as output : (1) normal quantile-quantile plot (specify "QQPlot" in option type) for checking normally distributed errors.; (2) residual plot (specify "ResidualPlot" in option type) for checking constant variance among different features.}

\usage{
modelBasedQCPlots(data,type,
	axis.size=10,dot.size=3,text.size=7,legend.size=7,
	width=10, height=10,
	featureName=TRUE,feature.QQPlot="all",which.Protein="all",address="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{result data.frame called "ModelQC" in output from function dataProcess.}
  \item{type}{choice of visualization. "QQPlots" represents normal quantile-quantile plot for each protein after fitting models. "ResidualPlots" represents a plot of residuals versus fitted values for each protein in the dataset.}
  \item{axis.size}{size of axes labels. Default is 10.}
 \item{dot.size}{size of points in the graph for residual plots and QQ plots. Default is 3.}
  \item{text.size}{size of labeling for feature names only in normal quantile-quantile plots separately for each feature. Default is 7.}
  \item{legend.size}{size of legend for feature names only in residual plots. Default is 7.}
 \item{featureName}{for "ResidualPlot" only, TRUE show feature labeling and FALSE means no feature legend printing.}
 \item{feature.QQPlot}{"all"(Default) means that one normal quantile-quantile plot will be generated with across features of a protein. "byFeature" will generate normal quantile-quantile plots separately for each feature of a protein.}
  \item{width}{width of the saved file. Default is 10.}
  \item{height}{height of the saved file. Default is 10.}
  \item{which.Protein}{Protein list to draw plots. List can be names of Proteins or order numbers of Proteins from levels(data$PROTEIN). Default is "all", which generates all plots for each protein.}
  \item{address}{the name of folder that will store the results. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory. If type="residualPlots" or "QQPlots", "ResidualPlots.pdf" or "QQPlots.plf" will be generated. The command address can help to specify where to store the file as well as how to modify the beginning of the file name. If address=FALSE, plot will be not saved as pdf file but showed in window.}
}

\details{
	Results based on statistical models are accurate as long as the assumptions of the model are met. The model assumes that the measurement errors are normally distributed with mean 0 and constant variance. The assumption of a constant variance can be checked by examining the residuals from the model.
\itemize{
\item{QQPlots : a normal quantile-quantile plot for each protein is generated in order to check whether the errors are well approximated by a normal distribution. If points fall approximately along a straight line, then the assumption is appropriate for that protein. Only large deviations from the line are problematic.}
\item{ResidualPlots : The plots of residuals against predicted(fitted) values. If it shows a random scatter, then the assumption is appropriate.}

}

The input of this function is "ModelQC" in the results from function (\code{\link{dataProcess}}).
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\examples{

QuantData<-dataProcess(SRMRawData, summaryMethod="linear", censoredInt=NULL)

head(QuantData$ModelQC)

# normal quantile-quantile plots
modelBasedQCPlots(data=QuantData$ModelQC,type="QQPlots",address="")

# residual plots
modelBasedQCPlots(data=QuantData$ModelQC,type="ResidualPlots",address="")

}

