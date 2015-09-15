\name{dataProcessPlots}
\alias{dataProcessPlots}
\title{Visualization for explanatory data analysis}
\description{To illustrate the quantitative data after data-preprocessing and quality control of MS runs, dataProcessPlots takes the quantitative data from function (\code{\link{dataProcess}}) as input and automatically generate three types of figures in pdf files as output : (1) profile plot (specify "ProfilePlot" in option type), to identify the potential sources of variation for each protein; (2) quality control plot (specify "QCPlot" in option type), to evaluate the systematic bias between MS runs; (3) mean plot for conditions (specify "ConditionPlot" in option type), to illustrate mean and variability of each condition per protein. }
\usage{
dataProcessPlots(data=data,type=type,featureName="Transition",
	ylimUp=FALSE,ylimDown=FALSE,scale=FALSE,interval="CI",
	x.axis.size=10,y.axis.size=10,
	text.size=4,text.angle=0,legend.size=7,dot.size.profile=2,dot.size.condition=3,
	width=10, height=10, which.Protein="all",
	originalPlot=TRUE, summaryPlot=TRUE, address="")
}

\arguments{
  \item{data}{name of the (output of dataProcess function) data set.}
  \item{type}{choice of visualization. "ProfilePlot" represents profile plot of log intensities across MS runs. "QCPlot" represents quality control plot of log intensities across MS runs. "ConditionPlot" represents mean plot of log ratios (Light/Heavy) across conditions.}
  \item{featureName}{for "ProfilePlot" only, "Transition" (default) means printing feature legend in transition-level; "Peptide" means printing feature legend in peptide-level; "NA" means no feature legend printing.} 
  \item{ylimUp}{upper limit for y-axis in the log scale. FALSE(Default) for Profile Plot and QC Plot is 30. FALSE(Default) for Condition Plot is maximum of log ratio + SD or CI.}
  \item{ylimDown}{lower limit for y-axis in the log scale. FALSE(Default) for Profile Plot and QC Plot is 0. FALSE(Default) for Condition Plot is minumum of log ratio - SD or CI.}
  \item{scale}{for "ConditionPlot" only, FALSE(default) means each conditional level is not scaled at x-axis according to its actual value (equal space at x-axis). TRUE means each conditional level is scaled at x-axis according to its actual value (unequal space at x-axis).}
  \item{interval}{for "ConditionPlot" only, "CI"(default) uses confidence interval with 0.95 significant level for the width of error bar. "SD" uses standard deviation for the width of error bar.}
  \item{x.axis.size}{size of x-axis labeling for "Run" in Profile Plot and QC Plot, and "Condition" in Condition Plot. Default is 10.}
  \item{y.axis.size}{size of y-axis labels. Default is 10.}
  \item{text.size}{size of labels represented each condition at the top of graph in Profile Plot and QC plot. Default is 4.}
    \item{text.angle}{angle of labels represented each condition at the top of graph in Profile Plot and QC plot or x-axis labeling in Condition plot. Default is 0.}
  \item{legend.size}{size of feature legend (transition-level or peptide-level) above graph in Profile Plot. Default is 7.}
  \item{dot.size.profile}{size of dots in profile plot. Default is 2.}
  \item{dot.size.condition}{size of dots in condition plot. Default is 3.}
  \item{width}{width of the saved file. Default is 10.}
  \item{height}{height of the saved file. Default is 10.}
  \item{which.Protein}{Protein list to draw plots. List can be names of Proteins or order numbers of Proteins from levels(data$PROTEIN). Default is "all", which generates all plots for each protein.}
  \item{originalPlot}{TRUE(default) draws original profile plots.}
  \item{summaryPlot}{TRUE(default) draws profile plots with summarization for run levels.}
  \item{address}{the name of folder that will store the results. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory. An output pdf file is automatically created with the default name of "ProfilePlot.pdf" or "QCplot.pdf" or "ComparisonPlot.pdf". The command address can help to specify where to store the file as well as how to modify the beginning of the file name. If address=FALSE, plot will be not saved as pdf file but showed in window.}
}

\details{

\itemize{
\item{Profile Plot : identify the potential sources of variation of each protein. X-axis is run. Y-axis is log-intensities of transitions. Reference/endogenous signals are in the left/right panel. Line colors indicate peptides and line types indicate transitions.}
\item{QC Plot : illustrate the systematic bias between MS runs. After normalization, the reference signals for all proteins should be stable across MS runs. X-axis is run. Y-axis is log-intensities of transition. Reference/endogenous signals are in the left/right panel. The pdf file contains (1) QC plot for all proteins and (2) QC plots for each protein separately.}
\item{Condition Plot : illustrate the systematic difference between conditions. X-axis is condition. Y-axis is log ratio of endogenous over reference. For label-free, Y-axis is log intensity of endogenous. If scale is TRUE, the levels of conditions is scaled according to its actual values at x-axis. Red points indicate the mean of log ratio for each condition. If interval is "CI", blue error bars indicate the confidence interval with 0.95 significant level for each condition. If interval is "SD", blue error bars indicate the standard deviation for each condition.The interval is not related with model-based analysis}
}
The input of this function is the quantitative data from function (\code{\link{dataProcess}}). 
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
	
Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements."  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\examples{
#Consider quantitative data (i.e. QuantData) from a yeast study with ten time points of interests, three biological replicates, and no technical replicates which is a time-course experiment. 
#The goal is to provide pre-analysis visualization by automatically generate two types of figures in two separate pdf files. 
#Protein IDHC (gene name IDP2) is differentially expressed in time point 1 and time point 7, whereas, Protein PMG2 (gene name GPM2) is not.

QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

# Profile plot

dataProcessPlots(data=QuantData,type="ProfilePlot")

# Quality control plot 

dataProcessPlots(data=QuantData,type="QCPlot")	

# Quantification plot for conditions

dataProcessPlots(data=QuantData,type="ConditionPlot")

}

