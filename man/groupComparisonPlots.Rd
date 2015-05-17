\name{groupComparisonPlots}
\alias{groupComparisonPlots}
\title{Visualization for model-based analysis and summarizing differentially abundant proteins}
\description{To summarize the results of log-fold changes and adjusted p-values for differentially abundant proteins, groupComparisonPlots takes testing results from function (\code{\link{groupComparison}}) as input and automatically generate three types of figures in pdf files as output : (1) volcano plot (specify "VolcanoPlot" in option type) for each comparison separately; (2) heatmap (specify "Heatmap" in option type) for multiple comparisons ; (3) comparison plot (specify "ComparisonPlot" in option type) for multiple comparisons per protein.}

\usage{
groupComparisonPlots(data=data,type=type,sig=0.05,FCcutoff=FALSE,
	logBase.pvalue=10,ylimUp=FALSE,ylimDown=FALSE,xlimUp=FALSE,
	x.axis.size=10,y.axis.size=10,dot.size=3,text.size=4,legend.size=7,
	ProteinName=TRUE,ProteinNameLoc=1, numProtein=100, clustering="both", 
	width=10, height=10, which.Comparison="all", address="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{'ComparisonResult' in testing output from function groupComparison.}
  \item{type}{choice of visualization. "VolcanoPlot" represents volcano plot of log fold changes and adjusted p-values for each comparison separately. "Heatmap" represents heatmap of adjusted p-values for multiple comparisons. "ComparisonPlot" represents comparison plot of log fold changes for multiple comparisons per protein.}
  \item{sig}{FDR cutoff for the adjusted p-values in heatmap and volcano plot. level of significance for comparison plot. 100(1-sig)\% confidence interval will be drawn.  sig=0.05 is default.}
  \item{FCcutoff}{for volcano plot or heatmap, whether involve fold change cutoff or not. FALSE (default) means no fold change cutoff is applied for significance analysis. FCcutoff = specific value means specific fold change cutoff is applied.}
  \item{logBase.pvalue}{for volcano plot or heatmap, (-) logarithm transformation of adjusted p-value with base 2 or 10(default).}
  \item{ylimUp}{for all three plots, upper limit for y-axis. FALSE (default) for volcano plot/heatmap use maximum of -log2 (adjusted p-value) or -log10 (adjusted p-value). FALSE (default) for comparison plot uses maximum of log-fold change + CI.}
  \item{ylimDown}{for all three plots, lower limit for y-axis. FALSE (default) for volcano plot/heatmap use minimum of -log2 (adjusted p-value) or -log10 (adjusted p-value). FALSE (default) for comparison plot uses minimum of log-fold change - CI.}
  \item{xlimUp}{for Volcano plot, the limit for x-axis. FALSE (default) for  use maximum for absolute value of log-fold change or 3 as default if maximum for absolute value of log-fold change is less than 3.}
  \item{x.axis.size}{size of axes labels, e.g. name of the comparisons in heatmap, and in comparison plot. Default is 10.}
  \item{y.axis.size}{size of axes labels, e.g. name of targeted proteins in heatmap. Default is 10.}
  \item{dot.size}{size of dots in volcano plot and comparison plot. Default is 3.}
 \item{text.size}{size of ProteinName label in the graph for Volcano Plot. Default is 4.}
  \item{legend.size}{size of legend for color at the bottom of volcano plot.  Default is 7.}
  \item{ProteinName}{for volcano plot only, whether display protein names or not. TRUE (default) means protein names are displayed next to the testing results. FALSE means no protein names are displayed.}
  \item{ProteinNameLoc}{for volcano plot only, assign the distance between the point and the displayed protein name. Default is 1.}
   \item{numProtein}{The number of proteins which will be presented in each heatmap. Default is 100. Maximum possible number of protein for one heatmap is 180. }
   \item{clustering}{Determines how to order proteins and comparisons. Hierarchical cluster analysis with Ward method(minimum variance) is performed. 'protein' means that protein dendrogram is computed and reordered based on protein means (the order of row is changed). 'comparison' means comparison dendrogram is computed and reordered based on comparison means (the order of comparison is changed). 'both' means to reorder both protein and comparison. Default is 'protein'.}
  \item{width}{width of the saved file. Default is 10.}
  \item{height}{height of the saved file. Default is 10.}
   \item{which.Comparison}{list of comparisons to draw plots. List can be labels of comparisons or order numbers of comparisons from levels(data$Label), such as levels(testResultMultiComparisons$ComparisonResult$Label). Default is "all", which generates all plots for each protein.}
  \item{address}{the name of folder that will store the results. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory. An output pdf file is automatically created with the default name of "VolcanoPlot.pdf" or "Heatmap.pdf" or "ComparisonPlot.pdf". The command address can help to specify where to store the file as well as how to modify the beginning of the file name. If address=FALSE, plot will be not saved as pdf file but showed in window.}
}

\details{
\itemize{
\item{Volcano plot : illustrate actual log-fold changes and adjusted p-values for each comparison separately with all proteins. The x-axis is the log fold change. The base of logarithm transformation is the same as specified in "logTrans" from \code{\link{dataProcess}}. The y-axis is the negative log2 or log10 adjusted p-values. The horizontal dashed line represents the FDR cutoff. The points below the FDR cutoff line are non-significantly abundant proteins (colored in black). The points above the FDR cutoff line are significantly abundant proteins (colored in red/blue for up-/down-regulated). If fold change cutoff is specified (FCcutoff = specific value), the points above the FDR cutoff line but within the FC cutoff line are non-significantly abundant proteins (colored in black)/}

\item{Heatmap : illustrate up-/down-regulated proteins for multiple comparisons with all proteins. Each column represents each comparison of interest. Each row represents each protein. Color red/blue represents proteins in that specific comparison are significantly up-regulated/down-regulated proteins with FDR cutoff and/or FC cutoff. The color scheme shows the evidences of significance. The darker color it is, the stronger evidence of significance it has. Color gold represents proteins are not significantly different in abundance.}

\item{Comparison plot : illustrate log-fold change and its variation of multiple comparisons for single protein. X-axis is comparison of interest. Y-axis is the log fold change. The red points are the estimated log fold change from the model. The blue error bars are the confidence interval with 0.95 significant level for log fold change. This interval is only based on the standard error, which is estimated from the model. }
}

The input of this function is "ComparisonResult" in the testing results from function (\code{\link{groupComparison}}). 
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

QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

## based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
comparison<-rbind(comparison1,comparison2, comparison3)
row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")

testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=QuantData)

testResultMultiComparisons$ComparisonResult

# Volcano plot with FDR cutoff = 0.05 and no FC cutoff
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,type="VolcanoPlot",logBase.pvalue=2,address="Ex1_")

# Volcano plot with FDR cutoff = 0.05, FC cutoff = 70, upper y-axis limit = 100, and no protein name displayed
# FCcutoff=70 is for demonstration purpose
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,type="VolcanoPlot",FCcutoff=70, logBase.pvalue=2, ylimUp=100, ProteinName=FALSE,address="Ex2_")

# show only 'T3-T1' comparisons
# Volcano plot with FDR cutoff = 0.05, FC cutoff = 70, upper y-axis limit = 100, and no protein name displayed
# FCcutoff=70 is for demonstration purpose
# groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,type="VolcanoPlot",FCcutoff=70, logBase.pvalue=2, ylimUp=100, ProteinName=FALSE,which.Comparison="T3-T1",address="Ex3_")

# Heatmap with FDR cutoff = 0.05
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,type="Heatmap", logBase.pvalue=2, address="Ex1_")

# Heatmap with FDR cutoff = 0.05 and FC cutoff = 70
# FCcutoff=70 is for demonstration purpose
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,type="Heatmap",FCcutoff=70, logBase.pvalue=2, address="Ex2_")

# Comparison Plot
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,type="ComparisonPlot",address="Ex1_")

# Comparison Plot
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,type="ComparisonPlot",ylimUp=8,ylimDown=-1,address="Ex2_")


}

