\name{dataProcess}
\alias{dataProcess}
\title{Data pre-processing and quality control of MS runs of raw data}
\description{Data pre-processing and quality control of MS runs of the original raw data into quantitative data for model fitting and group comparison. Log transformation is automatically applied and additional variables are created in columns for model fitting and group comparison process. Three options of data pre-processing and quality control of MS runs in dataProcess are (1) Transformation: logarithm transformation with base 2 or 10; (2) Normalization: to remove systematic bias between MS runs.}

\usage{
dataProcess(raw, logTrans=2, normalization="constant",nameStandards=NULL,betweenRunInterferenceScore=FALSE,address="")	
}

\arguments{
  \item{raw}{name of the raw (input) data set.}
  \item{logTrans}{logarithm transformation with base 2(default) or 10.}
  \item{normalization}{normalization to remove systematic bias between MS runs. There are three different normalizations supported. 'constant'(default) represents constant normalization based on reference signals is performed. 'quantile' represents quantile normalization based on reference signals is performed. 'globalStandards' represents normalization with global standards proteins. FALSE represents no normalization is performed.}
  \item{nameStandards}{vector of global standard protein names. only for normalization with global standard proteins. }
 \item{betweenRunInterferenceScore}{interference is detected by a between-run-interference score. TRUE means the scores are generated automatically and stored in a .csv file. FALSE(default) means no scores are generated.}
\item{address}{the name of folder that will store the results. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory. An output csv file is automatically created with the default name of "BetweenRunInterferenceFile.csv". The command address can help to specify where to store the file as well as how to modify the beginning of the file name.}
}

\details{
\itemize{
\item{raw : See \code{\link{RawData}} for the required data structure of raw (input) data.}
\item{logTrans : if logTrans=2, the measurement of Variable ABUNDANCE is log-transformed with base 2. Same apply to logTrans=10.}
\item{normalization : if normalization=TRUE and logTrans=2, the measurement of Variable ABUNDANCE is log-transformed with base 2 and normalized. Same as for logTrans=10.}

}
}

\section{Warning}{When a transition is missing completely in a condition or a MS run, a warning message is sent to the console notifying the user of the missing transitions.  
	
	The types of experiment that MSstats can analyze are LC-MS, SRM, DIA(SWATH) with label-free or labeled synthetic peptides. MSstats does not support for metabolic labeling or iTRAQ experiments.
}

\author{Ching-Yun Chang, Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{choi67@purdue.edu})
}
\references{
Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements"  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\examples{
# Consider a raw data (i.e. RawData) for a label-based SRM experiment from a yeast study with ten time points (T1-T10) of interests and three biological replicates. 
# It is a time course experiment. The goal is to detect protein abundance changes across time points.

head(RawData)

# Log2 transformation and normalization are applied (default)
QuantData<-dataProcess(RawData)
head(QuantData)

# Log10 transformation and normalization are applied
QuantData1<-dataProcess(RawData, logTrans=10)
head(QuantData1)

# Log2 transformation and no normalization are applied
QuantData2<-dataProcess(RawData,normalization=FALSE)
head(QuantData2)


}

