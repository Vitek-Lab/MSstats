\name{dataProcess}
\alias{dataProcess}
\title{Data pre-processing and quality control of MS runs of raw data}
\description{Data pre-processing and quality control of MS runs of the original raw data into quantitative data for model fitting and group comparison. Log transformation is automatically applied and additional variables are created in columns for model fitting and group comparison process. Two options of data pre-processing and quality control of MS runs in dataProcess are (1) Transformation: logarithm transformation with base 2 or 10; (2) Normalization: to remove systematic bias between MS runs.}

\usage{
dataProcess(raw,
            logTrans=2,
            normalization="equalizeMedians",
            nameStandards=NULL,
            address="",
            fillIncompleteRows=TRUE,
            featureSubset="all",
            remove_uninformative_feature_outlier=FALSE,
            n_top_feature=3,
            summaryMethod="TMP",
            equalFeatureVar=TRUE,
            censoredInt="NA",
            cutoffCensored="minFeature",
            MBimpute=TRUE,
            remove50missing=FALSE,
            maxQuantileforCensored=0.999,
            clusters=NULL)
}

\arguments{
  \item{raw}{name of the raw (input) data set.}
  \item{logTrans}{logarithm transformation with base 2(default) or 10.}
  \item{normalization}{normalization to remove systematic bias between MS runs. There are three different normalizations supported. 'equalizeMedians'(default) represents constant normalization (equalizing the medians) based on reference signals is performed. 'quantile' represents quantile normalization based on reference signals is performed. 'globalStandards' represents normalization with global standards proteins. FALSE represents no normalization is performed.}
  \item{nameStandards}{vector of global standard peptide names. only for normalization with global standard peptides. }
 \item{fillIncompleteRows}{If the input dataset has incomplete rows, TRUE(default) adds the rows with intensity value=NA for missing peaks. FALSE reports error message with list of features which have incomplete rows. }
  \item{featureSubset}{"all"(default) uses all features that the data set has. "top3" uses top 3 features which have highest average of log2(intensity) across runs. "topN" uses top N features which has highest average of log2(intensity) across runs. It needs the input for n_top_feature option. "highQuality" flags uninformative feature and outliers.}
  \item{remove_uninformative_feature_outlier}{It only works after users used featureSubset="highQuality" in dataProcess. TRUE allows to remove 1) the features are flagged in the column, feature_quality="Uninformative" which are features with bad quality, 2) outliers that are flagged in the column, is_outlier=TRUE, for run-level summarization. FALSE (default) uses all features and intensities for run-level summarization.}
  \item{n_top_feature}{The number of top features for featureSubset='topN'. Default is 3, which means to use top 3 features.}
\item{summaryMethod}{"TMP"(default) means Tukey's median polish, which is robust estimation method. "linear" uses linear mixed model.}
\item{equalFeatureVar}{only for summaryMethod="linear". default is TRUE. Logical variable for whether the model should account for heterogeneous variation among intensities from different features. Default is TRUE, which assume equal variance among intensities from features. FALSE means that we cannot assume equal variance among intensities from features, then we will account for heterogeneous variation from different features.}
\item{censoredInt}{Missing values are censored or at random. 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. '0' uses zero intensities as censored intensity. In this case, NA intensities are missing at random. The output from Skyline should use '0'. Null assumes that all NA intensites are randomly missing.}
\item{cutoffCensored}{Cutoff value for censoring. only with censoredInt='NA' or '0'. Default is 'minFeature', which uses minimum value for each feature.'minFeatureNRun' uses the smallest between minimum value of corresponding feature and minimum value of corresponding run. 'minRun' uses minumum value for each run. }
\item{MBimpute}{only for summaryMethod="TMP" and censoredInt='NA' or '0'. TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by Accelated failure model. FALSE uses the values assigned by cutoffCensored.}
\item{remove50missing}{only for summaryMethod="TMP". TRUE removes the runs which have more than 50\% missing values. FALSE is default.}
\item{address}{the name of folder that will store the results. Default folder is the current working directory. The other assigned folder has to be existed under the current working directory. An output csv file is automatically created with the default name of "BetweenRunInterferenceFile.csv". The command address can help to specify where to store the file as well as how to modify the beginning of the file name.}
\item{maxQuantileforCensored}{Maximum quantile for deciding censored missing values. default is 0.999}
\item{clusters}{a user specified number of clusters. default is NULL, which does not use cluster.}

}

\details{
\itemize{
\item{raw : See \code{\link{SRMRawData}} for the required data structure of raw (input) data.}
\item{logTrans : if logTrans=2, the measurement of Variable ABUNDANCE is log-transformed with base 2. Same apply to logTrans=10.}
\item{normalization : if normalization=TRUE and logTrans=2, the measurement of Variable ABUNDANCE is log-transformed with base 2 and normalized. Same as for logTrans=10.}
\item{featureSubset : After the data was normalized, we deeply looked at each single feature (which is a precursor in DDA, a fragment in DIA, and a transition in SRM) and quantify its un-explainable variation. Ultimately, we remove the features with interference.}
\item{equalFeatureVar : If the unequal variation of error for different peptide features is detected, then a possible solution is to account for the unequal error variation by means of a procedure called iteratively re-weighted least squares. equalFeatureVar=FALSE performs an iterative fitting procedure, in which features are weighted inversely proportionaly to the variation in their intensities, so that feature with large variation are given less importance in the estimation of parameters in the model.}

\item{LogIntensities in RunlevelData from dataProcess : Estimated protein abundance for aa`RUN` and `Protein`.}
\item{NumMeasuredFeature in RunlevelData from dataProcess : Number of measured features in a `RUN` and `Protein`. It counts feature intensities greater than 1, after selecting all/topN/highquality features and removing features with only one measurement across MS runs.}
\item{Missing percentage in RunlevelData from dataProcess : Number of measured features in a `RUN` and `Protein` /total number of features in a `Proteins`.}
\item{NumImputedFeature in RunlevelData from dataProcess : Number of imputed intensities (censored intensities) in a `RUN` and `Protein`. This column is shown only if users allow to impute the missing value.}

}
}

\value{A list of data.frame \emph{ProcessedData} is the data.frame of reformatted input of dataProcess including extra columns, such as log2-transformed and normalized intensities (abundance column); \emph{RunlevelData} is the the data.frame for run-level summarized data.}

\section{Warning}{When a transition is missing completely in a condition or a MS run, a warning message is sent to the console notifying the user of the missing transitions.

	The types of experiment that MSstats can analyze are LC-MS, SRM, DIA(SWATH) with label-free or labeled synthetic peptides. MSstats does not support for metabolic labeling or iTRAQ experiments.
}

\author{Meena Choi, Tsung-Heng Tsai, Olga Vitek.

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}
\references{
Meena Choi, Ching-Yun Chang, Timothy Clough, Daniel Broudy, Trevor Killeen, Brendan MacLean and Olga Vitek. "MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments" \emph{Bioinformatics}, 30(17):2524-2526, 2014.

Ching-Yun Chang, Paola Picotti, Ruth Huttenhain, Viola Heinzelmann-Schwarz, Marko Jovanovic, Ruedi Aebersold, Olga Vitek. "Protein significance analysis in selected reaction monitoring (SRM) measurements"  \emph{Molecular & Cellular Proteomics}, 11:M111.014662, 2012.

Timothy Clough, Safia Thaminy, Susanne Ragg, Ruedi Aebersold, Olga Vitek. "Statistical protein quantification and significance analysis in label-free LC-M experiments with complex designs" \emph{BMC Bioinformatics}, 13:S16, 2012.
}

\examples{
# Consider a raw data (i.e. SRMRawData) for a label-based SRM experiment from a yeast study
# with ten time points (T1-T10) of interests and three biological replicates.
# It is a time course experiment. The goal is to detect protein abundance changes
# across time points.

head(SRMRawData)

# Log2 transformation and normalization are applied (default)
QuantData<-dataProcess(SRMRawData)
head(QuantData$ProcessedData)

# Log10 transformation and normalization are applied
QuantData1<-dataProcess(SRMRawData, logTrans=10)
head(QuantData1$ProcessedData)

# Log2 transformation and no normalization are applied
QuantData2<-dataProcess(SRMRawData,normalization=FALSE)
head(QuantData2$ProcessedData)


}

