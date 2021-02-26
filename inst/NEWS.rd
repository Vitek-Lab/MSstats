\name{MSstatsnews}
\title{News for package, \pkg{MSstats}}
\encoding{UTF-8}

\section{Version 3.21.4 (2021-02-25)}{
    \itemize{
        \item Fix the bug for feature selection
    }
}

\section{Version 3.21.2 (2020-08-19)}{
    \itemize{
        \item Fix the bugs for converter.(dcast updated)
        \item Remove CV in sample size calculation plot.
    }
}

\section{Version 3.18.5 (2020-02-26)}{
    \itemize{
        \item Fix the error for feature selection option in dataProcess.
    }
}

\section{Version 3.18.4 (2019-12-17)}{
    \itemize{
        \item dataProcessPlots : bug to draw profile plot. The option for scale_colour_manual is updated.
    }
}

\section{Version 3.18.3 (2019-12-10)}{
    \itemize{
        \item groupComparisonPlots : add 'text.angle' option.
    }
}

\section{Version 3.18.2 (2019-11-22)}{
    \itemize{
        \item Fixed nonlinear_quantlim error.
    }
}

\section{Version 3.18.1 (2019-10-31)}{
    \itemize{
        \item Fixed message: transitions are completely missing in one condition.-> transitions are completely missing in at least one of the conditions.
        \item Fixed the bug : NumMeasuredFeature and NumImputedFeature in RunlevelData from dataProcess.
        \item change NEWS file to rd file.
    }
}

\section{Version 3.18.0 and 3.17.2 Bioconductor 3.10 Release (2019-10-30)}{
    \itemize{
        \item Fixed warning message: Warning message:replacing previous import ‘MASS::select’ by ‘dplyr::select’ when loading ‘MSstats’, is fixed.
        \item RETIRED FEATURES : designSampleSizeClassification and designSampleSizeClassificationPlots are retired. The new package, MSstatsSampleSize, is released in Bioconductor 3.10.
    }
}

\section{Version 3.16.2 (2019-08-15)}{
    \itemize{
        \item dataProcess, featureSubset='top3' or 'topN' with n_top_feature : fix the bug for featureSubset = 'top3' and 'topN' option and improve the memory consumption.
        \item groupComparison : remove the warning message for singularity issue.
    }
}

\section{Version 3.16.1 (2019-05-07)}{
    \itemize{
        \item MaxQtoMSstatsFormat : fix the bug when no + in Only.identified.by.site column in proteinGroups input.
    }
}

\section{Version 3.16.0 (2019-05-03)}{
    \itemize{
        \item UPDATED FEATURES : featureSubset="highquality", remove_uninformative_feature_outlier=TRUE available.
    }
}

\section{Version 3.13.5 (2018-07-30)}{
    \itemize{
        \item MaxQtoMSstatsFormat : 'Leading.razor.protein' in evidence can be used for protein IDs.
        \item MaxQtoMSstatsFormat : 'Run' column is not required in annotation.
    }
}
     
\section{Version 3.13.4 (2018-07-02)}{
    \itemize{
        \item UPDATED FEATURES : designSampleSizeClassification and designSampleSizeClassificationPlots. 
        \item UPDATED FEATURES : Converter functions check the required column of annootation file.
    }
}

\section{Version 3.13.3 (2018-06-01)}{
    \itemize{
        \item PDtoMSstatsFormat : adjust for different column ids from PD.
    }
}	 

\section{Version 3.13.2 (2018-05-29)}{
    \itemize{
        \item dataProcess : fix some reports.
    }
}

\section{Version 3.13.1 (2018-05-21)}{
    \itemize{
        \item RETIRED FEATURES : transformMSnSetToMSstats, transformMSstatsToMSnSet are retired.
        \item dataProcess : report summary of proteins with single feature.
        \item dataProcess : cluster=NULL is default. There is the issue for makecluster.
        \item dataProcess : fix the bug for fraction case.
    }
}

\section{Version 3.11.6 (2018-04-23)}{
    \itemize{
        \item SkylinetoMSstatsFormat : fix the inconsistency of column name from Skyline output.
    }
}	 

\section{Version 3.11.5 (2018-02-22)}{
    \itemize{
        \item add the package, stringr, for DIAUmpiretoMSstatsFormat function.
    }
}

\section{Version 3.11.4 (2018-02-19)}{
    \itemize{
        \item add set.seed for sample size calculation of classification.
    }
}

\section{Version 3.11.3 (2018-02-15)}{
    \itemize{
        \item nonlinear_quantlim : fix the bug for the resampling of the blank sample, increase the default number of bootstrap samples.
        \item designSampleSize : fix the bug.
        \item NEW FEATURES : new function : designSampleSizeClassification, designSampleSizeClassificationPlots - Calculate the optimal size of training data for classification problem by simulation.
        \item NEW FEATIRES : new converter functions : DIAUmpiretoMSstatsFormat, OpenMStoMSstatsFormat.
    }
}

\section{Version 3.10.5 (2018-01-10)}{
    \itemize{
        \item SpectronauttoMSstatsFormat : TRUE or FALSE are allowed for the values of the column, F.ExcludedFromQuantification. Check the value for this column.
    }
}

\section{Version 3.10.4 (2017-12-22)}{
    \itemize{
        \item MaxQtoMSstatsFormat : 'fewmeasurements' bug fixed.
    }
}

\section{Version 3.10.2 (2017-11-27)}{
    \itemize{
        \item make error messages for QQ plot and residual plot, if the protein couldn't be fitted by linear mixed effect model.
        \item ProgenesistoMSstatsFormat : make more generalization for different format.
    }
}	
	
\section{Version 3.9.7 (2017-10-23)}{
    \itemize{
        \item NEW FEATURES : cluster (default=1) is no longer available for groupComparison function, due to memory issue.
    }
}

\section{Version 3.8.6 (2017-09-26)}{
    \itemize{
        \item NEW FEATURES : can cluster (default=1) for dataProcess and groupComparison function (Thanks John!!).
    }
}

\section{Version 3.8.5 (2017-09-25)}{
    \itemize{
        \item PDtoMSstatsFormat : three options are added for outputs from different versions of PD. (Thanks to Felipe!) - which.quantification, which.proteinid, which.sequence
    }
}

\section{Version 3.8.4 (2017-08-28)}{
    \itemize{
        \item SkylinetoMSstatsFormat : DDA case lost ‘StandardType’ column after summing peaks. Fixed. (Thanks, Nick)
    }
}

\section{Version 3.8.3 (2017-07-13)}{
    \itemize{
        \item NEW FEATURES - SpectronauttoMSstatsFormat : if PG.Qvalue is available, filter out the data with greater than 0.01.
        \item NEW FEATURES - dataProcessPlots, groupComparisonPlots, modelBasedPlots : with address=FALSE option, one plot a time can be drawn in the panel and won't be saved as in pdf.
        \item designSampleSize : fix the calculation of variance (Thanks, Tsung-Heng)
        \item skylinetoMSstatsFormat : when Condition and BioReplicate columns are NA, there was issue for merge with annotation.
        \item SkylinetoMSstatsFormat : fix the bug to recognize the protein with one peptide only for the option: 'removeProtein_with1Peptide = TRUE'
        \item dataProcess : when cutoff.lower is negative, with maxQuantileforCensored option + censoredInt='0', zero log2 endogenous intensity should be censored.
        \item ProgenesistoMSstatsFormat : handle inputs with some limited columns. such as no Spectral.counts columns.
    }
}

\section{Version 3.8.2 (2017-04-21)}{
    \itemize{
        \item NEW FEATURES - required ‘Fraction’ information in annotation for pre-processing.
        \item NEW FEATURES - dataProcess function is updated for merge fractions.
        \item warning message during dataProcessPlots for profile plot is not shown anymore.
    }
}

\section{Version 3.7.4 (2017-04-21)}{
    \itemize{
        \item quantile normalization in dataProcess
        \item Show the progress for comparison plots in groupComparisonPlots.
        \itemize{
            \item Input check : whether any annotation information is missing or not.
            \item SpectronauttoMSstatsFormat function
        }
    }
}

\section{Version 3.7.3 (2017-01-16)}{
    \itemize{
        \item NEW FEATURES: ProgenesistoMSstatsFormat function : required input for this function is changed.
        \item NEW FEATURES: SpectronauttoMSstatsFormat : required input for this function is changed.
        \item NEW FEATURES: PDtoMSstatsFormat : new function for converting Proteome Discoverer output to MSstats required format.
        \item groupComparisionPlots function : when assigning FCcutoff, there is error sometimes. (It was due to 'logFC' vs 'log2FC' vs 'log10FC' in output of groupComparison. groupComparison function is updated.)
        \item dataProcess function : error with normalization='quantile' is fixed. (It was due to absense of 'Import(preprocessCore)' in the namespace and assigning 'originalRun' column.)
    }
}

\section{Version 3.5.5 (2016-09-30)}{
    \itemize{
        \item SIGNIFICANT CHANGES FOR METHOD
        \itemize{
            \item applied the algorithm for deciding the threshold of censoring.
            \item Method for calculation of the LOB/LOD is changed. LOQ is not calculated anymore. Please check help files for details.
            \item summaryMethod='logOfSum' option in dataProcess is retired.
            \item modelBasedQCPlots work with output from groupComparison in order to check the normality assumption of linear mixed effect model for whole plot level inference.
        }
        \item NEW FEATURES: Add ‘originalRUN’ column in xx$ProcessedData after dataProcess function.
        \item NEW FEATURES: Profile plot from dataProcessPlot distinguish censored missing data or not with different symbol.
        \item dataProcess with fractionation sample when filling incomplete rows. Especially, not balanced fractionation for heavy and light, (heavy in one fractionation, no heavy in other fractionation)
        \item groupComparison function : fix the issue with different columns from different summary Methods.
        \item MaxQtoMSstats function : option removeMpeptides=FALSE are now available. (Thanks, Danielle)
        \item In case of multiple injections per sample (from fractionation or multiple injections with different range of m/x), normalization is performed separately and multiple injections are merged in each sample.
    }
}

\section{Version 3.5.1}{
    \itemize{
        \item Fix bug : summaryMethod=‘logOfSum’, redesign for result table.
    }
}

\section{Version 3.3.11}{
    \itemize{
        \item New functionalities : calculation of the LOD and LOQ, 1) linear_quantlim, 2) nonlinear_quantlim, 3) plot_quantlim, and two example datasets, SpikeInDataLinear, SpikeInDataNonLinear are available.
        \item Update for featureSelection =‘highQuality’ in dataProcess.
        \item allow colon(“:”) in the peptide sequence
        \item fix the bug for ‘fill in incomplete rows’ in dataProcess. If there are only one feature has incomplete rows, the issue for getting run and feature ID in dataProcess and not show the list. Now, it works.
        \item change the default for blimp in dataProcessPlots for profile plot and QC plot. The upper limit of y-axis with ylimUp=FALSE is calculated by maximum log2(intensity) across all proteins after normalization + 3 and then rounding off to the nearest integer.
    }
}

\section{Version 3.3.10}{
    \itemize{
        \item dataProcess : When the number of proteins for $ProcessedData and $RunlevelData are different,  the bug happened for calculating the number of missing and imputation.
        \item groupComparison : When one of condition is completely missing or other special case, .fit.model.single can handle and output of .fit.model.single is not try-error. Then output for fitted and residual should be updated.
    }
}

\section{Version 3.3.9}{
    \itemize{
        \item Condition plot from dataProcessPlots : Now condition plots are drawn with run-level summarized intensities per condition.
        \item ComparisonResult from groupComparison 
        \itemize{
            \item flag about missingness and imputation : Calculation for MissingPercentage and ImputationPercentage columns is changed.
            \itemize{
                \item MissingPercentage : number of measured intensities/ total number of intensities (which is the number of features * the number of runs in a protein) in the conditions used for comparison (from ‘Label’ column) by protein. Therefore different comparisons(Label in the output) from the same protein can have the different percentage of missingness.
                \item ImputationPercentage : number of imputed intensities/total number of intensities in the conditions used for comparison (from ‘Label’ column) by protein. Therefore different comparisons(Label in the output) from the same protein can have the different percentage of imputation.
            }
            \item new column, ‘issue’, shows special cases, such as completely missing in a condition or all conditions for comparisons.
        }
        \item VolcanoPlot - flag the proteins which have the condition with completely missing. On the left of protein name, ‘*’ willz be appeared in Volcano plot
    }
}

\section{Version 3.3.8}{
    \itemize{
        \item normalization : overall median -> median of medians. For all workflow for MSstats, the result should not be changed. But, log(sum) will have slightly different result.
        \item flag about missingness and imputation
        \itemize{
            \item RunlevelData from dataProcess include two or three more columns
            \itemize{
                \item NumMeasuredFeature : number of measured features in a run
                \item Missing percentage : number of measured features/total number of features by run
                \item NumImputedFeature : number of imputed intensities in a run. This column is shown only if users allow to impute the missing value.
            }
            \item ComparisonResult from groupComparison : one or two columns will be added.
            \itemize{
                \item MissingPercentage : number of measured intensities/ total number of intensities (which is the number of features * the number of runs in a protein) by protein
                \item ImputationPercentage : number of imputed intensities/total number of intensities by protein
            }
        }
    }
}


\section{Version 3.3.4}{
    \itemize{
        \item fix the bug for featureSubset=‘highQuality’ with label-based experiment.
    }
}

\section{Version 3.3.3}{
    \itemize{
        \item add new option, remove_proteins_with_interference=FALSE (default), in dataProcess. whether it allows to remove proteins if deem interfered or not.
    }
}

\section{Version 3.3.2}{
    \itemize{
        \item ProteinName=TRUE in groupComparisonPlots shows only the name of significant proteins, adjusting location. ggrepel package is used.
        \item change featureSubset option in ‘dataProcess’ for selecting high quality features. featureSubset=‘highQuality’
        \item Fix the bug for ‘fillIncompleteRows=TRUE’ for label-based experiments.
        \item change ‘quantification’ function. use run summarization from dataProcess. If there are technical replicates, use median run summarization for each subject.
    }
}

\section{Version 3.3.1}{
    \itemize{
        \item fix the bug for volcano plot in groupComparisonPlots, with logbase=2.
        \item update all plots for ggplot2
        \item Change the default for ‘cutoffCensored’. Now the default is “minFeature”.
        \item for imputing the censored peak intensities, remove the features which has only 1 measurement for survreg function.
    }
}

\section{Version 3.2.3}{
    \itemize{
        \item bug for normalization=‘globanstandards’ in dataProcess. Even though wrong normalization options, data processing is kept going. Make the process stop if wrong input for normalization
        \item fix the bug with nameStandards=‘protein name’. (only peptide names worked with 3.2.2.)
    }
}

\section{Version 3.0.12}{
    \itemize{
        \item remove vignetter folder to remove install and build error in Bioconductor
    }
}

\section{Version 3.0.9}{
    \itemize{
        \item dataProcess
        \itemize{
            \item add options for ‘cutoffCensored=“minFeatureNRun”’.
            \item summaryMethods=“TMP” : output will have ‘more50missing’column.
            \item remove50missing=FALSE option : remove runs which has more than 50\% of missing measurement. It will be affected for TMP, with censored option.
            \item MBimpute : impute censored by survival model (AFT) with cutoff censored value
            \item featureSubset option : “all”,”top3”, “highQuality”
            \item change the default.
        }
        \item groupComparisonPlots : heatmap, for logBase=10, fix the bug for setting breaks.
    }
}

\section{Version 3.0.8}{
    \itemize{
        \item dataProcess : when censoredInt=“0”, intensity=0 works even though skylineReport=FALSE.
        \item dataProcess, with censored=“0” or “NA” : fix the bug for certain run has completely missing.
        \item cutoffCensored=“minRun” or “minFeature” : cutoff for each Run or each feature is little less (99\%) than minimum abundance.
        \item summaryMethod=“TMP”, censored works. censoredInt=NA or 0, and cutoffCensored=0, minFeature, minRun
    }
}

\section{Version 3.0.3}{
    \itemize{
        \item dataProcess : new option, skylineReport. for skyline MSstats report, there is ‘Truncated’ column. If Truncated==True, remove them. and keep zero value for summaryMethod=“skyline”.
        \item groupComparison : for skyline option, t.test, val.equal=TRUE, which is no adjustment for degree of freedom, just pooled variance.
    }
}

\section{Version 2.1.4}{
    \itemize{
        \item Change the numbering of file name for plots. If the file names are already existing in the folder, automatically next number will be add at the end of file name. Therefore file name will not be overwritten, and we can keep all pdf files for plots.
    }
}


\section{Version 2.1.3}{
    \itemize{
        \item fix the groupComparison for label-free experiments.
        \item automatically generate progress report as .txt files
        \item add progress message for groupComparison and dataProcessPlots function.
    }
}

\section{Version 2.1.1}{
    \itemize{
        \item fix the bug in Condition plot : 1. for label-based : match reference and endogenous, 2. for label-free : when there is one observation in each group, SD=NA. make it zero.
        \item fix the bug in heatmap and comparison plots : remove NA result for plotting
        \item fix the bug for label-free groupComparison :  how to get subject_nested parameter in make.contrast.free for unequal number per group
        \item fix the bug in group quantification : make.contrast.group.quantification fixed for subject_nested parameter
    }
}

\section{Version 1.99.1}{
    \itemize{
        \item fixed several NOTES, added .Rbuildignore, compacted vignettes
        \item removed warn -1
        \item added validity check when returning MSnSet
        \item used inherits/is. for class testing
    }
}

\section{Version 1.99.0}{
    \itemize{
        \item improve efficiency for computing groupComparison and quantification <2012-12-21>
        \item add .rnw <2012-12-03>
        \item update groupComparision for label-free time-course experiment with single Feature and with or without technical replicates <2013-04-08>
        \item add option for saving QQ plot and Residual plot in order to checkin the normality assumption in groupComparison function. <2013-04-08>
        \item use ggplot2 package for all plots. <2013-07-11>
        \item fix bug for volcano plot : different color labeling <2013-07-12>
        \item add power plot in sample size calculation plot <2013-07-12>
        \item add 'interference=TRUE or FALSE' in sample size calculation <2013-07-15>
        \item add 'legend.size=7'for size of feature name legends in dataProcessPlots <2013-07-23> 
        \item add 'text.angle=0'for angle of condition labeling in dataProcessPlots <2013-07-23>
        \item fix bug for quantification : when there are missing values in endogenous intensities, but values in reference intensities. <2013-07-24>
        \item fix bug for groupComparison : when there are missing values in endogenous intensities, but values in reference intensities, .make.constast.based or .free sub function were changed. <2013-07-25>
        \item two function for transformation between required input for MSstats and MSnSet class <2013-09-04>
        \item flexibility for visualization : save as pdf files or show in window with selected proteins or all proteins. <2013-09-04>
        \item handle unequal variance for feature in groupComparison function with featureVar=TRUE <2013-09-04>
        \item Add 'missing.action' for impute missing values in group comparison stage. <2013-09-20>
    }
}