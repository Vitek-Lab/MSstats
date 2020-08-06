#' Import converters from MSstatsConvert package 
#' @importFrom MSstatsConvert MSstatsImport MSstatsClean MSstatsPreprocess
#' @importFrom utils getFromNamespace
.documentFunction = utils::getFromNamespace(".documentFunction", "MSstatsConvert")
.makeAnnotation = utils::getFromNamespace(".makeAnnotation", "MSstatsConvert")
.updateColnames = utils::getFromNamespace(".updateColnames", "MSstatsConvert")

#' Import DIA-Umpire files 
#' 
#' @inheritParams .documentFunction
#' @param raw.frag name of FragSummary_date.xls data, which includes feature-level data.
#' @param raw.pep name of PeptideSummary_date.xls data, which includes selected fragments information.
#' @param raw.pro name of ProteinSummary_date.xls data, which includes selected peptides information.
#' @param annotation name of annotation data which includes Condition, BioReplicate, Run information.
#' @param useSelectedFrag TRUE will use the selected fragment for each peptide. 'Selected_fragments' column is required.
#' @param useSelectedPep TRUE will use the selected peptide for each protein. 'Selected_peptides' column is required.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame with the required format of MSstats.
#'
#' @author Meena Choi, Olga Vitek 
#'
#' @export
#' 
DIAUmpiretoMSstatsFormat = function(
    raw.frag, raw.pep, raw.pro, annotation, useSelectedFrag = TRUE,
    useSelectedPep = TRUE, fewMeasurements = "remove",
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, ...
) {
    input = MSstatsConvert::MSstatsImport(list(Fragments = raw.frag, Peptides = raw.pep, 
                               Proteins = raw.pro), 
                          type = "MSstats", tool = "DIAUmpire", ...)
    input = MSstatsConvert::MSstatsClean(input, use_frag = useSelectedFrag, 
                         use_pept = useSelectedPep)
    annotation = .makeAnnotation(input, annotation)
    input = MSstatsConvert::MSstatsPreprocess(input, annotation,
                              c("PeptideSequence", "FragmentIon"),
                              remove_shared_peptides = TRUE, 
                              remove_single_feature_proteins = removeProtein_with1Feature,
                              feature_cleaning = list(
                                  handle_features_with_few_measurements = fewMeasurements,
                                  summarize_multiple_psms = summaryforMultipleRows
                              ),
                              columns_to_fill = list("PrecursorCharge" = NA,
                                                     "ProductCharge" = NA,
                                                     "IsotopeLabelType" = "L"))
    input
}


#' Import MaxQuant files
#' 
#' @inheritParams .documentFunction
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Raw.file, Condition, BioReplicate, Run, IsotopeLabelType information.
#' @param proteinGroups name of 'proteinGroups.txt' data. It needs to matching protein group ID. If proteinGroups=NULL, use 'Proteins' column in 'evidence.txt'.
#' @param proteinID 'Proteins'(default) or 'Leading.razor.protein' for Protein ID.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @note Warning: MSstats does not support for metabolic labeling or iTRAQ experiments.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
MaxQtoMSstatsFormat = function(
    evidence, annotation, proteinGroups, proteinID = "Proteins", 
    useUniquePeptide = TRUE, summaryforMultipleRows = max, 
    fewMeasurements = "remove", removeMpeptides = FALSE,
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE, ...
) {
    input = MSstatsConvert::MSstatsImport(list(evidence = evidence, protein_groups = proteinGroups), 
                          type = "MSstats", tool = "MaxQuant", ...)
    input = MSstatsConvert::MSstatsClean(input, protein_id_col = proteinID)
    annotation = .makeAnnotation(input, annotation, "Run" = "Rawfile")
    
    m_filter = list(col_name = "PeptideSequence", pattern = "M", 
                    filter = removeMpeptides, drop_column = FALSE)
    oxidation_filter = list(col_name = "Modifications", pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, drop_column = TRUE)
    
    input = MSstatsConvert::MSstatsPreprocess(input, annotation,
                              feature_columns = c("PeptideSequence", "PrecursorCharge"),
                              remove_shared_peptides = TRUE, 
                              remove_single_feature_proteins = removeProtein_with1Peptide,
                              pattern_filtering = list(oxidation = oxidation_filter),
                              feature_cleaning = list(
                                  handle_features_with_few_measurements = fewMeasurements,
                                  summarize_multiple_psms = summaryforMultipleRows
                              ),
                              columns_to_fill = list("PrecursorCharge" = NA,
                                                     "ProductCharge" = NA,
                                                     "IsotopeLabelType" = "L"))
    input
}


#' Generate MSstatsTMT required input format from MaxQuant output
#' 
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param proteinGroups name of 'proteinGroups.txt' data.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mq' for the meaning of each column.
#' @param which.proteinid Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' or 'Gene.names' can be used instead to get the protein ID with single protein. However, those can potentially have the shared peptides.
#' @param rmProt_Only.identified.by.site TRUE will remove proteins with '+' in 'Only.identified.by.site' column from proteinGroups.txt, which was identified only by a modification site. FALSE is the default.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun TRUE will remove PSM with any missing value within each Run. Defaut is FALSE.
#' @param rmPSM_withfewMea_withinRun only for rmPSM_withMissing_withinRun = FALSE. TRUE(default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame of class "MSstatsTMT"
#' 
#' @export
MaxQtoMSstatsTMTFormat = function(
    evidence, proteinGroups, annotation, which.proteinid = 'Proteins',
    rmProt_Only.identified.by.site = FALSE, useUniquePeptide = TRUE,
    rmPSM_withMissing_withinRun = FALSE, rmPSM_withfewMea_withinRun = TRUE,
    rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum, ...
) {
    input = MSstatsConvert::MSstatsImport(list(evidence = evidence, protein_groups = proteinGroups), 
                          "MSstatsTMT", "MaxQuant", ...)
    input = MSstatsConvert::MSstatsClean(input, protein_id_col = which.proteinid, 
                         remove_by_site = rmProt_Only.identified.by.site,
                         channel_columns = "Reporterintensitycorrected")
    annotation = .makeAnnotation(input, annotation)
    few_measurements = ifelse(rmPSM_withfewMea_withinRun, "remove", "keep")
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        list(handle_features_with_few_measurements = few_measurements,
             summarize_multiple_psms = summaryforMultipleRows,
             remove_psms_with_any_missing = rmPSM_withMissing_withinRun)
    )
    colnames(input) = .updateColnames(input, "PrecursorCharge", "Charge")
    input = input[, c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", 
                      "TechRepMixture", "Run", "Channel", "BioReplicate", "Condition", "Intensity")]
    input
}


#' Import OpenMS files
#' 
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from OpenMS, which includes feature(peptide ion)-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
OpenMStoMSstatsFormat = function(
    input, annotation = NULL, useUniquePeptide = TRUE, fewMeasurements = "remove",
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstats", "OpenMS", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = .makeAnnotation(input, annotation)
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge", 
                            "FragmentIon", "ProductCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        list(handle_features_with_few_measurements = fewMeasurements,
             summarize_multiple_psms = summaryforMultipleRows)
    )
    input
}


#' Generate MSstatsTMT required input format for OpenMS output
#' @param input MSstatsTMT report from OpenMS
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun TRUE will remove PSM with any missing value within each Run. Defaut is FALSE.
#' @param rmPSM_withfewMea_withinRun only for rmPSM_withMissing_withinRun = FALSE. TRUE(default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultiplePSMs sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return `data.frame` of class `MSstatsTMT`.
#' 
#' @export
#' 
OpenMStoMSstatsTMTFormat = function(
    input, useUniquePeptide = TRUE, rmPSM_withMissing_withinRun = FALSE,
    rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE,
    summaryforMultiplePSMs = sum, ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstats", "OpenMS", ...)
    input = MSstatsConvert::MSstatsClean(input)
    few_measurements = ifelse(rmPSM_withfewMea_withinRun, "remove", "keep")
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, NULL, 
        feature_columns = c("PeptideSequence", "PrecursorCharge", "Reference", "RetentionTime"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        list(handle_features_with_few_measurements = few_measurements,
             summarize_multiple_psms = summaryforMultiplePSMs,
             remove_psms_with_any_missing = rmPSM_withMissing_withinRun)
    )
    colnames(input) = .updateColnames(input, "PrecursorCharge", "Charge")
    cols = c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture",  "Fraction",
             "TechRepMixture", "Run", "Channel", "Condition", "BioReplicate", "Intensity")
    input = input[, intersect(cols, colnames(input))]
    input
}


#' Import OpenSWATH files
#' 
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from OpenSWATH, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' @param filter_with_mscore TRUE(default) will filter out the features that have greater than mscore_cutoff in m_score column. Those features will be removed.
#' @param mscore_cutoff Cutoff for m_score. Default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#'  
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
OpenSWATHtoMSstatsFormat = function(
    input, annotation, filter_with_mscore = TRUE, mscore_cutoff = 0.01,
    useUniquePeptide = TRUE, fewMeasurements = "remove",
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstats", "OpenSWATH", ...)
    # TODO: check the existence of m_score column earlier
    input = MSstatsConvert::MSstatsClean(input)
    annotation = .makeAnnotation(input, .getDataTable(annotation))
    
    m_score_filter = list(score_column = "m_score", score_threshold = mscore_cutoff, 
                          direction = "smaller", behavior = "remove", 
                          handle_na = "remove", fill_value = NA,
                          filter = TRUE, drop_column = TRUE)
    decoy_filter = list(col_name = "decoy", filter_symbols = 1, 
                        filter = TRUE, drop_column = TRUE)
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        c("PeptideSequence", "PrecursorCharge", "FragmentIon"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(ms_filter = m_score_filter),
        exact_filtering = list(decoy = decoy_filter),
        columns_to_fill = c("ProductCharge" = NA, "IsotopeLabelType" = "L"))
    input
}


#' Import Progenesis files
#' 
#' @inheritParams .documentFunction
#' @param input name of Progenesis output, which is wide-format. 'Accession', 'Sequence', 'Modification', 'Charge' and one column for each run are required.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, Run information. It will be matched with the column name of input for MS runs.
#' @param ... additional parameters to `data.table::fread`.
#'
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek, Ulrich Omasits
#' 
#' @export
#' 
ProgenesistoMSstatsFormat = function(
    input, annotation, useUniquePeptide = TRUE, summaryforMultipleRows = max,
    fewMeasurements = "remove", removeOxidationMpeptides = FALSE, 
    removeProtein_with1Peptide = FALSE, ... 
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstats", "Progenesis", ...)
    annotation = .makeAnnotation(input, annotation)
    input = MSstatsConvert::MSstatsClean(input, unique(annotation$Run), TRUE)
    
    oxidation_filter = list(col_name = "PeptideSequence", pattern = "Oxidation", 
                            filter = TRUE, drop_column = FALSE)
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        c("PeptideSequence", "PrecursorCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        pattern_filtering = list(oxidation = oxidation_filter),
        columns_to_fill = list("FragmentIon" = NA, "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    colnames(input) = .updateColnames(input, "PeptideSequence", "PeptideModifiedSequence")
    input
}


#' Import Proteome Discoverer files
#' 
#' @inheritParams .documentFunction
#' @param input PD report or a path to it.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, 
#' Run information. 'Run' will be matched with 'Spectrum.File'.
#' @param useNumProteinsColumn TRUE removes peptides which have more than 1 in # Proteins column of PD output.
#' @param which.quantification Use 'Precursor.Area'(default) column for quantified intensities. 'Intensity' or 'Area' can be used instead.
#' @param which.proteinid Use 'Protein.Accessions'(default) column for protein name. 'Master.Protein.Accessions' can be used instead.
#' @param which.sequence Use 'Sequence'(default) column for peptide sequence. 'Annotated.Sequence' can be used instead.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
PDtoMSstatsFormat = function(
    input, annotation, useNumProteinsColumn = FALSE, useUniquePeptide = TRUE,
    summaryforMultipleRows = max, fewMeasurements = "remove",
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE,
    which.quantification = 'Precursor.Area', 
    which.proteinid = 'Protein.Group.Accessions', which.sequence = 'Sequence',
    ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstats", "ProteomeDiscoverer", ...)
    input = MSstatsConvert::MSstatsClean(input, quantification_column = which.quantification, 
                         protein_id_column = which.proteinid,
                         sequence_column = which.sequence, 
                         remove_shared = useNumProteinsColumn)
    annotation = .makeAnnotation(input, annotation)
    
    oxidation_filter = list(col_name = "PeptideSequence", pattern = "Oxidation", 
                            filter = TRUE, drop_column = FALSE)
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        pattern_filtering = list(oxidation = oxidation_filter),
        columns_to_fill = list("FragmentIon" = NA, "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    colnames(input) = .updateColnames(colnames(input), "PeptideSequence",
                                      "PeptideModifiedSequence")
    input
}


#' Convert Proteome Discoverer output to MSstatsTMT format.
#' 
#' @param input PD report or a path to it.
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column.
#' @param which.proteinid Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' or 'Gene.names' can be used instead to get the protein ID with single protein. However, those can potentially have the shared peptides.
#' @param useNumProteinsColumn logical, if TRUE, shared peptides will be removed.
#' @param useUniquePeptide lgl, if TRUE (default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun lgl, if TRUE, will remove PSM with any missing value within each Run. Default is FALSE.
#' @param rmPSM_withfewMea_withinRun lgl, only for rmPSM_withMissing_withinRun = FALSE. TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @export
#' 
PDtoMSstatsTMTFormat <- function(
    input, annotation, which.proteinid = 'Protein.Accessions', 
    useNumProteinsColumn = TRUE, useUniquePeptide = TRUE, 
    rmPSM_withMissing_withinRun = FALSE, rmPSM_withfewMea_withinRun = TRUE, 
    rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum, ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstatsTMT", "ProteomeDiscoverer", ...)
    input = MSstatsConvert::MSstatsClean(input, protein_id_column = which.proteinid,
                         remove_shared = useNumProteinsColumn)
    annotation = .makeAnnotation(input, annotation)
    
    few_measurements = ifelse(rmPSM_withfewMea_withinRun, "remove", "keep")
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        feature_cleaning = list(handle_features_with_few_measurements = few_measurements,
                                summarize_multiple_psms = summaryforMultipleRows,
                                remove_psms_with_any_missing = rmPSM_withMissing_withinRun)
    )
    colnames(input) = .updateColnames(input, "PrecursorCharge", "Charge")
    input = input[, c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", 
                      "TechRepMixture", "Run", "Channel", "Condition", "BioReplicate", "Intensity")] # unique?
    input
}


#' Import Skyline files
#'
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from Skyline, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Skyline, use annotation=NULL (default). It will use the annotation information from input.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which are labeld 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in DetectionQValue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
SkylinetoMSstatsFormat = function(
    input, annotation = NULL, removeiRT = TRUE, filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01, useUniquePeptide = TRUE, fewMeasurements = "remove",
    removeOxidationMpeptides = FALSE, removeProtein_with1Feature = FALSE, ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstats", "Skyline", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = .makeAnnotation(input, annotation, Run = "Rawfile")
    
    decoy_filter = list(col_name = "ProteinName",
                        filter_symbols = c("DECOY", "Decoys"),
                        filter = TRUE, drop_column = FALSE)
    irt_filter = list(col_name = "StandardType", filter_symbols = "iRT",
                      filter = removeiRT, drop_column = FALSE)
    oxidation_filter = list(col_name = "PeptideSequence",
                            pattern = "\\+16", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    truncated_filter = list(score_column = "Truncated", score_threshold = 0,
                            direction = "smaller", behavior = "fill",
                            fill_value = NA_real_, handle_na = "keep",
                            filter = TRUE, drop_column = TRUE)
    qval_filter = list(score_column = "DetectionQValue", 
                       score_threshold = qvalue_cutoff, direction = "smaller",
                       behavior = "fill", fill_value = 0, handle_na = "keep",
                       filter = filter_with_Qvalue, drop_column = TRUE)
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        score_filtering = list(truncated = truncated_filter, qval = qval_filter),
        pattern_filtering = list(oxidation = oxidation_filter),
        exact_filtering = list(decoy = decoy_filter, irt = irt_filter),
        aggregate_isotopic = TRUE,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = sum)
    )
    input
}


#' Import data from SpectroMine
#'
#' @param input data name of SpectroMine PSM output. Read PSM sheet.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mine' for the meaning of each column.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with NA and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun TRUE will remove PSM with any missing value within each Run. Defaut is FALSE.
#' @param rmPSM_withfewMea_withinRun only for rmPSM_withMissing_withinRun = FALSE. TRUE(default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @export 
#' 
SpectroMinetoMSstatsTMTFormat <- function(
    input, annotation, filter_with_Qvalue = TRUE, qvalue_cutoff = 0.01,
    useUniquePeptide = TRUE, rmPSM_withMissing_withinRun = FALSE,
    rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE,
    summaryforMultipleRows = sum, ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstatsTMT", "SpectroMine", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = .makeAnnotation(input, annotation)
    
    few_measurements = ifelse(rmPSM_withfewMea_withinRun, "remove", "keep")
    pq_filter = list(score_column = "PGQValue", score_threshold = 0.01, 
                     direction = "smaller", behavior = "fill", 
                     handle_na = "keep", fill_value = NA,
                     filter = TRUE, drop_column = TRUE)
    qval_filter = list(score_column = "Qvalue", score_threshold = qvalue_cutoff, 
                       direction = "smaller", behavior = "fill", 
                       handle_na = "keep", fill_value = NA,
                       filter = filter_with_Qvalue, drop_column = TRUE)
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        score_filtering = list(pgq = pq_filter, psm_q = qval_filter),
        feature_cleaning = list(handle_features_with_few_measurements = few_measurements,
                                summarize_multiple_psms = summaryforMultipleRows,
                                remove_psms_with_any_missing = rmPSM_withMissing_withinRun)
    )
    colnames(input) = .updateColnames(input, "PrecursorCharge", "Charge")
    input = input[, c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", 
                      "TechRepMixture", "Run", "Channel", "BioReplicate", "Condition", "Intensity")]
    input
}


#' Import Spectronaut files
#' 
#' @param input name of Spectronaut output, which is long-format. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity, F.ExcludedFromQuantification are required. Rows with F.ExcludedFromQuantification=True will be removed.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Spectronaut, use annotation=NULL (default). It will use the annotation information from input.
#' @param intensity 'PeakArea'(default) uses not normalized peak area. 'NormalizedPeakArea' uses peak area normalized by Spectronaut.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .documentFunction
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
SpectronauttoMSstatsFormat = function(
    input, annotation = NULL, intensity = 'PeakArea', filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01, useUniquePeptide = TRUE, fewMeasurements="remove",
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, ...
) {
    input = MSstatsConvert::MSstatsImport(list(input = input), "MSstats", "Spectronaut", ...)
    input = MSstatsConvert::MSstatsClean(input, intensity = intensity)
    annotation = .makeAnnotation(input, annotation, "Run" = "RFileName", 
                                 "Condition" = "RCondition", "BioReplicate" = "RReplicate")
    
    pq_filter = list(score_column = "PGQvalue", score_threshold = 0.01, 
                     direction = "smaller", behavior = "fill", 
                     handle_na = "keep", fill_value = NA,
                     filter = TRUE, drop_column = TRUE)
    qval_filter = list(score_column = "Qvalue", score_threshold = qvalue_cutoff, 
                       direction = "smaller", behavior = "fill", 
                       handle_na = "keep", fill_value = 0, 
                       filter = filter_with_Qvalue, drop_column = TRUE)
    
    input = MSstatsConvert::MSstatsPreprocess(
        input, annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(pgq = pq_filter, psm_q = qval_filter),
        columns_to_fill = list("IsotopeLabelType" = "L"))
    input
}
