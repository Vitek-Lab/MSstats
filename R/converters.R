#' A dummy function to store shared documentation items.
#' 
#' @import data.table
#' @importFrom MSstatsConvert MSstatsImport MSstatsClean MSstatsPreprocess 
#' MSstatsBalancedDesign MSstatsMakeAnnotation MSstatsSaveSessionInfo
#' MSstatsLogsSettings
#' 
#' @param fewMeasurements 'remove'(default) will remove the features that have 1 or 2 measurements across runs.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned for more than one proteins. 
#' We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have only 1 feature, which is the combination of peptide, precursor charge, fragment and charge. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.
#' @param removeOxidationMpeptides TRUE will remove the peptides including 'oxidation (M)' in modification. FALSE is default.
#' @param removeMpeptides TRUE will remove the peptides including 'M' sequence. FALSE is default.
#' @param use_log_file logical. If TRUE, information about data processing
#' will be saved to a file.
#' @param append logical. If TRUE, information about data processing will be added
#' to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be printed
#' to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. 
#' If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' @param sesion_info_path character. Optional path to a file to which session 
#' information will be saved.
#' 
#' @keywords internal
#' 
.documentFunction = function() {
    NULL
}

standard_columns = c("ProteinName", "PeptideSequence", 
                     "PeptideModifiedSequence", "PrecursorCharge", 
                     "FragmentIon", "ProductCharge", "IsotopeLabelType",
                     "Condition", "BioReplicate", "Run", "StandardType", 
                     "Fraction", "DetectionQValue", "Intensity")

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
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo(session_info_path, append = TRUE)
    
    input = MSstatsImport(list(Fragments = raw.frag, 
                               Peptides = raw.pep, 
                               Proteins = raw.pro), 
                          type = "MSstats", tool = "DIAUmpire", ...)
    input = MSstatsClean(input, 
                         use_frag = useSelectedFrag, 
                         use_pept = useSelectedPep)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "FragmentIon")
    input = MSstatsPreprocess(
        input, 
        annotation,
        feature_columns,
        remove_shared_peptides = TRUE, 
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("PrecursorCharge" = NA,
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsBalancedDesign(input, feature_columns)
    input[, intersect(standard_columns, colnames(input))]
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
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo()
    
    input = MSstatsImport(list(evidence = evidence, 
                               protein_groups = proteinGroups), 
                          type = "MSstats", tool = "MaxQuant", ...)
    input = MSstatsClean(input, 
                         protein_id_col = proteinID, 
                         remove_by_site = TRUE)
    annotation = MSstatsMakeAnnotation(input, 
                                       annotation, 
                                       "Run" = "Rawfile")
    
    m_filter = list(col_name = "PeptideSequence", 
                    pattern = "M", 
                    filter = removeMpeptides, 
                    drop_column = FALSE)
    
    oxidation_filter = list(col_name = "Modifications", 
                            pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, 
                            drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsPreprocess(
        input, 
        annotation,
        feature_columns,
        remove_shared_peptides = useUniquePeptide, 
        remove_single_feature_proteins = removeProtein_with1Peptide,
        pattern_filtering = list(oxidation = oxidation_filter,
                                 m = m_filter),
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("FragmentIon" = NA,
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsBalancedDesign(input, feature_columns)
    input[, intersect(standard_columns, colnames(input))]
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
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo()
    
    input = MSstatsImport(list(input = input), 
                          "MSstats", "OpenMS", ...)
    input = MSstatsClean(input)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", 
                        "FragmentIon", "ProductCharge")
    input = MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows))
    input = MSstatsBalancedDesign(input, feature_columns)
    input[, intersect(standard_columns, colnames(input))]
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
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo()
    
    input = MSstatsImport(list(input = input), 
                          "MSstats", "OpenSWATH", ...)
    input = MSstatsClean(input)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    m_score_filter = list(score_column = "m_score", 
                          score_threshold = mscore_cutoff, 
                          direction = "smaller", 
                          behavior = "remove", 
                          handle_na = "remove", 
                          fill_value = NA,
                          filter = filter_with_mscore, 
                          drop_column = TRUE)
    
    decoy_filter = list(col_name = "decoy", 
                        filter_symbols = 1, 
                        behavior = "remove",
                        fill_value = NULL,
                        filter = TRUE, 
                        drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon")
    input = MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(ms_filter = m_score_filter),
        exact_filtering = list(decoy = decoy_filter),
        columns_to_fill = c("ProductCharge" = NA, 
                            "IsotopeLabelType" = "L"))
    input = MSstatsBalancedDesign(input, feature_columns, 
                                  fix_missing = "na_to_zero")
    input[, intersect(standard_columns, colnames(input))]
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
    removeProtein_with1Peptide = FALSE, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo()
    
    input = MSstatsImport(list(input = input), 
                          "MSstats", "Progenesis", ...)
    input = MSstatsClean(input, 
                         unique(as.character(annotation$Run)), 
                         fix_colnames = TRUE)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    oxidation_filter = list(col_name = "PeptideSequence", 
                            pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        pattern_filtering = list(oxidation = oxidation_filter),
        columns_to_fill = list("FragmentIon" = NA, 
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsBalancedDesign(input, feature_columns)
    data.table::setnames(input, "PeptideSequence", "PeptideModifiedSequence",
                         skip_absent = TRUE)
    input[, intersect(standard_columns, colnames(input))]
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
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo()
    
    input = MSstatsImport(list(input = input), 
                          "MSstats", "ProteomeDiscoverer", ...)
    input = MSstatsClean(input, 
                         quantification_column = which.quantification, 
                         protein_id_column = which.proteinid,
                         sequence_column = which.sequence, 
                         remove_shared = useNumProteinsColumn)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    oxidation_filter = list(col_name = "PeptideSequence", 
                            pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        pattern_filtering = list(oxidation = oxidation_filter),
        columns_to_fill = list("FragmentIon" = NA, 
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsBalancedDesign(input, feature_columns)
    data.table::setnames(input, "PeptideSequence", "PeptideModifiedSequence",
                         skip_absent = TRUE)
    input[, intersect(standard_columns, colnames(input))]
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
    removeOxidationMpeptides = FALSE, removeProtein_with1Feature = FALSE,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo()
    
    input = MSstatsImport(list(input = input), 
                          "MSstats", "Skyline", ...)
    input = MSstatsClean(input)
    annotation = MSstatsMakeAnnotation(input, 
                                       annotation, 
                                       Run = "FileName")
    
    decoy_filter = list(col_name = "ProteinName",
                        pattern = c("DECOY", "Decoys"),
                        filter = TRUE, 
                        drop_column = FALSE)
    
    irt_filter = list(col_name = "StandardType", 
                      filter_symbols = "iRT",
                      filter = removeiRT, 
                      behavior = "remove",
                      fill_value = NULL,
                      drop_column = FALSE)
    
    oxidation_filter = list(col_name = "PeptideSequence",
                            pattern = "\\+16", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    
    truncated_filter = list(col_name = "Truncated", 
                            filter_symbols = "TRUE",
                            behavior = "fill",
                            fill_value = NA_real_,
                            filter = TRUE, 
                            drop_column = TRUE)
    
    qval_filter = list(score_column = "DetectionQValue", 
                       score_threshold = qvalue_cutoff, 
                       direction = "smaller",
                       behavior = "fill", 
                       fill_value = 0, 
                       handle_na = "keep",
                       filter = filter_with_Qvalue, 
                       drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")
    input = MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        score_filtering = list(qval = qval_filter),
        pattern_filtering = list(decoy = decoy_filter, 
                                 oxidation = oxidation_filter),
        exact_filtering = list(irt = irt_filter,
                               truncated = truncated_filter),
        aggregate_isotopic = TRUE,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = sum))
    input = MSstatsBalancedDesign(input, feature_columns)
    input[, intersect(standard_columns, colnames(input))]
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
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path)
    MSstatsSaveSessionInfo()
    
    input = MSstatsImport(list(input = input), 
                          "MSstats", "Spectronaut", ...)
    input = MSstatsClean(input, intensity = intensity)
    annotation = MSstatsMakeAnnotation(input, 
                                       annotation, 
                                       "Run" = "RFileName", 
                                       "Condition" = "RCondition", 
                                       "BioReplicate" = "RReplicate")
    
    pq_filter = list(score_column = "PGQvalue", 
                     score_threshold = 0.01, 
                     direction = "smaller", 
                     behavior = "fill", 
                     handle_na = "keep", 
                     fill_value = NA,
                     filter = TRUE, 
                     drop_column = TRUE)
    qval_filter = list(score_column = "Qvalue", 
                       score_threshold = qvalue_cutoff, 
                       direction = "smaller", 
                       behavior = "fill", 
                       handle_na = "keep", 
                       fill_value = 0, 
                       filter = filter_with_Qvalue, 
                       drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")
    input = MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(handle_features_with_few_measurements = fewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(pgq = pq_filter, 
                               psm_q = qval_filter),
        columns_to_fill = list("IsotopeLabelType" = "L"))
    input = MSstatsBalancedDesign(input, feature_columns)
    input[, intersect(standard_columns, colnames(input))]
}
