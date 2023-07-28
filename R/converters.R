#' A dummy function to store shared documentation items.
#' 
#' @import data.table
#' @importFrom MSstatsConvert MSstatsImport MSstatsClean MSstatsPreprocess 
#' MSstatsBalancedDesign MSstatsMakeAnnotation MSstatsLogsSettings
#' 
#' @param removeFewMeasurements TRUE (default) will remove the features that have 1 or 2 measurements across runs.
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
#' 
#' @keywords internal
#' 
.documentFunction = function() {
    NULL
}


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
#' @return data.frame in the MSstats required format.
#'
#' @author Meena Choi, Olga Vitek 
#'
#' @export
#' 
#' @examples 
#' diau_frag = system.file("tinytest/raw_data/DIAUmpire/dia_frag.csv", 
#'                              package = "MSstatsConvert")
#' diau_pept = system.file("tinytest/raw_data/DIAUmpire/dia_pept.csv", 
#'                              package = "MSstatsConvert")
#' diau_prot = system.file("tinytest/raw_data/DIAUmpire/dia_prot.csv", 
#'                              package = "MSstatsConvert")
#' annot = system.file("tinytest/annotations/annot_diau.csv", 
#'                     package = "MSstats")
#' diau_frag = data.table::fread(diau_frag) 
#' diau_pept = data.table::fread(diau_pept) 
#' diau_prot = data.table::fread(diau_prot) 
#' annot = data.table::fread(annot)
#' diau_frag = diau_frag[, lapply(.SD, function(x) if (is.integer(x)) as.numeric(x) else x)]
#' # In case numeric columns are not interpreted correctly
#' 
#' diau_imported = DIAUmpiretoMSstatsFormat(diau_frag, diau_pept, diau_prot, 
#'                                          annot, use_log_file = FALSE)
#' head(diau_imported)
#' 
DIAUmpiretoMSstatsFormat = function(
    raw.frag, raw.pep, raw.pro, annotation, useSelectedFrag = TRUE,
    useSelectedPep = TRUE, removeFewMeasurements = TRUE,
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(Fragments = raw.frag, 
                                               Peptides = raw.pep, 
                                               Proteins = raw.pro), 
                                          type = "MSstats", 
                                          tool = "DIAUmpire", ...)
    input = MSstatsConvert::MSstatsClean(input, 
                                         use_frag = useSelectedFrag, 
                                         use_pept = useSelectedPep)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "FragmentIon")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation,
        feature_columns,
        remove_shared_peptides = TRUE, 
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("PrecursorCharge" = NA,
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
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
#' @return data.frame in the MSstats required format.
#'  
#' @note Warning: MSstats does not support for metabolic labeling or iTRAQ experiments.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
#' @examples 
#' mq_ev = data.table::fread(system.file("tinytest/raw_data/MaxQuant/mq_ev.csv",
#'                                       package = "MSstatsConvert"))
#' mq_pg = data.table::fread(system.file("tinytest/raw_data/MaxQuant/mq_pg.csv",
#'                                       package = "MSstatsConvert"))
#' annot = data.table::fread(system.file("tinytest/raw_data/MaxQuant/annotation.csv",
#'                                       package = "MSstatsConvert"))
#' maxq_imported = MaxQtoMSstatsFormat(mq_ev, annot, mq_pg, use_log_file = FALSE)
#' head(maxq_imported)
#' 
MaxQtoMSstatsFormat = function(
    evidence, annotation, proteinGroups, proteinID = "Proteins", 
    useUniquePeptide = TRUE, summaryforMultipleRows = max, 
    removeFewMeasurements = TRUE, removeMpeptides = FALSE,
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(evidence = evidence, 
                                               protein_groups = proteinGroups), 
                                          type = "MSstats",
                                          tool = "MaxQuant", ...)
    input = MSstatsConvert::MSstatsClean(input, 
                                         protein_id_col = proteinID, 
                                         remove_by_site = TRUE)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, 
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
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation,
        feature_columns,
        remove_shared_peptides = useUniquePeptide, 
        remove_single_feature_proteins = removeProtein_with1Peptide,
        pattern_filtering = list(oxidation = oxidation_filter,
                                 m = m_filter),
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("FragmentIon" = NA,
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
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
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
#' 
#' @examples 
#' openms_raw = data.table::fread(system.file("tinytest/raw_data/OpenMS/openms_input.csv", 
#'                                            package = "MSstatsConvert"))
#' openms_imported = OpenMStoMSstatsFormat(openms_raw, use_log_file = FALSE)
#' head(openms_imported)
#' 
OpenMStoMSstatsFormat = function(
    input, annotation = NULL, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "OpenMS", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", 
                        "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
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
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
#' @examples 
#' os_raw = system.file("tinytest/raw_data/OpenSWATH/openswath_input.csv", 
#'                              package = "MSstatsConvert")
#' annot = system.file("tinytest/annotations/annot_os.csv", 
#'                     package = "MSstats")
#' os_raw = data.table::fread(os_raw) 
#' annot = data.table::fread(annot)
#' 
#' os_imported = OpenSWATHtoMSstatsFormat(os_raw, annot, use_log_file = FALSE)
#' head(os_imported)
#' 
OpenSWATHtoMSstatsFormat = function(
    input, annotation, filter_with_mscore = TRUE, mscore_cutoff = 0.01,
    useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "OpenSWATH", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
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
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(ms_filter = m_score_filter),
        exact_filtering = list(decoy = decoy_filter),
        columns_to_fill = c("ProductCharge" = NA, 
                            "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns, 
                                                  fix_missing = "na_to_zero",
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}


#' Import Diann files
#' 
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from Diann, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' @param MBR True if analysis was done with match between runs
#' @param global_qvalue_cutoff The global qvalue cutoff
#' @param qvalue_cutoff local qvalue cutoff for library
#' @param pg_qvalue_cutoff local qvalue cutoff for protein groups Run should be the same as filename.
#' @param useUniquePeptide should unique pepties be removed
#' @param removeFewMeasurements should proteins with few measurements be removed
#' @param removeOxidationMpeptides should peptides with oxidation be removed
#' @param removeProtein_with1Feature should proteins with a single feature be removed
#' @param ... additional parameters to `data.table::fread`.
#'  
#' @return data.frame in the MSstats required format.
#' 
#' @author Elijah Willie
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' input = fread('diann_pooled_report.tsv')
#' annot = fread('Annotation.csv')
#' colnames(annot) = c('Condition', 'Run', 'BioReplicate')
#' input = DIANNtoMSstatsFormat(input, annotation = annot, MBR = F)
#' head(input)
#' }
DIANNtoMSstatsFormat = function(input, annotation = NULL,
                                 global_qvalue_cutoff = 0.01,
                                 qvalue_cutoff = 0.01, 
                                 pg_qvalue_cutoff = 0.01,
                                 useUniquePeptide = TRUE, 
                                 removeFewMeasurements = TRUE,
                                 removeOxidationMpeptides = TRUE, 
                                 removeProtein_with1Feature = TRUE,
                                 use_log_file = TRUE, append = FALSE, 
                                 verbose = TRUE, log_file_path = NULL,
                                 MBR = TRUE,...) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path)
  
  input = MSstatsConvert::MSstatsImport(list(input = input),
                                        "MSstats", "DIANN")
  input = MSstatsConvert::MSstatsClean(input, MBR = MBR)
  annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)

  decoy_filter = list(col_name = "ProteinName",
                      pattern = c("DECOY", "Decoys"),
                      filter = T, 
                      drop_column = FALSE)
  oxidation_filter = list(col_name = "PeptideSequence",
                          pattern = "\\(UniMod\\:35\\)",
                          filter = removeOxidationMpeptides,
                          drop_column = FALSE)
  
  msg = paste0('** Filtering on Global Q Value < ', global_qvalue_cutoff)
  getOption("MSstatsLog")("INFO", msg)
  getOption("MSstatsMsg")("INFO", msg)
  
  input = input[DetectionQValue < global_qvalue_cutoff, ]
  if (MBR) {
    msg = '** MBR was used to analyze the data. Now setting names and filtering'
    msg_1_mbr = paste0('-- LibPGQValue < ', pg_qvalue_cutoff)
    msg_2_mbr = paste0('-- LibQValue < ', qvalue_cutoff)
    input = input[LibPGQValue < pg_qvalue_cutoff, ]
    input = input[LibQValue < qvalue_cutoff, ]
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    getOption("MSstatsLog")("INFO", msg_1_mbr)
    getOption("MSstatsMsg")("INFO", msg_1_mbr)
    getOption("MSstatsLog")("INFO", msg_2_mbr)
    getOption("MSstatsMsg")("INFO", msg_2_mbr)
    # getOption("MSstatsLog")("INFO", "\n")
  } else{
    msg = '** MBR was not used to analyze the data. Now setting names and filtering'
    msg_1 = paste0('-- Filtering on GlobalPGQValue < ', pg_qvalue_cutoff)
    msg_2 = paste0('-- Filtering on GlobalQValue < ', qvalue_cutoff)
    input = input[GlobalPGQValue < pg_qvalue_cutoff, ]
    input = input[GlobalQValue < qvalue_cutoff, ]
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    getOption("MSstatsLog")("INFO", msg_1)
    getOption("MSstatsMsg")("INFO", msg_1)
    getOption("MSstatsLog")("INFO", msg_2)
    getOption("MSstatsMsg")("INFO", msg_2)
    # getOption("MSstatsLog")("INFO", "\n")
  }
  
  feature_columns = c("PeptideSequence", "PrecursorCharge",
                      "FragmentIon", "ProductCharge")
  input = MSstatsConvert::MSstatsPreprocess(
    input, 
    annotation, 
    feature_columns,
    remove_shared_peptides = useUniquePeptide,
    remove_single_feature_proteins = removeProtein_with1Feature,
    exact_filtering = NULL,
    pattern_filtering = list(decoy = decoy_filter, 
                             oxidation = oxidation_filter),
    aggregate_isotopic = FALSE,
    feature_cleaning = list(
      remove_features_with_few_measurements = removeFewMeasurements,
      summarize_multiple_psms = max),
    columns_to_fill = list(Fraction = 1,
                           IsotopeLabelType = "Light"))
  
  input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns, 
                                                fill_incomplete = FALSE,
                                                handle_fractions = FALSE,
                                                remove_few = removeFewMeasurements
  )
  
  msg_final = paste("** Finished preprocessing. The dataset is ready",
                    "to be processed by the dataProcess function.")
  getOption("MSstatsLog")("INFO", msg_final)
  getOption("MSstatsMsg")("INFO", msg_final)
  getOption("MSstatsLog")("INFO", "\n")
  input
}




#' Import Progenesis files
#' 
#' @inheritParams .documentFunction
#' @param input name of Progenesis output, which is wide-format. 'Accession', 'Sequence', 'Modification', 'Charge' and one column for each run are required.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, Run information. It will be matched with the column name of input for MS runs.
#' @param ... additional parameters to `data.table::fread`.
#'
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek, Ulrich Omasits
#' 
#' @export
#' 
#' @examples 
#' progenesis_raw = system.file("tinytest/raw_data/Progenesis/progenesis_input.csv", 
#'                              package = "MSstatsConvert")
#' annot = system.file("tinytest/raw_data/Progenesis/progenesis_annot.csv", 
#'                     package = "MSstatsConvert")
#' progenesis_raw = data.table::fread(progenesis_raw) 
#' annot = data.table::fread(annot)
#' 
#' progenesis_imported = ProgenesistoMSstatsFormat(progenesis_raw, annot,
#'                                                 use_log_file = FALSE)
#' head(progenesis_imported)
#' 
ProgenesistoMSstatsFormat = function(
    input, annotation, useUniquePeptide = TRUE, summaryforMultipleRows = max,
    removeFewMeasurements = TRUE, removeOxidationMpeptides = FALSE, 
    removeProtein_with1Peptide = FALSE, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose,
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Progenesis", ...)
    input = MSstatsConvert::MSstatsClean(input, 
                                         unique(as.character(annotation$Run)), 
                                         fix_colnames = TRUE)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    oxidation_filter = list(col_name = "PeptideSequence", 
                            pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        pattern_filtering = list(oxidation = oxidation_filter),
        columns_to_fill = list("FragmentIon" = NA, 
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    data.table::setnames(input, "PeptideSequence", "PeptideModifiedSequence",
                         skip_absent = TRUE)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
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
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
#' @examples 
#' 
#' pd_raw = system.file("tinytest/raw_data/PD/pd_input.csv", 
#'                      package = "MSstatsConvert")
#' annot = system.file("tinytest/annotations/annot_pd.csv", package = "MSstats")
#' pd_raw = data.table::fread(pd_raw)
#' annot = data.table::fread(annot)
#' 
#' pd_imported = PDtoMSstatsFormat(pd_raw, annot, use_log_file = FALSE)
#' head(pd_imported)
#' 
PDtoMSstatsFormat = function(
    input, annotation, useNumProteinsColumn = FALSE, useUniquePeptide = TRUE,
    summaryforMultipleRows = max, removeFewMeasurements = TRUE,
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE,
    which.quantification = 'Precursor.Area', 
    which.proteinid = 'Protein.Group.Accessions', which.sequence = 'Sequence',
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "ProteomeDiscoverer", ...)
    input = MSstatsConvert::MSstatsClean(
        input, 
        quantification_column = which.quantification, 
        protein_id_column = which.proteinid,
        sequence_column = which.sequence, 
        remove_shared = useNumProteinsColumn)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    oxidation_filter = list(col_name = "PeptideSequence", 
                            pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        pattern_filtering = list(oxidation = oxidation_filter),
        columns_to_fill = list("FragmentIon" = NA, 
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    data.table::setnames(input, "PeptideSequence", "PeptideModifiedSequence",
                         skip_absent = TRUE)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}


#' Import Skyline files
#'
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from Skyline, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Skyline, use annotation=NULL (default). It will use the annotation information from input.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which are labeled 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in DetectionQValue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
#' @examples 
#' skyline_raw = system.file("tinytest/raw_data/Skyline/skyline_input.csv",
#'                           package = "MSstatsConvert")
#' skyline_raw = data.table::fread(skyline_raw)
#' skyline_imported = SkylinetoMSstatsFormat(skyline_raw)
#' head(skyline_imported)
#' 
SkylinetoMSstatsFormat = function(
    input, annotation = NULL, removeiRT = TRUE, filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
    removeOxidationMpeptides = FALSE, removeProtein_with1Feature = FALSE,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Skyline", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, 
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
    
    feature_columns = c("IsotopeLabelType", "PeptideSequence", "PrecursorCharge", 
                        "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
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
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = sum))
    input = MSstatsBalancedDesign(input, c("PeptideSequence", "PrecursorCharge", 
                                           "FragmentIon", "ProductCharge"),
                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
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
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
#' @examples 
#' spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
#'                               package = "MSstatsConvert")
#' spectronaut_raw = data.table::fread(spectronaut_raw)
#' spectronaut_imported = SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
#' head(spectronaut_imported)
#' 
SpectronauttoMSstatsFormat = function(
    input, annotation = NULL, intensity = 'PeakArea', filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01, useUniquePeptide = TRUE, removeFewMeasurements=TRUE,
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Spectronaut", ...)
    input = MSstatsConvert::MSstatsClean(input, intensity = intensity)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    pq_filter = list(score_column = "PGQvalue", 
                     score_threshold = 0.01, 
                     direction = "smaller", 
                     behavior = "fill", 
                     handle_na = "keep", 
                     fill_value = NA,
                     filter = TRUE, 
                     drop_column = TRUE)
    qval_filter = list(score_column = "EGQvalue", 
                       score_threshold = qvalue_cutoff, 
                       direction = "smaller", 
                       behavior = "fill", 
                       handle_na = "keep", 
                       fill_value = 0, 
                       filter = filter_with_Qvalue, 
                       drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = removeFewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(pgq = pq_filter, 
                               psm_q = qval_filter),
        columns_to_fill = list("IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}

#' Import FragPipe files
#' 
#' @param input name of FragPipe msstats.csv export. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity are required.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .documentFunction
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Devon Kohler
#' 
#' @export
#' 
#' @examples 
#' fragpipe_raw = system.file("tinytest/raw_data/FragPipe/fragpipe_input.csv",
#'                               package = "MSstatsConvert")
#' fragpipe_raw = data.table::fread(fragpipe_raw)
#' fragpipe_imported = FragPipetoMSstatsFormat(fragpipe_raw, use_log_file = FALSE)
#' head(fragpipe_imported)
#' 
FragPipetoMSstatsFormat = function(
        input, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "FragPipe", ...)
    input = MSstatsConvert::getInputFile(input, "input")
    
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, NULL)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = removeFewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}
