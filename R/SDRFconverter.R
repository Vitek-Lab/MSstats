#' Convert SDRF experimental design file into an MSstats annotation file
#' 
#' Takes an SDRF file and outputs an MSstats annotation file. Note
#' the information in the SDRF file must be correctly annotated for MSstats so 
#' that MSstats can identify the experimental design. In particular the 
#' biological replicates must be correctly annotated, with group comparison 
#' experiments having a unique ID for each BioReplicate. For more information 
#' on this please see the Supplementary of the most recent  
#' \href{https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00834}{MSstats paper}
#' 
#' @param data SDRF annotation file
#' @param run_name Column name in SDRF file which contains the name of the MS 
#' run. The information in this column must match exactly with the run names in
#' the PSM file
#' @param condition_name Column name in SDRF file which contains information on 
#' the conditions in the data.
#' @param biological_replicate Column name in SDRF file which contains the 
#' identifier for the biological replicte. Note MSstats uses this column to 
#' determine if the experiment is a repeated measure design. BioReplicte IDs 
#' should only be reused if the replicate was measured multiple times.
#' @param fraction Column name in SDFT file which contains information on the 
#' fractionation in the data. Only required if data contains fractions. Default 
#' is `NULL`
#' 
#' @importFrom data.table setDT
#' 
#' @export
#' 
#' @examples
#' head(example_SDRF)
#' 
#' msstats_annotation = SDRFtoAnnotation(example_SDRF)
#' 
#' head(msstats_annotation)
SDRFtoAnnotation = function(
        data, 
        run_name = "comment[data file]",
        condition_name = "characteristics[disease]",
        biological_replicate = "characteristics[biological replicate]",
        fraction = NULL){
    
    data = data.table::setDT(data)
    
    extract_cols = c(run_name, condition_name, biological_replicate)
    if (!is.null(fraction)){
        extract_cols = c(extract_cols, fraction)
    }
    colnames(data) = MSstatsConvert:::.standardizeColnames(colnames(data))
    extract_cols = MSstatsConvert:::.standardizeColnames(extract_cols)
    
    data = data[, ..extract_cols]
    if (length(colnames(data)) < length(extract_cols)){
        stop("ERROR: One or more of the column passed in the parameters were not found in the data. Please ensure that the column names are correct.")
    }
    data.table::setnames(data, extract_cols, 
                         c("Run", "Condition", "BioReplicate"))
    
    return(data)
}

#' Extract experimental design from MSstats format into SDRF format
#' 
#' @param data MSstats formatted data that is the output of a dedicated 
#' converter, such as `MaxQtoMSstatsFormat`, `SkylinetoMSstatsFormat`, ect.
#' @param run_name Run column name in SDRF data
#' @param condition_name Condition column name in SDRF data
#' @param biological_replicate Biological replicate column name in SDRF data
#' @param fraction Fraction column name in SDRF data (if applicable). Default is
#' `NULL`. If there are no fractions keep `NULL`.
#' @param meta_data A data.frame including any additional meta data for the SDRF 
#' file that is not included in MSstats. This meta data will be added into the 
#' final SDRF file. Please ensure the run names in the meta data matches the 
#' run names in the MSstats data.
#' 
#' @importFrom data.table as.data.table
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
#' SDRF_file = extractSDRF(maxq_imported)
extractSDRF = function(
        data, 
        run_name = "comment[data file]",
        condition_name = "characteristics[disease]",
        biological_replicate = "characteristics[biological replicate]",
        fraction = NULL,
        meta_data = NULL){
    
    extract_cols = c("Condition", "BioReplicate", "Run", "Fraction")
    data = as.data.table(data)
    data = data[, ..extract_cols]
    data = unique(data)
    
    if (is.null(fraction)){
        data$Fraction = NULL
        data.table::setnames(data, c("Condition", "BioReplicate", "Run"), 
                 c(run_name, condition_name, biological_replicate))
    } else {
        data.table::setnames(data, extract_cols, 
                 c(run_name, condition_name, biological_replicate, fraction))
    }
    
    if (!is.null(meta_data)){
        meta_data = data.table::setDT(meta_data)
        data = merge(data, meta_data, all.x = TRUE, all.y = TRUE, by = run_name)
    }
    
    return(data)
}