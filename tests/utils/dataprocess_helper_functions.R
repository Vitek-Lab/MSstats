invoke_dataprocess_feature_subset_all<-function(input_file, 
                                                summary_method="TMP", 
                                                mb_impute=FALSE,
                                                censored_int="0",
                                                feature_subset="all"){
  #
  #against stable version
  output_stable = MSstats::dataProcess(as.data.frame(unclass(input_file)), 
                                       summaryMethod=summary_method, 
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int,
                                       featureSubset = feature_subset)
  #against dev version
  output_dev = MSstatsdev::dataProcess(input_file, 
                                       summaryMethod=summary_method, 
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int,
                                       featureSubset = feature_subset)
  return(list(stable=output_stable, dev=output_dev))
}

invoke_dataprocess_feature_subset_topn<-function(input_file, 
                                                 summary_method="TMP", 
                                                 mb_impute=FALSE,
                                                 censored_int="0",
                                                 feature_subset="topN",
                                                 n_top_feature=5){
  # data process function parameterized on 'featureSubset' argument to 'topN' 
  # and 'n_top_feature' is defaulted to 5
  #against stable version
  output_stable = MSstats::dataProcess(as.data.frame(unclass(input_file)), 
                                       summaryMethod=summary_method, 
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int,
                                       n_top_feature = n_top_feature,
                                       featureSubset = feature_subset)
  #against dev version
  output_dev = MSstatsdev::dataProcess(input_file, 
                                       summaryMethod=summary_method, 
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int,
                                       n_top_feature = n_top_feature,
                                       featureSubset = feature_subset)
  return(list(stable=output_stable, dev=output_dev))
}

invoke_dataprocess_feature_subset_high_quality<-function(input_file, 
                                                         summary_method="TMP", 
                                                         mb_impute=FALSE, 
                                                         censored_int="0",
                                                         feature_subset = "highQuality",
                                                         remove_uninformative_feature_outlier=TRUE){
  # dataprocess function parameterized on 'featureSubset' argument to 'highQuality'
  # remove_uninformative_feature_outlier is defaulted to TRUE
  #against stable version
  output_stable = MSstats::dataProcess(as.data.frame(unclass(input_file)), 
                                       summaryMethod=summary_method,
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int,
                                       featureSubset = feature_subset,
                                       remove_uninformative_feature_outlier=remove_uninformative_feature_outlier)
  #against dev version
  output_dev = MSstatsdev::dataProcess(input_file, 
                                       summaryMethod=summary_method,
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int,
                                       featureSubset = feature_subset,
                                       remove_uninformative_feature_outlier=remove_uninformative_feature_outlier)
  return(list(stable=output_stable, dev=output_dev))
}


compare_values_processed_data <- function(input_df){
  #this logic is overkilled for now. Need to handle columns more flexibly
  input_df$match.orgint <- input_df$INTENSITY.x == input_df$INTENSITY.y
  input_df$match.abn <- input_df$ABUNDANCE.x == input_df$ABUNDANCE.y
  input_df$match.censored <- input_df$censored.x == input_df$censored.y
  input_df$match.ftrslct <- input_df$feature_quality.x == input_df$feature_quality.y
  input_df$match.otr <- input_df$is_outlier.x == input_df$is_outlier.y
  return(input_df)
}

compare_values_runlevel_data <- function(input_df){
  #this logic is overkilled for now. Need to handle columns more flexibly
  input_df$match.int <- input_df$LogIntensities.x == input_df$LogIntensities.y
  input_df$match.perc <- input_df$MissingPercentage.x == input_df$MissingPercentage.y
  input_df$match.numimpf <- input_df$NumImputedFeature.x == input_df$NumImputedFeature.y
  input_df$match.nummsrf <- input_df$NumMeasuredFeature.x == input_df$NumMeasuredFeature.y
  return(input_df)
}

handle_na_values_processed_data <- function(input_df ){
  ## if both NA, they also match.
  input_df[is.na(
    input_df$INTENSITY.x) & is.na(input_df$INTENSITY.y), 'match.orgint'] <- TRUE
  input_df[is.na(
    input_df$ABUNDANCE.x) & is.na(input_df$ABUNDANCE.y), 'match.abn'] <- TRUE
  input_df[is.na(
    input_df$censored.x) & is.na(input_df$censored.y), 'match.censored'] <- TRUE
  input_df[is.na(
    input_df$feature_quality.x) & is.na(input_df$feature_quality.y), 'match.ftrslct'] <- TRUE
  input_df[is.na(
    input_df$is_outlier.x) & is.na(input_df$is_outlier.y), 'match.otr'] <- TRUE
  
  ## others, should be FALSE : if one of them NA, match column NA -> change to FALSE
  input_df[is.na(input_df$match.orgint), 'match.orgint'] <- FALSE
  input_df[is.na(input_df$match.abn), 'match.abn'] <- FALSE
  input_df[is.na(input_df$match.censored), 'match.censored'] <- FALSE
  input_df[is.na(input_df$match.ftrslct), 'match.ftrslct'] <- FALSE
  input_df[is.na(input_df$match.otr), 'match.otr'] <- FALSE
  return(input_df)
}

handle_na_values_runlevel_data <- function(input_df){
  input_df[is.na(input_df$LogIntensities.x) & is.na(
    input_df$LogIntensities.y), 'match.int'] <- TRUE
  input_df[is.na(input_df$MissingPercentage.x) & is.na(
    input_df$MissingPercentage.y), 'match.perc'] <- TRUE
  input_df[is.na(input_df$NumImputedFeature.x) & is.na(
    input_df$NumImputedFeature.y), 'match.numimpf'] <- TRUE
  input_df[is.na(input_df$NumMeasuredFeature.x) & is.na(
    input_df$NumMeasuredFeature.y), 'match.nummsrf'] <- TRUE
  
  ## others, should be FALSE : if one of them NA, match column NA -> change to FALSE
  input_df[is.na(input_df$match.int), 'match.int'] <- FALSE
  input_df[is.na(input_df$match.perc), 'match.perc'] <- FALSE
  input_df[is.na(input_df$match.numimpf), 'match.numimpf'] <- FALSE
  input_df[is.na(input_df$match.nummsrf), 'match.nummsrf'] <- FALSE
  return(input_df)
}

run_comparisons <- function(dataprocess_output, master_df, notes, summary_method, dataset_path){
  dataprocess_output_v3 <- dataprocess_output$stable
  dataprocess_output_v4 <- dataprocess_output$dev
  
  # checking output on processed data
  processed_data_v3 <- as.data.table(dataprocess_output_v3$ProcessedData)
  processed_data_v4 <- as.data.table(dataprocess_output_v4$ProcessedData)
  
  #check columns that have null values
  processed_data_v3[, which(colnames(
    processed_data_v3) %in% c('GROUP', 'SUBJECT_NESTED', 'SUBJECT')):=NULL]
  #rename data frame column names
  processed_data_v3 <- rename_column(
    processed_data_v3, old_names = c('GROUP_ORIGINAL', 'SUBJECT_ORIGINAL'), 
    new_names = c('GROUP', 'SUBJECT'))
  
  #merge two dataframes
  compare_processed <- merge(
    processed_data_v3, processed_data_v4, by = setdiff(
      colnames(processed_data_v3),c("GROUP", "SUBJECT", "INTENSITY", "ABUNDANCE", 
                                    "censored", 'predicted', 'remove', 
                                    'feature_quality', 'is_outlier')),
    all.x = T, all.y = T)
  
  
  ## flag the difference : TRUE - matched, FALSE - issue
  compare_processed <- compare_values_processed_data(compare_processed)
  
  ## if both NA, they also match.
  compare_processed <- handle_na_values_processed_data(compare_processed)
  
  ## report this number
  processed.report <- data.frame(
    summary_method_type = summary_method,
    data_type = "Processed data",
    parameters=notes,
    mismatched_intensity = sum(!compare_processed$match.orgint),
    mismatched_abundance = sum(!compare_processed$match.abn),
    mismatched_censored = sum(!compare_processed$match.censored),
    mismatched_feature_quality = sum(!compare_processed$match.ftrslct),
    mismatched_outlier = sum(!compare_processed$match.otr),
    s3_dataset_path = dataset_path
  )
  
  master_df$master_processed_data <-rbind(master_df$master_processed_data, 
                                          processed.report)
  # checking output on runlevel data
  runlevel_data_v3 <- as.data.table(dataprocess_output_v3$RunlevelData)
  runlevel_data_v4 = as.data.table(dataprocess_output_v4$RunlevelData)
  
  runlevel_data_v3[, which(colnames(runlevel_data_v3) %in% c('GROUP', 'SUBJECT_NESTED', 'SUBJECT')):=NULL]
  runlevel_data_v4[, GROUP := as.factor(as.character(GROUP))]
  runlevel_data_v4[, SUBJECT := as.factor(as.character(SUBJECT))]
  
  summarized_runlevel <- merge(
    runlevel_data_v3, runlevel_data_v4, by = c("RUN", "Protein"),
    all.x = T, all.y = T
  )
  
  ## flag the difference : TRUE - matched, FALSE - issue
  summarized_runlevel <- compare_values_runlevel_data(
    summarized_runlevel)
  
  ## if both NA, they also match.
  summarized_runlevel <- handle_na_values_runlevel_data(
    summarized_runlevel)
  
  ## report this number
  runlevel.report <- data.frame(
    summary_method_type = summary_method,
    data_type = "RunLevel data",
    parameters=notes,
    mismatched_log_intensity = sum(!summarized_runlevel$match.int),
    mismatched_missing_percentage = sum(!summarized_runlevel$match.perc),
    mismatched_imputed_feature_count = sum(!summarized_runlevel$match.numimpf),
    mismatched_measured_feature_count = sum(!summarized_runlevel$match.nummsrf),
    s3_dataset_path = dataset_path
  )
  master_df$master_run_level_data <-rbind(master_df$master_run_level_data, 
                                          runlevel.report)
  return(master_df)
}