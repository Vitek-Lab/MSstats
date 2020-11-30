#!/usr/bin/env Rscript

library(data.table)

source("/home/rstudio/code/deployment/tests/utils/s3_helper_functions.R")
source("/home/rstudio/code/deployment/tests/utils/dataprocess_helper_functions.R")
source("/home/rstudio/code/deployment/tests/utils/constants.R")
source("/home/rstudio/code/deployment/tests/utils/generic_utils.R")

######################## Navarro DIAUmpire dataset ############################

# `load()` R objects from S3
# frags <- get_file_from_s3(s3_file_path = datasets$navarro_diaumpire_frags,
#   local_file_name = "frags.RDS")
# pepts <- get_file_from_s3(s3_file_path = datasets$navarro_diaumpire_pepts,
#   local_file_name = "pepts.RDS")
# prots <- get_file_from_s3(s3_file_path = datasets$navarro_diaumpire_prots,
#   local_file_name = "prots.RDS")
# dia_annot <- get_file_from_s3(s3_file_path = datasets$navarro_diaumpire_dia_annot,
#   local_file_name = "dia_annot.RDS")

###############################################################################



############################### Navarro OpenSWATH ############################
# `load()` R objects from S3
open_swath_input <- get_file_from_s3(s3_file_path = datasets$navarro_openswath_input,
  local_file_name = "open_swath_input.RDS")
open_swath_annot <- get_file_from_s3(s3_file_path = datasets$navarro_openswath_annot,
  local_file_name = "open_swath_annot.RDS")

#converter
open_swath = MSstatsdev::OpenSWATHtoMSstatsFormat(open_swath_input, open_swath_annot)

#runs dataprocess function on stable and dev versions
navarro_open_swath_output <- invoke_dataprocess(open_swath)
navarro_open_swath_output_v3 <- navarro_open_swath_output$stable
navarro_open_swath_output_v4 <- navarro_open_swath_output$dev

# checking output on processed data
navarro_open_swath_processed_data_v3 <- as.data.table(navarro_open_swath_output_v3$ProcessedData)
navarro_open_swath_processed_data_v4 <- as.data.table(navarro_open_swath_output_v4$ProcessedData)

#check columns that have null values
navarro_open_swath_processed_data_v3[, which(colnames(
  navarro_open_swath_processed_data_v3) %in% c('GROUP', 'SUBJECT_NESTED', 'SUBJECT')):=NULL]
#rename data frame column names
navarro_open_swath_processed_data_v3 <- rename_column(navarro_open_swath_processed_data_v3, 
              old_names = c('GROUP_ORIGINAL', 'SUBJECT_ORIGINAL'), 
              new_names = c('GROUP', 'SUBJECT'))

#merge two dataframes
navarro_open_swath_compare_processed <- merge_dataframes(
  df1 = navarro_open_swath_processed_data_v3, 
  df2 = navarro_open_swath_processed_data_v4,
  col_names = c("GROUP", "SUBJECT", "INTENSITY", "ABUNDANCE", "censored", 
                'predicted', 'remove', 'feature_quality', 'is_outlier'))


## flag the difference : TRUE - matched, FALSE - issue
navarro_open_swath_compare_processed <- compare_values_processed_data(
  navarro_open_swath_compare_processed)

## if both NA, they also match.
navarro_open_swath_compare_processed <- handle_na_values_processed_data(
  navarro_open_swath_compare_processed
)

## report this number
navarro.processed.report <- data.frame(
  dataset = "Navarro OpenSWATH",
  s3_dataset_path = datasets$navarro_openswath_input,
  args = "linear summary method",
  notes = "Processed data",
  match.orgint = sum(!navarro_open_swath_compare_processed$match.orgint),
  match.abn = sum(!navarro_open_swath_compare_processed$match.abn),
  match.censored = sum(!navarro_open_swath_compare_processed$match.censored),
  match.ftrslct = sum(!navarro_open_swath_compare_processed$match.ftrslct),
  match.otr = sum(!navarro_open_swath_compare_processed$match.otr)
  )

# checking output on runlevel data
navarro_open_swath_runlevel_data_v3 <- as.data.table(navarro_open_swath_output_v3$RunlevelData)
navarro_open_swath_runlevel_data_v4 = as.data.table(navarro_open_swath_output_v4$RunlevelData)

navarro_open_swath_runlevel_data_v3 <- navarro_open_swath_runlevel_data_v3[, -which(colnames(navarro_open_swath_runlevel_data_v3) %in% c('GROUP', 'SUBJECT_NESTED', 'SUBJECT'))]

navarro_open_swath_runlevel_data_v4 <- as.data.table(nanavarro_open_swath_runlevel_data_v4)
navarro_open_swath_runlevel_data_v4[, GROUP := as.factor(as.character(GROUP))]
navarro_open_swath_runlevel_data_v4[, SUBJECT := as.factor(as.character(SUBJECT))]

navarro_open_swath_compare_runlevel_summarized = merge_dataframes(
  df1 = navarro_open_swath_runlevel_data_v3,
  df2 = navarro_open_swath_runlevel_data_v4,
  col_names = c("RUN", "Protein")
)

## flag the difference : TRUE - matched, FALSE - issue
navarro_open_swath_compare_runlevel_summarized <- compare_values_runlevel_data(
  navarro_open_swath_compare_runlevel_summarized)

## if both NA, they also match.
navarro_open_swath_compare_runlevel_summarized <- handle_na_values_runlevel_data(
  navarro_open_swath_compare_runlevel_summarized)

## report this number
navarro.runlevel.report <- data.frame(
  dataset = "Navarro OpenSWATH",
  s3_dataset_path = datasets$navarro_openswath_input,
  args = "linear summary method",
  notes = "RunLevel data",
  match.int = sum(!navarro_open_swath_compare_runlevel_summarized$match.int),
  match.perc = sum(!navarro_open_swath_compare_runlevel_summarized$match.perc),
  match.numimpf = sum(!navarro_open_swath_compare_runlevel_summarized$match.numimpf),
  match.nummsrf = sum(!navarro_open_swath_compare_runlevel_summarized$match.nummsrf)
)
#uploading the results to s3 as csv
store_csv_file_to_s3(s3_path = results$navarro_openswath_results, 
                 local_file_name = "report.xlsx", upload_file=list(
                   navarro.processed.report, navarro.runlevel.report))
################################################################################


####### ##################### cleaning up #####################################
closeAllConnections()
rm(list=ls(all=TRUE))
###############################################################################
