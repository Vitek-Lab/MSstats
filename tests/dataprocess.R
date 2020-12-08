#!/usr/bin/env Rscript

library(data.table)

source("/home/rstudio/code/deployment/tests/utils/s3_helper_functions.R")
source("/home/rstudio/code/deployment/tests/utils/dataprocess_helper_functions.R")
source("/home/rstudio/code/deployment/tests/utils/constants.R")
source("/home/rstudio/code/deployment/tests/utils/generic_utils.R")

######################## get required files ############################

path_to_datasets = "datasets/"
#get metadata from s3
metadata_s3 <- get_file_from_s3(s3_file_path = datasets$metadata_rds,
                                local_file_name = "metadata.RDS")
datasets_s3 <- get_file_from_s3(s3_file_path = datasets$datasets_rds,
                                local_file_name = "datasets.RDS")
metadata_s3 = metadata_s3[!(sapply(metadata_s3, function(x) x$name) %in% c("Azimifa2014"))]

###############################################################################



######################### define dataprocess function ##########################
run_dataprocess <- function(data, 
                            dataset_path, master_result_df, 
                            mb_impute=FALSE,censored_int="0"){
  
  ################# parameterized dataprocess run 1 ####################
  # summaryMethod='TMP' + MBimpute=T or F + censoredInt= 'NA' or '0' + featureSub='All'
  feature_sub_all_dataprocess_output <- invoke_dataprocess_feature_subset_all(
    data, summary_method="TMP", mb_impute, censored_int, feature_subset="all")
  master_result_df <- run_comparisons(feature_sub_all_dataprocess_output,
                                      master_df=master_result_df,
                                      notes="summaryMethod='TMP' + MBimpute=T/F + censoredInt= 'NA'/'0' + featureSub='All'",
                                      summary_method="TMP", dataset_path)
  print("summaryMethod='TMP' + MBimpute=T or F + censoredInt= 'NA' or '0' + featureSub='All'")
  
  ############### parameterized dataprocess run 2 #########################
  # summaryMethod='TMP' + MBimpute=T + censoredInt= 'NA' or '0' + featureSub='topN' + n_top_feature=5
  top_n_dataprocess_output <-invoke_dataprocess_feature_subset_topn(data, summary_method="TMP",
                                                                    mb_impute,
                                                                    censored_int,
                                                                    feature_subset="topN",
                                                                    n_top_feature=5)
  master_result_df <- run_comparisons(top_n_dataprocess_output, master_df=master_result_df,
                                      notes = "summaryMethod='TMP' + MBimpute=T/F + censoredInt= 'NA'/'0' + featureSub='topN' + n_top_feature=5",
                                      summary_method="TMP", dataset_path)
  print("summaryMethod='TMP' + MBimpute=T + censoredInt= 'NA' or '0' + featureSub='topN' + n_top_feature=5")
  
  ############### parameterized dataprocess run 3#########################
  # summaryMethod='TMP' + MBimpute=T/F + censoredInt= 'NA'/'0' + featureSub='highQuality' + remove_uninformative_feature_outlier = T/ F
  hq_dataprocess_output <- invoke_dataprocess_feature_subset_high_quality(
    data, summary_method="TMP", mb_impute, censored_int,
    feature_subset = "highQuality",
    remove_uninformative_feature_outlier=TRUE)
  master_result_df <- run_comparisons(hq_dataprocess_output, master_df=master_result_df,
                                      notes = "summaryMethod='TMP' + MBimpute=T/F + censoredInt= 'NA'/'0' + featureSub='highQuality' + remove_uninformative_feature_outlier = T/ F",
                                      summary_method="TMP", dataset_path)
  print("summaryMethod='TMP' + MBimpute=T/F + censoredInt= 'NA'/'0' + featureSub='highQuality' + remove_uninformative_feature_outlier = T/ F")
  
  ############### parameterized dataprocess run 4#########################
  # summaryMethod='Linear' + MBimpute=F + censoredInt= 'NA' or '0' + featureSub='All'
  linear_feature_sub_all_dataprocess_output <- invoke_dataprocess_feature_subset_all(
    data, summary_method="linear", mb_impute, censored_int, feature_subset="All")
  master_result_df <- run_comparisons(linear_feature_sub_all_dataprocess_output, master_df=master_result_df,
                                      notes = "summaryMethod='Linear' + MBimpute=F + censoredInt= 'NA' or '0' + featureSub='All'",
                                      summary_method="linear", dataset_path)
  print("summaryMethod='Linear' + MBimpute=F + censoredInt= 'NA' or '0' + featureSub='All'")
  
  ############### parameterized dataprocess run 5#########################
  # summaryMethod='Linear' + MBimpute=F + censoredInt= 'NA' or '0' + featureSub='topN'
  linear_feature_sub_all_dataprocess_output_topn <- invoke_dataprocess_feature_subset_topn(
    data, summary_method="linear", mb_impute, censored_int, feature_subset="topN",n_top_feature=5)
  master_result_df <- run_comparisons(linear_feature_sub_all_dataprocess_output_topn, master_df=master_result_df,
                                      notes = "parameterized on remove_uninformative_feature_outlier = FALSE",
                                      summary_method="linear", dataset_path)
  print("summaryMethod='Linear' + MBimpute=F + censoredInt= 'NA' or '0' + featureSub='topN'")
  
  return(master_result_df)
}
################################################################################



######################### WIDER TESTING ########################################

run_wider_testing <- function(metadata,
                              remove_single_feature = FALSE, 
                              remove_few = TRUE) {
  
  remove_few_lf = ifelse(remove_few, "remove", "keep")
  master_processed_data <- data.frame()
  master_run_level_data <- data.frame()
  master_results <- list(master_processed_data,master_run_level_data)
  
  for (dataset in metadata) {
    tryCatch({
      names = names(dataset)
      dataset = lapply(dataset, function(x) gsub(path_to_datasets, "", x))
      names(dataset) = names
      if (dataset$tool == "MaxQuant") {
        #######################################################################
        ################# constructing dataset paths on s3 ###################
        s3_evidence_path <- gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$evidence_path))))
        s3_protein_group_path <-gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$protein_groups))))
        s3_annotation_path <- gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$annotation))))
        ################# Loading object to memory ###################
        evidence <- get_file_from_s3(s3_file_path = s3_evidence_path,
                                     local_file_name = "evidence.RDS")
        protein_groups <- get_file_from_s3(s3_file_path = s3_protein_group_path,
                                           local_file_name = "protein_group.RDS")
        annotation <- get_file_from_s3(s3_file_path = s3_annotation_path,
                                       local_file_name = "annotation.RDS")
        #######################################################################
        if (dataset$type == "TMT"){
          try({
            max_quant_conv = MSstatsTMTdev::MaxQtoMSstatsTMTFormat(
              evidence, protein_groups, annotation,
              rmPSM_withfewMea_withinRun = remove_few,
              rmProtein_with1Feature = remove_single_feature)
          })
        }else{
          try({
            max_quant_conv = MSstatsdev::MaxQtoMSstatsFormat(
              evidence, annotation, protein_groups,
              rmPSM_withfewMea_withinRun = remove_few,
              rmProtein_with1Feature = remove_single_feature)
          })
        }
        try({
          master_results <- run_dataprocess(data=max_quant_conv, 
                                            dataset_path=s3_evidence_path, 
                                            master_result_df = master_results, 
                                            mb_impute=TRUE,
                                            censored_int="NA")
        })
      } 
      else if (dataset$tool == "DIAUmpire") {
        #######################################################################
        ################# constructing dataset paths on s3 ###################
        s3_fragment_path <- gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$fragment_path))))
        s3_peptides_path <- gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$peptide_path))))
        s3_proteins_path <- gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$protein_path))))
        s3_annotation_path <- gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$annotation))))
        ################# Loading object to memory ###########################
        fragments <- get_file_from_s3(s3_file_path = s3_fragment_path,
                                      local_file_name = "fragments.RDS")
        peptides <- get_file_from_s3(s3_file_path = s3_peptides_path,
                                     local_file_name = "peptides.RDS")
        proteins <- get_file_from_s3(s3_file_path = s3_proteins_path,
                                     local_file_name = "proteins.RDS")
        annotation <- get_file_from_s3(s3_file_path = s3_annotation_path,
                                       local_file_name = "annotation.RDS")
        #######################################################################
        try({
          diau_conv = MSstatsdev::DIAUmpiretoMSstatsFormat(
            fragments, peptides, proteins, annotation,
            removeProtein_with1Feature = remove_single_feature,
            fewMeasurements = remove_few_lf)
          master_results <- run_dataprocess(data=diau_conv,
                                            dataset_path=s3_fragment_path,
                                            master_result_df=master_results,
                                            mb_impute=is_impute,
                                            censored_int="NA")
        })
      }else if (dataset$tool == "OpenMS") {
        #######################################################################
        ################# constructing dataset paths on s3 ###################
        s3_input_path = gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$input_path))))
        ################# Loading object to memory ###################
        input <- get_file_from_s3(s3_file_path = s3_input_path,
                                  local_file_name = "input.RDS")
        #######################################################################
        if (dataset$type == "TMT") {
          try({
            open_ms_converter = MSstatsTMTdev::OpenMStoMSstatsTMTFormat(
              input,rmPSM_withfewMea_withinRun = remove_few,
              rmProtein_with1Feature = remove_single_feature)})
        } else {
          try({
            open_ms_converter = MSstatsdev::OpenMStoMSstatsFormat(
              input, rmPSM_withfewMea_withinRun = remove_few,
              rmProtein_with1Feature = remove_single_feature)})
        }
        master_results <- run_dataprocess(data=open_ms_converter,
                                          dataset_path=s3_input_path,
                                          master_result_df=master_results,
                                          mb_impute=is_impute,
                                          censored_int="NA")
      } else {
        #######################################################################
        ################# constructing dataset paths on s3 ###################
        s3_input_data_path <- gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "", dataset$input_path))))
        s3_annotation_path = gsub(path_to_datasets, "", paste0(
          datasets$processed_data, gsub("tsv|csv|xls|txt", "RDS", sub(
            ".", "",dataset$annotation_path))))
        ################# Loading object to memory ###################
        input <- get_file_from_s3(s3_file_path = s3_input_data_path,
                                  local_file_name = "input.RDS")
        annotation <- get_file_from_s3(s3_file_path = s3_annotation_path,
                                       local_file_name = "annotation.RDS")
        input = as.data.frame(input)
        annotation = as.data.frame(annotation)
        #######################################################################
        
        if (dataset$tool == "PD") {
          if (dataset$type == "TMT") {
            pd_converter = MSstatsTMTdev::PDtoMSstatsTMTFormat(
              input, annotation, rmPSM_withfewMea_withinRun = remove_few,
              rmProtein_with1Feature = remove_single_feature)
          } else {
            pd_converter = MSstatsdev::PDtoMSstatsFormat(
              input, annotation,removeProtein_with1Peptide = remove_single_feature,
              fewMeasurements = remove_few_lf)
          }
          master_results <- run_dataprocess(
            data=pd_converter, dataset_path=s3_input_data_path,
            master_result_df=master_results,
            mb_impute=is_impute, censored_int="0")
        }
        if (dataset$tool == "OpenSWATH") {
          os_converter = MSstatsdev::OpenSWATHtoMSstatsFormat(
            input, annotation,removeProtein_with1Feature = remove_single_feature,
            fewMeasurements = remove_few_lf)
          master_results <- run_dataprocess(
            data=os_converter, dataset_path=s3_input_data_path,
            master_result_df=master_results,
            mb_impute=is_impute,censored_int="0")
        }
        if (dataset$tool == "Progenesis") {
          progenesis_converter = MSstatsdev::ProgenesistoMSstatsFormat(
            input, annotation,
            removeProtein_with1Peptide = remove_single_feature,
            fewMeasurements = remove_few_lf)
          master_results <- run_dataprocess(
            data=progenesis_converter, dataset_path=s3_input_data_path,
            master_result_df=master_results,
            mb_impute=is_impute,censored_int="0")
        }
        if (dataset$tool == "Skyline") {
          skyline_converter = MSstatsdev::SkylinetoMSstatsFormat(
            input, annotation,removeProtein_with1Feature = remove_single_feature,
            fewMeasurements = remove_few_lf)
          master_results <- run_dataprocess(
            data=skyline_converter, dataset_path=s3_input_data_path,
            master_result_df=master_results,
            mb_impute=is_impute,censored_int="0")
        }
        if (dataset$tool == "SpectroMine") {
          spectromine_converter = MSstatsTMTdev::SpectroMinetoMSstatsTMTFormat(
            input, annotation, rmPSM_withfewMea_withinRun = remove_few,
            rmProtein_with1Feature = remove_single_feature)
          master_results <- run_dataprocess(
            data=spectromine_converter, dataset_path=s3_input_data_path,
            master_result_df=master_results,
            mb_impute=is_impute,censored_int="0")
        }
        if (dataset$tool == "Spectronaut") {
          spectronaut_converter = MSstatsdev::SpectronauttoMSstatsFormat(
            input, annotation,removeProtein_with1Feature = remove_single_feature,
            fewMeasurements = remove_few_lf)
          master_results <- run_dataprocess(
            data=spectronaut_converter, dataset_path=s3_input_data_path,
            master_result_df=master_results,
            mb_impute=is_impute,censored_int="0")
        }
      }
    }, error = function(e) {
      problems <<- c(problems, dataset$folder_path)
      print(e)
    })
  }
  return(master_results)
}

problems = list()
result_df <-run_wider_testing(metadata_s3)

#uploading the results to s3 as csv
store_csv_file_to_s3(s3_path = results$code_deploy_results,
                     local_file_name = "report.xlsx", upload_file=result_df)
####### ##################### cleaning up #####################################
closeAllConnections()
rm(list=ls(all=TRUE))
###############################################################################
