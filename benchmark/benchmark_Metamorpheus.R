library(MSstatsConvert)
library(MSstats)
library(parallel)
library(stringr)
library(jsonlite)
library(dplyr)

source("metamorpheus_Process.R")
config <- fromJSON("scriptController.json", simplifyVector = FALSE)

runBenchmarkForMetaMorpheusData <- function(datasetPath, config) {

  dataset_config <- config$datasets[[datasetPath]]
  dataset_config <- as.list(dataset_config)

  cat("Processing Dataset:", dataset_config$name, "\n")

  filePath <- file.path(dataset_config$parent, dataset_config$data)
  annotPath <- dataset_config$parent

  input = data.table::fread(file.path(filePath, "QuantifiedPeaks.tsv"))
  annot = data.table::fread(file.path(annotPath, "annotation.csv"))


  cat("Dataset File Path:", filePath, "\n")
  cat("Annotation File Path:", annotPath, "\n")

  input = input %>% filter(!str_detect(`Protein Group`, ";")) # remove multiple protein group in same cell
  input = input %>% filter(!str_detect(`Protein Group`, "DECOY")) # remove decoys

  protein_mappings = data.table::fread(file.path(filePath, "QuantifiedProteins.tsv"))
  valid_organisms <- unique(input$Organism)

  protein_mappings = protein_mappings %>% filter(Organism %in% valid_organisms)

  input = input %>% filter(`Protein Group` %in% protein_mappings$`Protein Groups`)

  output = MetamorpheusToMSstatsFormat(input, annot)

  data_process_tasks <- list(
    list(
      label = "Data process with Normalized Data",
      result = function() dataProcess(output, featureSubset = "topN", n_top_feature = 20)
    ),
    list(
      label = "Data process with Normalization and MBImpute False",
      result = function() dataProcess(output, featureSubset = "topN", n_top_feature = 20, MBimpute = FALSE)
    ),
    list(
      label = "Data process without Normalization",
      result = function() dataProcess(output, featureSubset = "topN", normalization = "FALSE", n_top_feature = 20)
    ),
    list(
      label = "Data process without Normalization with MBImpute False",
      result = function() dataProcess(output, featureSubset = "topN", normalization = "FALSE", n_top_feature = 20, MBimpute = FALSE)
    ),
    list(
      label = "Data process without Normalization and Imputation On for all features",
      result = function() dataProcess(output, featureSubset = "all", normalization = "FALSE", MBimpute = FALSE)
    ),
    list(
      label = "Data process without Normalization and Imputation On for top3 features",
      result = function() dataProcess(output, featureSubset = "top3", normalization = "FALSE", MBimpute = FALSE)
    )
  )
      
  start_time <- Sys.time()

  num_cores <- detectCores() - 1 

  summarized_results <- mclapply(data_process_tasks, function(task) {
    list(label = task$label, summarized = task$result())
  }, mc.cores = num_cores)	


  results_list <- mclapply(summarized_results, function(res) {
    calculate_Metrics(res$summarized, protein_mappings, res$label)
  }, mc.cores = num_cores)


  final_results <- do.call(rbind, results_list)
  end_time <- Sys.time()
  total_time <- end_time - start_time
  print(final_results)
  print(paste("Total Execution Time:", total_time))

}



runBenchmarkForMetaMorpheusData("DDA-Solivais2024-Metamorpheus_MBR_LFQ", config)
runBenchmarkForMetaMorpheusData("DDA-Solivais2024-Metamorpheus_NoMBR_LFQ", config)