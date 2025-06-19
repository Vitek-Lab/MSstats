library(MSstatsConvert)
library(MSstats)
library(stringr)
library(parallel)
library(jsonlite)

source("metamorpheus_Process.R")

config <- fromJSON("scriptController.json", simplifyVector = FALSE)

dataset_config <- config$datasets[["DDA-Solivais2024-Metamorpheus_NoMBR_LFQ"]]
dataset_config <- as.list(dataset_config)

cat("Processing Dataset:", dataset_config$name, "\n")

filePath <- file.path(dataset_config$parent, dataset_config$data)
annotPath <- dataset_config$parent

input = data.table::fread(file.path(filePath, "QuantifiedPeaks.tsv"))
annot = data.table::fread(file.path(annotPath, "annotation.csv"))


cat("Dataset File Path:", filePath, "\n")
cat("Dataset File Path:", annotPath, "\n")

input = input %>% filter(!str_detect(`Protein Group`, ";")) # remove multiple protein group in same cell
input = input %>% filter(!str_detect(`Protein Group`, "DECOY")) # remove decoys

protein_mappings = data.table::fread(file.path(filePath, "QuantifiedProteins.tsv"))
protein_mappings = protein_mappings %>% filter(Organism %in% c("Escherichia coli (strain K12)", "Homo sapiens"))

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

cat("Summarized Results:\n")
print(str(summarized_results))
flush.console()


results_list <- lapply(summarized_results, function(res) {
  cat("Processing result for:", res$label, "\n")
  flush.console()
  out <- tryCatch({
    calculate_Metrics(res$summarized, protein_mappings, res$label)
  }, error = function(e) {
    message("Error in calculate_Metrics for ", res$label, ": ", e$message)
    NULL
  })
  print(str(out))
  flush.console()
  out
})

cat("Results List structure:\n")
print(str(results_list))
flush.console()

final_results <- tryCatch({
  do.call(rbind, results_list)
}, error = function(e) {
  message("Error during rbind: ", e$message)
  NULL
})


end_time <- Sys.time()
total_time <- end_time - start_time
print(final_results)
print(paste("Total Execution Time:", total_time))