library(MSstatsConvert)
library(MSstats)
library(stringr)
library(parallel)
library(jsonlite)

source("calculateMetrics.R")

config <- fromJSON("scriptController.json", simplifyVector = FALSE)

dataset_config <- config$datasets[["DDA-Dowell2021-HEqe408_LFQ"]]
dataset_config <- as.list(dataset_config)

cat("Processing Dataset:", dataset_config$name, "\n")
cat("Dataset File Path:", dataset_config$file, "\n")


start_time <- Sys.time()

fragpipe_raw <- data.table::fread(dataset_config$file)

head(fragpipe_raw)

msstats_format = MSstatsConvert::FragPipetoMSstatsFormat(fragpipe_raw, use_log_file = FALSE)


data_process_tasks <- list(
  list(
    label = "Data process with Normalized Data",
    result = function() dataProcess(msstats_format, featureSubset = "topN", n_top_feature = 20)
  ),
  list(
    label = "Data process with Normalization and MBImpute False",
    result = function() dataProcess(msstats_format, featureSubset = "topN", n_top_feature = 20, MBimpute = FALSE)
  ),
  list(
    label = "Data process without Normalization",
    result = function() dataProcess(msstats_format, featureSubset = "topN", normalization = "FALSE", n_top_feature = 20)
  ),
  list(
    label = "Data process without Normalization with MBImpute False",
    result = function() dataProcess(msstats_format, featureSubset = "topN", normalization = "FALSE", n_top_feature = 20, MBimpute = FALSE)
  )
)

start_time <- Sys.time()

num_cores <- detectCores() - 1 

summarized_results <- mclapply(data_process_tasks, function(task) {
  list(label = task$label, summarized = task$result())
}, mc.cores = num_cores)	


results_list <- mclapply(summarized_results, function(res) {
  calculateResult(res$summarized, res$label, dataset_config$samples)
}, mc.cores = num_cores)


final_results <- do.call(rbind, results_list)
end_time <- Sys.time()
total_time <- end_time - start_time
print(final_results)
print(paste("Total Execution Time:", total_time))