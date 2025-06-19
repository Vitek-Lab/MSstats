library(MSstatsConvert)
library(MSstats)
library(stringr)
library(parallel)
library(jsonlite)

source("metamorpheus_Process.R")

options(echo = TRUE)

config <- fromJSON("scriptController.json", simplifyVector = FALSE)
dataset_config <- config$datasets[["DDA-Solivais2024-Metamorpheus_NoMBR_LFQ"]]
dataset_config <- as.list(dataset_config)

cat("Processing Dataset:", dataset_config$name, "\n")
filePath <- file.path(dataset_config$parent, dataset_config$data)
annotPath <- dataset_config$parent
cat("File Path:", filePath, "\n")
cat("Annotation Path:", annotPath, "\n")
flush.console()

input <- data.table::fread(file.path(filePath, "QuantifiedPeaks.tsv"))
cat("Input loaded: ", nrow(input), " rows\n")
flush.console()

annot <- data.table::fread(file.path(annotPath, "annotation.csv"))
cat("Annotation loaded: ", nrow(annot), " rows\n")
flush.console()

# Filters
input <- input %>% filter(!str_detect(`Protein Group`, ";"))
cat("Post semicolon filter: ", nrow(input), " rows\n")
input <- input %>% filter(!str_detect(`Protein Group`, "DECOY"))
cat("Post DECOY filter: ", nrow(input), " rows\n")
flush.console()

protein_mappings <- data.table::fread(file.path(filePath, "QuantifiedProteins.tsv"))
cat("Protein mappings loaded: ", nrow(protein_mappings), " rows\n")
protein_mappings <- protein_mappings %>% filter(Organism %in% c("Escherichia coli (strain K12)", "Homo sapiens"))
cat("Post organism filter: ", nrow(protein_mappings), " rows\n")
flush.console()

input <- input %>% filter(`Protein Group` %in% protein_mappings$`Protein Groups`)
cat("Post protein group mapping filter: ", nrow(input), " rows\n")
flush.console()

output <- tryCatch({
  MetamorpheusToMSstatsFormat(input, annot)
}, error = function(e) {
  message("Error in MetamorpheusToMSstatsFormat: ", e$message)
  NULL
})
cat("MSstats format complete\n")
flush.console()

if (is.null(output)) stop("Output is NULL, aborting")

# DEBUG: check structure
cat("Structure of MSstats formatted output:\n")
print(str(output))
flush.console()

data_process_tasks <- list(
  list(label = "Normalized", result = function() dataProcess(output, featureSubset = "topN", n_top_feature = 20)),
  list(label = "Norm + MBimpute=FALSE", result = function() dataProcess(output, featureSubset = "topN", n_top_feature = 20, MBimpute = FALSE)),
  list(label = "No Norm", result = function() dataProcess(output, featureSubset = "topN", normalization = "FALSE", n_top_feature = 20)),
  list(label = "No Norm + MBimpute=FALSE", result = function() dataProcess(output, featureSubset = "topN", normalization = "FALSE", n_top_feature = 20, MBimpute = FALSE)),
  list(label = "No Norm + Impute all", result = function() dataProcess(output, featureSubset = "all", normalization = "FALSE", MBimpute = FALSE)),
  list(label = "No Norm + Impute top3", result = function() dataProcess(output, featureSubset = "top3", normalization = "FALSE", MBimpute = FALSE))
)

start_time <- Sys.time()
num_cores <- detectCores() - 1 

# Debug: Use lapply instead of mclapply
summarized_results <- lapply(data_process_tasks, function(task) {
  cat("Running task:", task$label, "\n")
  flush.console()
  list(label = task$label, summarized = task$result())
})

cat("Completed summarization tasks\n")
flush.console()

results_list <- lapply(summarized_results, function(res) {
  cat("Processing result for:", res$label, "\n")
  flush.console()
  tryCatch({
    out <- calculate_Metrics(res$summarized, protein_mappings, res$label)
    print(str(out))
    flush.console()
    out
  }, error = function(e) {
    message("Error in calculate_Metrics for ", res$label, ": ", e$message)
    NULL
  })
})

cat("All metrics calculated\n")
flush.console()

final_results <- tryCatch({
  do.call(rbind, results_list)
}, error = function(e) {
  message("Error during rbind: ", e$message)
  NULL
})

end_time <- Sys.time()
total_time <- end_time - start_time

cat("Final Results:\n")
print(final_results)
print(paste("Total Execution Time:", total_time))
flush.console()
