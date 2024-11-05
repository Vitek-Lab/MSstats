library(MSstatsConvert)
library(MSstats)
library(ggplot2)
library(dplyr)
library(stringr)
library(parallel)


calculateResult <- function(summarized, label){
  
  model = groupComparison("pairwise", summarized)
  comparisonResult <- model$ComparisonResult
  
  human_comparisonResult <- comparisonResult %>% filter(grepl("_HUMAN$", Protein))
  
  ecoli_comparisonResult <- comparisonResult %>% filter(grepl("_ECOLI$", Protein))
  
  yeast_comparisonResult <- comparisonResult %>% filter(grepl("_YEAST$", Protein))
  
  
  human_median <- median(human_comparisonResult$log2FC, na.rm = TRUE)
  ecoli_median <- median(ecoli_comparisonResult$log2FC, na.rm = TRUE) 
  yeast_median <- median(yeast_comparisonResult$log2FC, na.rm = TRUE)
  
  cat("Expected Log Change Human:", human_median, "\n")
  cat("Expected Log Change Ecoli:", ecoli_median, "\n")
  cat("Expected Log Change Yeast:", yeast_median, "\n")
  
  #calculate SD and mean
  
  
  # Kept the code for Individual Boxplots
  
  # boxplot(human_comparisonResult$log2FC,
  #         main = "Boxplot of log2FC for Human",
  #         ylab = "log2FC",
  #         col = "lightblue")
  # 
  # 
  boxplot(ecoli_comparisonResult$log2FC,
          main = "Boxplot of log2FC for E. coli",
          ylab = "log2FC",
          col = "lightgreen")
  # 
  # boxplot(yeast_comparisonResult$log2FC,
  #         main = "Boxplot of log2FC for Yeast",
  #         ylab = "log2FC",
  #         col = "lightpink")
  
  combined_data <- list(
    Human = human_comparisonResult$log2FC,
    Ecoli = ecoli_comparisonResult$log2FC,
    Yeast = yeast_comparisonResult$log2FC
  )
  
  
  unique_ecoli_proteins <- unique(ecoli_comparisonResult$Protein)
  unique_yeast_proteins <- unique(yeast_comparisonResult$Protein)
  
  all_proteins <- c(union(unique_ecoli_proteins, unique_yeast_proteins))  # find out the significant proteins in FragData
  
  extracted_proteins <- sapply(all_proteins, function(x) {
    split_string <- strsplit(x, "\\|")[[1]]  # Split the string by '|'
    if (length(split_string) >= 2) {
      return(split_string[2])  # Return the second element
    } else {
      return(NA)  # Return NA if there's no second element
    }
  })
  
  extracted_proteins <- unname(unlist(extracted_proteins))
  
  proteins <- c(extracted_proteins)
  
  
  TP <- comparisonResult %>% filter(grepl(paste(proteins, collapse = "|"), Protein) & adj.pvalue < 0.05) %>% nrow()
  
  
  FP <- comparisonResult %>% filter(!grepl(paste(proteins, collapse = "|"), Protein) & adj.pvalue < 0.05) %>% nrow()
  
  
  TN <- comparisonResult %>% filter(!grepl(paste(proteins, collapse = "|"), Protein) & adj.pvalue >= 0.05) %>% nrow()
  
  
  FN <- comparisonResult %>% filter(grepl(paste(proteins, collapse = "|"), Protein) & adj.pvalue >= 0.05) %>% nrow()
  
  cat("True Positives (Yeast and EColi): ", TP, "\n")
  cat("False Positives (Human Samples)", FP, "\n")
  cat("True Negatives", TN, "\n")
  cat("False Negatives", FN, "\n")
  
  FPR <- FP / (FP + TN)
  
  # Accuracy
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  # Recall
  recall <- TP / (TP + FN)
  
  results <- data.frame(
    Label = label,
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    FPR = FPR,
    Accuracy = accuracy,
    Recall = recall
  )
  
  return(results)
  
}

start_time <- Sys.time()

# Use fread directly to read the CSV
fragpipe_raw = data.table::fread("..//data//FragPipeMsStatsBenchmarking.csv")

head(fragpipe_raw)

fragpipe_raw$Condition = unlist(lapply(fragpipe_raw$Run, function(x){
  paste(str_split(x, "\\_")[[1]][4:5], collapse="_")
}))

fragpipe_raw$BioReplicate = unlist(lapply(fragpipe_raw$Run, function(x){
  paste(str_split(x, "\\_")[[1]][4:7], collapse="_")
}))

# Convert to MSstats format
msstats_format = MSstatsConvert::FragPipetoMSstatsFormat(fragpipe_raw, use_log_file = FALSE)


# Define the tasks with descriptive labels
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
    result = function() dataProcess(msstats_format, normalization = "FALSE", n_top_feature = 20)
  ),
  list(
    label = "Data process without Normalization with MBImpute False",
    result = function() dataProcess(msstats_format, normalization = "FALSE", n_top_feature = 20, MBimpute = FALSE)
  )
)

# Start the timer
start_time <- Sys.time()

# Use mclapply to run the dataProcess tasks in parallel
num_cores <- detectCores() - 1  # Use one less than the total cores available

# Run data processing tasks in parallel and collect results with labels
summarized_results <- mclapply(data_process_tasks, function(task) {
  list(label = task$label, summarized = task$result())
}, mc.cores = num_cores)

# Run calculateResult on each summarized result in parallel
results_list <- mclapply(summarized_results, function(res) {
  calculateResult(res$summarized, res$label)
}, mc.cores = num_cores)

# Combine all results into a single data frame
final_results <- do.call(rbind, results_list)

# End the timer
end_time <- Sys.time()
total_time <- end_time - start_time

# Display the final results and execution time
print(final_results)
print(paste("Total Execution Time:", total_time))