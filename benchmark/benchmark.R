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
  
  boxplot(ecoli_comparisonResult$log2FC,
          main = "Boxplot of log2FC for E. coli",
          ylab = "log2FC",
          col = "lightgreen")
  combined_data <- list(
    Human = human_comparisonResult$log2FC,
    Ecoli = ecoli_comparisonResult$log2FC,
    Yeast = yeast_comparisonResult$log2FC
  )
  
  
  unique_ecoli_proteins <- unique(ecoli_comparisonResult$Protein)
  unique_yeast_proteins <- unique(yeast_comparisonResult$Protein)
  
  all_proteins <- c(union(unique_ecoli_proteins, unique_yeast_proteins)) 
  
  extracted_proteins <- sapply(all_proteins, function(x) {
    split_string <- strsplit(x, "\\|")[[1]]  
    if (length(split_string) >= 2) {
      return(split_string[2])  
    } else {
      return(NA)  
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
  
  ecoli_comparisonResult %>% 
      filter(is.finite(log2FC)) %>% 
      ggplot() + geom_boxplot(aes(y=log2FC)) +
      geom_hline(aes(yintercept=-2), linetype="dashed", color="red", linewidth=2) +
      theme_bw() + 
      labs(title = "Boxplot of log2FC for E. coli",
           y = "log2FC")
  
  human_comparisonResult %>% 
      filter(is.finite(log2FC)) %>% 
      ggplot() + geom_boxplot(aes(y=log2FC)) +
      geom_hline(aes(yintercept=-2), linetype="dashed", color="blue", linewidth=2) +
      theme_bw() + 
      labs(title = "Boxplot of log2FC for Human",
           y = "log2FC")
  
  
  yeast_comparisonResult %>% 
      filter(is.finite(log2FC)) %>% 
      ggplot() + geom_boxplot(aes(y=log2FC)) +
      geom_hline(aes(yintercept=-2), linetype="dashed", color="green", linewidth=2) +
      theme_bw() + 
      labs(title = "Boxplot of log2FC for Yeast",
           y = "log2FC")
  
  FPR <- FP / (FP + TN)
  
  # Accuracy
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  # Recall
  recall <- TP / (TP + FN)

  # False Discover Rate (FDR)
  
  fdr <- FP/ (FP + TP)

  results <- data.frame(
    Label = label,
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    FPR = FPR,
    Accuracy = accuracy,
    Recall = recall,
    FDR = fdr
  )
  
  return(results)
  
}

start_time <- Sys.time()

fragpipe_raw = data.table::fread("/work/VitekLab/Data/MS/Benchmarking/DDA-Puyvelde2022/DDA-Puyvelde2022-HYE5600735_LFQ/FragPipe/TOP0/MSstats.csv")

head(fragpipe_raw)

fragpipe_raw$Condition = unlist(lapply(fragpipe_raw$Run, function(x){
  paste(str_split(x, "\\_")[[1]][4:5], collapse="_")
}))

fragpipe_raw$BioReplicate = unlist(lapply(fragpipe_raw$Run, function(x){
  paste(str_split(x, "\\_")[[1]][4:7], collapse="_")
}))

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
  calculateResult(res$summarized, res$label)
}, mc.cores = num_cores)

final_results <- do.call(rbind, results_list)
end_time <- Sys.time()
total_time <- end_time - start_time
print(final_results)
print(paste("Total Execution Time:", total_time))