library(ggplot2)
library(dplyr)
library(stringr)

calculateResult <- function(summarized, label, samples) {
  model <- groupComparison("pairwise", summarized)
  comparisonResult <- model$ComparisonResult
  
  TP <- 0
  FP <- 0
  TN <- 0
  FN <- 0
  
  for (sample_name in names(samples)) {
    sample <- samples[[sample_name]]
    is_significant <- sample$type == "significant"
    
    filtered_proteins <- comparisonResult %>% filter(grepl(sample$pattern, Protein))
    
    if (is_significant) {
      TP <- TP + nrow(filtered_proteins %>% filter(adj.pvalue < 0.05))
      FN <- FN + nrow(filtered_proteins %>% filter(adj.pvalue >= 0.05))
    } else {
      FP <- FP + nrow(filtered_proteins %>% filter(adj.pvalue < 0.05))
      TN <- TN + nrow(filtered_proteins %>% filter(adj.pvalue >= 0.05))
    }
  }
  
  FPR <- FP / (FP + TN)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  recall <- TP / (TP + FN)
  fdr <- FP / (FP + TP)
  
  cat("Metrics for Label:", label, "\n")
  cat("True Positives (TP):", TP, "\n")
  cat("False Positives (FP):", FP, "\n")
  cat("True Negatives (TN):", TN, "\n")
  cat("False Negatives (FN):", FN, "\n\n")
  
  comparisonResult %>%
    filter(is.finite(log2FC)) %>%
    ggplot(aes(y = log2FC)) +
    geom_boxplot() +
    geom_hline(yintercept = -2, linetype = "dashed", color = "red", linewidth = 1.5) +
    theme_bw() +
    labs(title = paste("Boxplot of log2FC for", label), y = "log2FC")
  
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
