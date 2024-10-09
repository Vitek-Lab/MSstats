library(MSstatsConvert)
library(MSstats)
library(ggplot2)
library(dplyr)
library(stringr)


# Use fread directly to read the CSV
fragpipe_raw = data.table::fread("..//data//FragPipeMsStatsBenchmarking.csv")

# Check the first few rows
head(fragpipe_raw)

fragpipe_raw$Condition = unlist(lapply(fragpipe_raw$Run, function(x){
    paste(str_split(x, "\\_")[[1]][4:5], collapse="_")
}))

fragpipe_raw$BioReplicate = unlist(lapply(fragpipe_raw$Run, function(x){
    paste(str_split(x, "\\_")[[1]][4:7], collapse="_")
}))

# Convert to MSstats format
msstats_format = MSstatsConvert::FragPipetoMSstatsFormat(fragpipe_raw, use_log_file = FALSE)


# summarized = dataProcess(msstats_format, featureSubset = "topN",
#                          n_top_feature = 20)

# MB_Impute
summarized = dataProcess(msstats_format, featureSubset = "topN",n_top_feature = 20, MBimpute = FALSE)

# use this 
# summarized = dataProcess(msstats_format, normalization = "FALSE")
# summarized = dataProcess(msstats_format, normalization = "FALSE", MBimpute = FALSE)

model = groupComparison("pairwise", summarized)

#TODO: with a contrast matrix


comparisonResult <- model$ComparisonResult

human_comparisonResult <- comparisonResult %>% filter(grepl("_HUMAN$", Protein))

ecoli_comparisonResult <- comparisonResult %>% filter(grepl("_ECOLI$", Protein))

yeast_comparisonResult <- comparisonResult %>% filter(grepl("_YEAST$", Protein))

# For each dataset, verify the expected log fold change for each set of proteins
# (e.g. in one paper, it claimed expected log fold changes to be 0, -1, and 2 for human, yeast, and ecoli respectively).

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

boxplot(combined_data,
        main = "Boxplots of expected Log fold for Human, E. coli, and Yeast",
        ylab = "log2FC",
        names = c("Human", "E. coli", "Yeast"),
        col = c("lightblue", "lightgreen", "lightpink"))


#Ecoli and Yeast are TP
#Humans is FP

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

cat("False Positive Rate: ", FPR, "\n")
cat("Accuracy: ", accuracy, "\n")
cat("Recall: ", recall, "\n")

significant_proteins <- comparisonResult %>% filter(adj.pvalue < 0.05)

upregulated_proteins <- significant_proteins %>% filter(log2FC > 0) %>% nrow()

downregulated_proteins <- significant_proteins %>% filter(log2FC < 0) %>% nrow()

difference_in_proteins <- upregulated_proteins - downregulated_proteins

cat("Upregulated proteins:", upregulated_proteins, "\n")
cat("Downregulated proteins:", downregulated_proteins, "\n")
cat("Difference (Up - Down):", difference_in_proteins, "\n")


# TODO should we do for a particular protein or for all?
# groupComparisonPlots(
#     model$ComparisonResult,
#     type="VolcanoPlot",
#     sig = 0.05,
#     FCcutoff = FALSE,
#     logBase.pvalue = 10,
#     ylimUp = FALSE,
#     ylimDown = FALSE,
#     xlimUp = FALSE,
#     x.axis.size = 10,
#     y.axis.size = 10,
#     dot.size = 3,
#     text.size = 4,
#     text.angle = 0,
#     legend.size = 13,
#     ProteinName = TRUE,
#     colorkey = TRUE,
#     numProtein = 100,
#     clustering = "both",
#     width = 800,
#     height = 600,
#     which.Protein = "all",
#     isPlotly = FALSE
# )