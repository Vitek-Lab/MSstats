library(MSstatsConvert)
library(MSstats)
library(tidyverse)


# No MBR
input_no_mbr = data.table::fread("/projects/VitekLab/Data/MS/Benchmarking/DDA-Solivais2024_Metamorpheus/Current/FlashLFQ_NoNormalization_NoPIP/QuantifiedPeaks.tsv")



annot_no_mbr = data.table::fread("/projects/VitekLab/Data/MS/Benchmarking/DDA-Solivais2024_Metamorpheus/Current/FlashLFQ_NoNormalization_NoPIP/annotation.csv")


input_no_mbr = input_no_mbr %>% filter(!str_detect(`Protein Group`, ";")) # remove multiple protein group in same cell
input_no_mbr = input_no_mbr %>% filter(!str_detect(`Protein Group`, "DECOY")) # remove decoys
protein_mappings = data.table::fread("/projects/VitekLab/Data/MS/Benchmarking/DDA-Solivais2024_Metamorpheus/Current/FlashLFQ_NoNormalization_NoPIP/QuantifiedProteins.tsv")
protein_mappings = protein_mappings %>% filter(Organism %in% c("Escherichia coli (strain K12)", "Homo sapiens"))
input_no_mbr = input_no_mbr %>% filter(`Protein Group` %in% protein_mappings$`Protein Groups`)


# input_no_mbr$`Protein Group` = ifelse(
#   input_no_mbr$`Protein Group` %in% ecoli$`Protein Groups`, 
#   paste(input_no_mbr$`Protein Group`, "|ECOLI", sep = ""), 
#   paste(input_no_mbr$`Protein Group`, "|HUMAN", sep = ""))
# write.csv(input_no_mbr, "QuantifiedPeaks.csv", row.names = FALSE)


output_no_mbr = MetamorpheusToMSstatsFormat(input_no_mbr, annot_no_mbr)
QuantData_no_mbr = dataProcess(output_no_mbr, normalization = FALSE)

dataProcessPlots(QuantData_no_mbr, "QCPlot", which.Protein = "allonly", address = FALSE, isPlotly = TRUE)

comparison <- matrix(c(-1,0,0,0,1, # 3x
                       -1,0,0,1,0, # 2.5x
                       -1,0,1,0,0, # 2x
                       -1,1,0,0,0),nrow=4,byrow = TRUE) # 1.5x


row.names(comparison) <- c("E-A", "D-A", "C-A", "B-A")
groups = levels(QuantData_no_mbr$ProteinLevelData$GROUP)
colnames(comparison) <- groups[order(as.numeric(groups))]
model_no_mbr <- groupComparison(contrast.matrix=comparison, data=QuantData_no_mbr,
                                use_log_file = FALSE)

library(tidyverse)
ecoli_proteins = protein_mappings %>% filter(Organism == "Escherichia coli (strain K12)")
model_no_mbr$ComparisonResult = model_no_mbr$ComparisonResult %>% mutate(ecoli = Protein %in% ecoli_proteins$`Protein Groups`)

e_group_no_mbr = model_no_mbr$ComparisonResult %>% filter(Label == "B-A") %>% filter(is.na(issue))
ecoli_no_mbr = e_group_no_mbr %>% filter(ecoli == TRUE)
hist(ecoli_no_mbr$log2FC)

ecoli_no_mbr = e_group_no_mbr %>% filter(adj.pvalue < 0.05) %>% filter(ecoli == TRUE)
human_no_mbr = e_group_no_mbr %>% filter(adj.pvalue < 0.05) %>% filter(ecoli == FALSE)
FDR_no_mbr = nrow(human_no_mbr) / (nrow(ecoli_no_mbr) + nrow(human_no_mbr))

cat("FDR no MBR", FDR_no_mbr, "\n")

# With MBR
library(MSstatsConvert)
library(MSstats)
library(tidyverse)
input = data.table::fread("/projects/VitekLab/Data/MS/Benchmarking/DDA-Solivais2024_Metamorpheus/Current/FlashLFQ_v1.0_NoNormalization_wPIP/QuantifiedPeaks.tsv")
annot = data.table::fread("/projects/VitekLab/Data/MS/Benchmarking/DDA-Solivais2024_Metamorpheus/Current/FlashLFQ_NoNormalization_NoPIP/annotation.csv")

input = input %>% filter(!str_detect(`Protein Group`, ";")) # remove multiple protein group in same cell
input = input %>% filter(!str_detect(`Protein Group`, "DECOY")) # remove decoys

protein_mappings = data.table::fread("/projects/VitekLab/Data/MS/Benchmarking/DDA-Solivais2024_Metamorpheus/Current/FlashLFQ_v1.0_NoNormalization_wPIP/QuantifiedProteins.tsv")
protein_mappings = protein_mappings %>% filter(Organism %in% c("Escherichia coli (strain K12)", "Homo sapiens"))
input = input %>% filter(`Protein Group` %in% protein_mappings$`Protein Groups`)

# input$`Protein Group` = ifelse(
#   input$`Protein Group` %in% ecoli$`Protein Groups`,
#   paste(input$`Protein Group`, "|ECOLI", sep = ""),
#   paste(input$`Protein Group`, "|HUMAN", sep = ""))
# write.csv(input, "QuantifiedPeaks-MBR.csv", row.names = FALSE)


output = MetamorpheusToMSstatsFormat(input, annot)
QuantData = dataProcess(output, normalization = FALSE)

dataProcessPlots(QuantData, "QCPlot", which.Protein = "allonly", address = FALSE, isPlotly = TRUE)

comparison <- matrix(c(-1,0,0,0,1, # 3x
                       -1,0,0,1,0, # 2.5x
                       -1,0,1,0,0, # 2x
                       -1,1,0,0,0),nrow=4,byrow = TRUE) # 1.5x
row.names(comparison) <- c("E-A", "D-A", "C-A", "B-A")
groups = levels(QuantData$ProteinLevelData$GROUP)
colnames(comparison) <- groups[order(as.numeric(groups))]
model <- groupComparison(contrast.matrix=comparison, data=QuantData,
                         use_log_file = FALSE)

library(tidyverse)
ecoli_proteins = protein_mappings %>% filter(Organism == "Escherichia coli (strain K12)")
model$ComparisonResult = model$ComparisonResult %>% mutate(ecoli = Protein %in% ecoli_proteins$`Protein Groups`)

e_group = model$ComparisonResult %>% filter(Label == "B-A") %>% filter(is.na(issue))
ecoli = e_group %>% filter(ecoli == TRUE)
hist(ecoli$log2FC)

ecoli = e_group %>% filter(adj.pvalue < 0.05) %>% filter(ecoli == TRUE)
human = e_group %>% filter(adj.pvalue < 0.05) %>% filter(ecoli == FALSE)
FDR = nrow(human) / (nrow(ecoli) + nrow(human))


cat("FDR MBR", FDR, "\n")
# FDR no MBR seems to be lower than that of FDR with MBR (except for E-A label), but it's not by a wide margin.  
# When normalization was enabled, FDR spiked to 38% without MBR and 58% with MBR.
# When we set adj.pvalue to 0.01, FDR without MBR does better, but not by much.
# Less proteins detected as significant with MBR disabled.