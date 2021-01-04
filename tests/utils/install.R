#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

if (!requireNamespace("paws", quietly = TRUE))
    install.packages("paws")

if (!requireNamespace("xlsx", quietly = TRUE))
    install.packages("xlsx")

BiocManager::install("MSstats")
devtools::install_github("Vitek-Lab/MSstats-dev", ref = "develop")
devtools::install_github("Vitek-Lab/MSstatsTMT", ref = "develop")
