#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

if (!requireNamespace("paws", quietly = TRUE))
    install.packages("paws")

if (!requireNamespace("xlsx", quietly = TRUE))
    install.packages("xlsx")

# update.packages(oldPkgs = old.packages(c("devtools", "paws", "xlsx", "BiocManager")),
#                 ask = FALSE,
#                 lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.0")

# BiocManager::install("MSstats", ask = FALSE,
#                      lib = "~/R/x86_64-pc-linux-gnu-library/4.0")
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")
devtools::install_github("mstaniak/MSstats", ref = "patch-2", 
                         dependencies = NA, quiet = TRUE,
                         upgrade = "never")
devtools::install_github("Vitek-Lab/MSstats-dev", ref = "develop", 
                         dependencies = NA, quiet = TRUE,
                         upgrade = "never")
devtools::install_github("Vitek-Lab/MSstatsTMT", ref = "develop", 
                         dependencies = NA, quiet = TRUE,
                         upgrade = "never")
