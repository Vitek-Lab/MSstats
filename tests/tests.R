#!/usr/bin/env Rscript

# helper functions to read  binary files from s3
source("/home/rstudio/code/deployment/tests/utils/s3_helper_functions.R")

# `load()` R objects from S3
frags <- get_file_from_s3(s3_file_path = "Navarro2016_DIA_DIAumpire_input_FragSummary.RDS",
                     local_file_name = "frags.rds")

pepts <- s3$get_object(Bucket = aws_bucket_name,
                       Key = "Navarro2016_DIA_DIAumpire_input_PeptideSummary.RDS")
pepts <- read_bin_files_s3(pepts$Body, "pepts.rds")

prots <- s3$get_object(Bucket = aws_bucket_name,
                       Key = "Navarro2016_DIA_DIAumpire_input_ProtSummary.RDS")
prots <- read_bin_files_s3(prots$Body, "prots.rds")

dia_annot <- s3$get_object(Bucket = aws_bucket_name,
                           Key = "Navarro2016_DIA_DIAumpire_input_annotation.RDS")
dia_annot <- read_bin_files_s3(dia_annot$Body, "dia_annot.rds")

diau <- MSstatsdev::DIAUmpiretoMSstatsFormat(frags, pepts, prots, dia_annot)
# diau_v3 <- MSstats::dataProcess(as.data.frame(unclass(diau)), MBimpute = TRUE)
# diau_v4 <- MSstatsdev::dataProcess(diau, MBimpute = TRUE)

# to check if we can store the rds object to S3
s3$put_object(Body = write_bin_files_s3(diau, "dia_annot.rds"), 
              Bucket = aws_bucket_name, Key = "from-ec2.RDS")

rm(list=ls(all=TRUE))