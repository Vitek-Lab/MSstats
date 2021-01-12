#!/usr/bin/env Rscript

#############################################################
# load all the config settings/credentials
source("/home/rstudio/code/config.R")
# set aws session variables
Sys.setenv(
  AWS_ACCESS_KEY_ID = aws_key,
  AWS_SECRET_ACCESS_KEY = aws_secret,
  AWS_REGION = aws_region,
  AWS_S3_BUCKET = aws_bucket_name
)
############################################################


# initialise s3 object
s3 <- paws::s3()
# quick smoke test to check bucket contents
print(s3$list_buckets())
print(getwd())

read_bin_files_s3 <- function(bin_file, file_name){
  # helper functions to read files from s3(paws package streams file as binary)
  writeBin(bin_file, con = file_name)
  rds_file <- readRDS(file_name)
  unlink(file_name)
  return (rds_file)
}

write_bin_files_s3 <- function(r_object, file_name, is_rds=TRUE){
  # helper functions to write binary files to s3(paws package streams bin file)
  if(is_rds){
    saveRDS(r_object, file = file_name)
  }
  # Load the file as a raw binary
  print(file_name)
  read_file <- file(file_name, "rb")
  bin_file <- readBin(read_file, "raw", n = file.size(file_name))
  unlink(file_name)
  return(bin_file)
}

get_file_from_s3 <- function(s3_file_path, local_file_name){
  # fetches file from aws s3
  s3_stream_file <- s3$get_object(Bucket = aws_bucket_name, Key=s3_file_path)
  return(read_bin_files_s3(s3_stream_file$Body, local_file_name))
}

store_csv_file_to_s3 <- function(s3_path, local_file_name, upload_file){
  # helper function to upload csv to s3
  print("uploading results to s3...")
  s3_file_path <- generate_s3_path(s3_path)
  # generate_xlsx(upload_file, local_file_name)
  write.csv(upload_file, file=local_file_name)
  s3$put_object(
    Body = write_bin_files_s3(local_file_name, local_file_name, is_rds = FALSE), 
    Bucket = aws_bucket_name, 
    Key = s3_file_path)
  print("upload to s3 finished. Results located in: ")
  print(s3_file_path)
  closeAllConnections()
  unlink(local_file_name)
}

generate_s3_path <- function(s3_path){
  # helper function that appends datetime to the filename(to create unique names)
  s3_file_path <- paste(s3_path, "report-", 
                        paste(
                          as.Date(format(Sys.time(), "%D"), format = "%m/%d/%y"), 
                          format(Sys.time(), "-%H-%M-%S"), ".xlsx",
                          sep = ""), 
                        sep = "")
  return(s3_file_path)
}

store_rds <- function(file_up, local_name, some_path){
  s3$put_object(
    Body = write_bin_files_s3(file_up, local_name, is_rds = T), 
    Bucket = aws_bucket_name, 
    Key = some_path)
}