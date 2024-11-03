#!/bin/bash
#SBATCH --job-name=msstats_benchmark_job
#SBATCH --output=job_output.txt
#SBATCH --error=job_error.txt
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --partition=short

# Set custom library path for R packages
export R_LIBS_USER=/home/raina.ans/library

# Load updated GCC and CMake if available
module load gcc/9.3.0  # Replace with an available version if needed
module load cmake/3.18.4  # Replace as needed

# Create the library directory if it doesn’t exist
mkdir -p $R_LIBS_USER

# Load R module
module load R

# Install packages
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org'); \
BiocManager::install('MSstats', lib = Sys.getenv('R_LIBS_USER'), update = FALSE); \
BiocManager::install('MSstatsConvert', lib = Sys.getenv('R_LIBS_USER'), update = FALSE); \
install.packages(c('dplyr', 'stringr', 'ggplot2'), lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org')"

# Run the R script
Rscript benchmark.R
