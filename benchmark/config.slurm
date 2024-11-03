#!/bin/bash
#SBATCH --job-name=msstats_benchmark_job
#SBATCH --output=job_output.txt
#SBATCH --error=job_error.txt
#SBATCH --time=01:00:00         # Set the maximum run time
#SBATCH --ntasks=1              # Number of tasks (one process)
#SBATCH --cpus-per-task=8       # Use 8 CPU cores for the task
#SBATCH --mem=256G              # Request 256GB of memory
#SBATCH --partition=short       # Use the 'short' partition (or change as needed)

export R_LIBS_USER=/home/raina.ans/library 

mkdir -p $R_LIBS_USER

module load R
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org'); \
BiocManager::install('MSstats', lib = Sys.getenv('R_LIBS_USER'), update = FALSE); \
BiocManager::install('MSstatsConvert', lib = Sys.getenv('R_LIBS_USER'), update = FALSE); \
install.packages(c('dplyr', 'stringr', 'ggplot2'), lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org')"
Rscript benchmark.R