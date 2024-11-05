#!/bin/bash
#SBATCH --job-name=msstats_benchmark_job
#SBATCH --output=job_output.txt
#SBATCH --error=job_error.txt
#SBATCH --time=01:00:00         # Set the maximum run time
#SBATCH --ntasks=1              # Number of tasks (one process)
#SBATCH --cpus-per-task=8      # Use 8 CPU cores for the task
#SBATCH --mem=128G              # Request 256GB of memory
#SBATCH --partition=short       # Use the 'short' partition (or change as needed)

module load R-geospatial


module load gcc/11.1.0
module load cmake/3.23.2

export LC_ALL=C
export R_LIBS_USER=/home/raina.ans/R/x86_64-pc-linux-gnu-library/4.2-geospatial


mkdir -p $R_LIBS_USER

module load R
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org'); \
BiocManager::install('MSstats', lib = Sys.getenv('R_LIBS_USER'), update = FALSE); \
BiocManager::install('MSstatsConvert', lib = Sys.getenv('R_LIBS_USER'), update = FALSE); \
install.packages(c('dplyr', 'stringr', 'ggplot2'), lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org')"

Rscript benchmark.R