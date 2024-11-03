#!/bin/bash
#SBATCH --job-name=msstats_benchmark_job
#SBATCH --output=job_output.txt
#SBATCH --error=job_error.txt
#SBATCH --time=01:00:00         # Set the maximum run time
#SBATCH --ntasks=1              # Number of tasks (one process)
#SBATCH --cpus-per-task=8       # Use 8 CPU cores for the task
#SBATCH --mem=256G              # Request 256GB of memory
#SBATCH --partition=short       # Use the 'short' partition (or change as needed)

module load R
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
BiocManager::install('MSstats', update = FALSE); \
BiocManager::install('MSstatsConvert', update = FALSE); \
install.packages('ggplot2', repos = 'https://cloud.r-project.org'); \
install.packages('stringr', repos = 'https://cloud.r-project.org'); \
install.packages('dplyr', repos = 'https://cloud.r-project.org')"
Rscript benchmark.R