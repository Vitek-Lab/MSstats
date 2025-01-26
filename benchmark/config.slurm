#!/bin/bash
#SBATCH --job-name=msstats_benchmark_job
#SBATCH --chdir=/work/VitekLab/Projects/Benchmarking/
#SBATCH --output=job_output.txt
#SBATCH --error=job_error.txt
#SBATCH --time=01:00:00         # Set the maximum run time
#SBATCH --ntasks=1              # Number of tasks (one process)
#SBATCH --cpus-per-task=8      # Use 8 CPU cores for the task
#SBATCH --mem=256G              # Request 256GB of memory
#SBATCH --partition=short       # Use the 'short' partition (or change as needed)

module load R-geospatial

module load gcc/11.1.0
module load cmake/3.23.2

export LC_ALL=C
export R_LIBS_USER=/home/raina.ans/R/x86_64-pc-linux-gnu-library/4.2-geospatial


mkdir -p $R_LIBS_USER

module load R
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org'); \
remotes::install_github('Vitek-Lab/MSstats', ref = 'devel', lib = Sys.getenv('R_LIBS_USER')); \
remotes::install_github('Vitek-Lab/MSstatsConvert', ref = 'master', lib = Sys.getenv('R_LIBS_USER')); \
install.packages(c('dplyr', 'stringr', 'ggplot2'), lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org')"

R_SCRIPTS=("benchmark.R" "benchmark1.R")

for script in $R_SCRIPTS; do
    echo "Executing script: $script"
    Rscript "$script"
done