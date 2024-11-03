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
Rscript R/benchmark.R