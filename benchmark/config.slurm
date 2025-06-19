#!/bin/bash
#SBATCH --job-name=msstats_benchmark_job_updated
#SBATCH --chdir=/projects/VitekLab/Projects/Benchmarking/benchmark
#SBATCH --output=job_output.txt
#SBATCH --error=job_error.txt
#SBATCH --open-mode=append
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=short

module load R
module load cmake/3.30.2

export LC_ALL=C
export R_LIBS_USER=/home/raina.ans/R/x86_64-pc-linux-gnu-library/4.4
mkdir -p $R_LIBS_USER

mkdir -p $HOME/lib_fix
ln -sf /shared/EL9/explorer/R/4.4.1/lib64/R/lib/libRlapack.so $HOME/lib_fix/libRlapack.so.3

export LD_LIBRARY_PATH=$HOME/lib_fix:/shared/EL9/explorer/R/4.4.1/lib64/R/lib:/usr/lib64:$LD_LIBRARY_PATH

echo "Running ldd on preprocessCore.so" >> job_output.txt
ldd $(find ~/R -name preprocessCore.so | head -n 1) >> job_output.txt 2>> job_error.txt
echo "Running ldd on nloptr.so" >> job_output.txt
ldd $(find ~/R -name nloptr.so | head -n 1) >> job_output.txt 2>> job_error.txt

Rscript -e '
.libPaths("/home/raina.ans/R/x86_64-pc-linux-gnu-library/4.4")
Sys.unsetenv("R_LIBS")
cat(".libPaths():\n"); print(.libPaths())

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")

install.packages("nloptr", type = "source", repos = "https://cloud.r-project.org")

for (pkg in c("dplyr", "stringr", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}
BiocManager::install(c("MSstatsConvert", "preprocessCore"), force = TRUE)
remotes::install_github("Vitek-Lab/MSstats", ref = "devel", force = TRUE)
'

R_SCRIPTS=("benchmark_Dowell2021-HEqe408_LFQ.R" "benchmark_Puyvelde2022-HYE5600735_LFQ.R" "benchmark_Metamorpheus.R")

for script in "${R_SCRIPTS[@]}"; do
    echo "Executing script: $script" >> job_output.txt
    stdbuf -oL -eL Rscript "$script" >> job_output.txt 2>> job_error.txt
    wait
    echo "Finished executing script: $script" >> job_output.txt
    echo -e "\n\n"
done