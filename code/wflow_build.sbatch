#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --partition=xinhe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=wflow_build
#SBATCH --output=wflow_build-%J.out
#SBATCH --error=wflow_build-%J.err

module load R/3.6.1

Rscript wflow_build.R