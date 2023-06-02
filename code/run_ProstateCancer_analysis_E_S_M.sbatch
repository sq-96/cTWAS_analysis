#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-xinhe
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=36:00:00

module load R/4.1.0

Rscript run_ProstateCancer_ctwas_rss_LDR_E_S_M.R $TRAIT_FILE $LD_DIR $WEIGHT_FILE $CONFIG_FILE $OUTNAME_E $OUTNAME $OUTDIR

