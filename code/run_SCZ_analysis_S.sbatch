#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-xinhe
#SBATCH --cpus-per-task=14
#SBATCH --mem=56G
#SBATCH --time=36:00:00

module load R/3.6.1

Rscript run_SCZ_ctwas_rss_LDR_S.R $TRAIT_FILE $LD_DIR $WEIGHT_FILE $CONFIG_FILE $OUTNAME_E $OUTNAME $OUTDIR

