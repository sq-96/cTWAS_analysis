#!/bin/bash

TRAIT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/scz_2018.RDS"
LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1_b37"
WEIGHT_FILE="/project2/xinhe/shengqian/fusion_twas-master/methylation_weights_3kb"
CONFIG_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b37.R"
OUTNAME_E="SCZ_ME_expr"
OUTNAME="SCZ_ME_ctwas"
OUTDIR="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/SCZ_2018_ME"
job_name="SCZ_2018_ME"
job_out="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/SCZ_2018_ME_out/"$job_name".out"
job_err="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/SCZ_2018_ME_out/"$job_name".err"
  
FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
if [ ! -f $FINAL_FILE ]
then 
	sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_SCZ_2018_analysis_ME.sbatch
fi

