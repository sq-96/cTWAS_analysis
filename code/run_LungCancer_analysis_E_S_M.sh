#!/bin/bash

for WEIGHT in Lung
do 
	TRAIT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/ieu-b-4955.vcf.gz"
	LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1"
	WEIGHT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/Predictive_Models/Lung_E.db;/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/Predictive_Models/Lung_S.db;/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/Predictive_Models/Lung_M.db"
	CONFIG_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b38.R"
	OUTNAME_E="LungCancer_"$WEIGHT"_expr"
	OUTNAME="LungCancer_"$WEIGHT"_ctwas"
	OUTDIR="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/LungCancer_E_S_M/"$WEIGHT
	job_name="LungCancer_"$WEIGHT
	job_out="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/LungCancer_E_S_M_out/"$job_name".out"
	job_err="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/LungCancer_E_S_M_out/"$job_name".err"
  
	FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
	if [ ! -f $FINAL_FILE ]
	then 
		sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_LungCancer_analysis_E_S_M.sbatch
	fi
done

