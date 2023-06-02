#!/bin/bash

for WEIGHT in WholeBlood
do 
	TRAIT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/ukb-d-30000_irnt.vcf.gz"
	LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1"
	WEIGHT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/predictive_models_protein_coding/WholeBlood_E.db;/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/predictive_models_protein_coding/WholeBlood_S.db;/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/predictive_models_protein_coding/WholeBlood_M.db"
	CONFIG_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b38.R"
	OUTNAME_E="WhiteBlood_"$WEIGHT"_expr"
	OUTNAME="WhiteBlood_"$WEIGHT"_ctwas"
	OUTDIR="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/WhiteBlood_E_S_M_PC/"$WEIGHT
	job_name="WhiteBlood_"$WEIGHT
	job_out="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/WhiteBlood_E_S_M_PC_out/"$job_name".out"
	job_err="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/WhiteBlood_E_S_M_PC_out/"$job_name".err"
  
	FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
	if [ ! -f $FINAL_FILE ]
	then 
		sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_WhiteBlood_analysis_E_S_M.sbatch
	fi
done

