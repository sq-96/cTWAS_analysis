#!/bin/bash

for WEIGHT in Liver
do 
	TRAIT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/ukb-d-30780_irnt.vcf.gz"
	LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1"
	WEIGHT_FILE="/project2/compbio/predictdb/mashr_models/mashr_Liver.db"
	CONFIG_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b38.R"
	OUTNAME_E="LDL_"$WEIGHT"_expr"
	OUTNAME="LDL_"$WEIGHT"_ctwas"
	OUTDIR="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/ctwas_rss_LDR/"$WEIGHT
	job_name="ctwas_rss_LDR_"$WEIGHT
	job_out="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_rss_LDR_out/"$job_name".out"
	job_err="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_rss_LDR_out/"$job_name".err"
  
	FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
	if [ ! -f $FINAL_FILE ]
	then 
		sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_ctwas_rss_LDR.sbatch
	fi
done

