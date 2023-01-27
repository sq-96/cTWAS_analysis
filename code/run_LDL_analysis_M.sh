#!/bin/bash

for WEIGHT in BreastMammary ColonTransverse KidneyCortex Lung MuscleSkeletal Ovary Prostate Testis WholeBlood
do 
	TRAIT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/ukb-d-30000_irnt.vcf.gz"
	LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1"
	WEIGHT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/mqtl/"$WEIGHT".db"
	CONFIG_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b38.R"
	OUTNAME_E="White_Blood_"$WEIGHT"_expr"
	OUTNAME="White_Blood_"$WEIGHT"_ctwas"
	OUTDIR="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/White_Blood_M/"$WEIGHT
	job_name="White_Blood_"$WEIGHT
	job_out="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/White_Blood_M_out/"$job_name".out"
	job_err="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/White_Blood_M_out/"$job_name".err"
  
	FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
	if [ ! -f $FINAL_FILE ]
	then 
		sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_LDL_analysis_M.sbatch
	fi
done

