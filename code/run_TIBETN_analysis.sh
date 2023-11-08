#!/bin/bash
for TRAIT in LB MC OxHb Pulse
do
	for WEIGHT in Artery_Coronary Artery_Tibial Heart_Atrial_Appendage Heart_Left_Ventricle Kidney_Cortex Liver Lung Pancreas Whole_Blood
	do 
		TRAIT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/"$TRAIT".txt"
		LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1"
		WEIGHT_FILE="/project2/compbio/predictdb/mashr_models/mashr_"$WEIGHT".db"
		CONFIG_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b38.R"
		OUTNAME_E=$TRAIT"_"$WEIGHT"_expr"
		OUTNAME=$TRAIT"_"$WEIGHT"_ctwas"
		OUTDIR="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/TIBETN/"$TRAIT"/"$WEIGHT
		job_name=$TRAIT"_"$WEIGHT
		job_out="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/TIBETN_out/"$job_name".out"
		job_err="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/TIBETN_out/"$job_name".err"
  
		FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
		if [ ! -f $FINAL_FILE ]
		then 
			sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_TIBETN_analysis.sbatch
		fi
	done
done

