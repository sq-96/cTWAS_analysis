#!/bin/bash
TRAIT_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/ukb-d-30780_irnt.vcf.gz"
LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1"
WEIGHT_FILE="/project2/compbio/predictdb/mashr_models/mashr_Liver.db;/project2/compbio/predictdb/mashr_models/mashr_Adipose_Subcutaneous.db;/project2/compbio/predictdb/mashr_models/mashr_Brain_Cerebellum.db;/project2/compbio/predictdb/mashr_models/mashr_Adipose_Visceral_Omentum.db;/project2/compbio/predictdb/mashr_models/mashr_Whole_Blood.db;/project2/compbio/predictdb/mashr_models/mashr_Lung.db;/project2/compbio/predictdb/mashr_models/mashr_Artery_Tibial.db;/project2/compbio/predictdb/mashr_models/mashr_Heart_Left_Ventricle.db;/project2/compbio/predictdb/mashr_models/mashr_Stomach.db;/project2/compbio/predictdb/mashr_models/mashr_Pancreas.db"
CONFIG_FILE="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b38.R"
OUTNAME_E="LDL_Liver_expr"
OUTNAME="LDL_Liver_ctwas"
OUTDIR="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/LDL_multi_tissue/"
job_name="LDL_multi_tissue"
job_out="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/LDL_multi_tissue_out/"$job_name".out"
job_err="/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/LDL_multi_tissue_out/"$job_name".err"
  
FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
if [ ! -f $FINAL_FILE ]
then 
	sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_LDL_analysis_multi_tissue.sbatch
fi

