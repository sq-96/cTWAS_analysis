library(ctwas)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args) <7) {
  stop(" 7 arguments:
       * zscore file name
       * ld R directory
       * weight
       * config file name (.R)
       * out expr z file name
       * out file name
       * outputdir", call.=FALSE)
}

print(args)

outputdir <- args[7]

dir.create(outputdir, showWarnings=F, recursive=T)

z_snp_stem <- unlist(strsplit(rev(unlist(strsplit(args[1], "/")))[1],"[.]"))[1]
z_snp_outfile <- paste0(outputdir, "/", z_snp_stem, ".RDS")

if (file.exists(z_snp_outfile)){
  z_snp <- readRDS(z_snp_outfile)
} else {
  z_snp <- fread("/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/OxHb.txt",fill=TRUE)
  z_snp <- na.omit(z_snp)
  z_snp$SNP <- sapply(z_snp$rs, function(x){unlist(strsplit(x, split="[:]"))[1]})
  z_snp$Z <- z_snp$beta/z_snp$se
  z_snp$SS <- 966
  z_snp <- z_snp[,c("SNP", "allele1", "allele0", "Z", "SS")]
  colnames(z_snp) <- c("id", "A1", "A2", "z", "ss")
  z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),] #drop multiallelic variants (id not unique)
  saveRDS(z_snp, file=z_snp_outfile)
}

ld_R_dir <- args[2]
weight <- args[3]
outname.e <- args[5]
outname <- args[6]

source(args[4]) # config

#preharmonize snp z score
if (file.exists(paste0(outputdir, "/", outname, "_z_snp.Rd"))){
  load(file = paste0(outputdir, "/", outname, "_z_snp.Rd"))
} else {
  res <- ctwas:::preharmonize_z_ld(z_snp=z_snp, 
                                   ld_R_dir=ld_R_dir, 
                                   outputdir=outputdir,
                                   outname=outname,
                                   harmonize_z=T, 
                                   strand_ambig_action_z="drop")
  z_snp <- res$z_snp
  save(z_snp, file = paste0(outputdir, "/", outname, "_z_snp.Rd"))
  rm(res)
}

#impute gene z-scores for both sets of prediction weights by chromosome
for (i in 1:22){
  if (!file.exists(paste0(outputdir, "/", outname, "_chr", i, ".expr.gz"))){
    res <- impute_expr_z(z_snp = z_snp,
                         weight = weight,
                         ld_R_dir = ld_R_dir,
                         outputdir = outputdir,
                         outname = outname,
                         harmonize_z = F,
                         harmonize_wgt = T,
                         strand_ambig_action_wgt="drop",
                         ncore=10, 
                         chrom=i)
  }
}

#combine the imputed gene z-scores
ld_exprfs <- paste0(outputdir, "/", outname, "_chr", 1:22, ".expr.gz")
z_gene <- list()
for (i in 1:22){
  load(paste0(outputdir, "/", outname, "_chr", i, ".exprqc.Rd"))
  z_gene[[i]] <- z_gene_chr
}
z_gene <- do.call(rbind, z_gene)
rownames(z_gene) <- NULL

rm(qclist, wgtlist, z_gene_chr)

#run ctwas_rss parameter estimation
if (file.exists(paste0(outputdir, "/", outname, ".s2.susieIrssres.Rd"))){
  print("skip parameter estimation")
  load(paste0(outputdir, "/", outname, ".s2.susieIrssres.Rd"))

  group_prior_rec <- group_prior_rec[,ncol(group_prior_rec)]
  group_prior_var_rec <- group_prior_var_rec[,ncol(group_prior_var_rec)]

  if (!file.exists(paste0(outputdir, "/", outname, ".susieIrss.txt"))){
    print("start fine mapping")
    ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = NULL, ld_R_dir = ld_R_dir, 
              ld_regions = ld_regions, ld_regions_version = ld_regions_version, 
              thin = 0.1, max_snp_region = max_snp_region, outputdir = outputdir, 
              outname = outname, ncore = 10, ncore.rerun = 1, prob_single = prob_single,
              merge=F, fine_map=T,
              group_prior = group_prior_rec, 
              group_prior_var = group_prior_var_rec, 
              estimate_group_prior = F, 
              estimate_group_prior_var = F,
              ncore_LDR=5)
  }
} else {
  print("parameter estimation")
  ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = NULL, ld_R_dir = ld_R_dir, 
            ld_regions = ld_regions, ld_regions_version = ld_regions_version, thin = 0.1, 
            max_snp_region = max_snp_region, outputdir = outputdir, outname = outname, 
            ncore = 10, ncore.rerun = 1, prob_single = prob_single,
            merge=F, 
            fine_map=F,
            ncore_LDR=5)
}

sessionInfo()