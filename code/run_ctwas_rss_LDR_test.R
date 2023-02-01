library(ctwas)

args = commandArgs(trailingOnly=TRUE)

args <- c("/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB/ukb-d-30000_irnt.vcf.gz",
          "/project2/mstephens/wcrouse/UKB_LDR_0.1",
          "/project2/compbio/predictdb/mashr_models/mashr_Liver.db;/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/mashr_sqtl/sqtl/mashr/mashr_Liver_Splicing.db;/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/mqtl/WholeBlood.db",
          "/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/code/ctwas_config_b38.R",
          "WhiteBlood_WholeBlood_expr",
          "WhiteBlood_WholeBlood_ctwas",
          "/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/WhiteBlood_E_S_M/WholeBlood")

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

####################

outputdir <- args[7]

dir.create(outputdir, showWarnings=F, recursive=T)

z_snp_stem <- unlist(strsplit(rev(unlist(strsplit(args[1], "/")))[1],"[.]"))[1]
z_snp_outfile <- paste0(outputdir, "/", z_snp_stem, ".RDS")

if (file.exists(z_snp_outfile)){
  z_snp <- readRDS(z_snp_outfile)
} else {
  z_snp <- VariantAnnotation::readVcf(args[1])
  z_snp <- as.data.frame(gwasvcf::vcf_to_tibble(z_snp))
  z_snp$Z <- z_snp$ES/z_snp$SE
  z_snp <- z_snp[,c("rsid", "ALT", "REF", "Z", "SS")]
  colnames(z_snp) <- c("id", "A1", "A2", "z", "ss")
  z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),] #drop multiallelic variants (id not unique)
  saveRDS(z_snp, file=z_snp_outfile)
}

ld_R_dir <- args[2]

weight <- args[3]
weight <- unlist(strsplit(weight, ";"))

outname.e <- args[5]
outname <- args[6]

source(args[4]) # config

# get gene z score
if (file.exists(paste0(outputdir, "/", outname.e, "_z_gene.Rd"))){
  ld_exprfs <- paste0(outputdir, "/", outname.e, "_chr", 1:22, ".expr.gz")
  load(file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
  load(file = paste0(outputdir, "/", outname.e, "_z_snp.Rd"))
} else {
  res <- impute_expr_z(z_snp, weight = weight, ld_R_dir = ld_R_dir,
                       method = NULL, outputdir = outputdir, outname = outname.e,
                       harmonize_z = F, harmonize_wgt = F,
                       strand_ambig_action_z = "none", recover_strand_ambig_wgt = T)
  z_gene <- res$z_gene
  ld_exprfs <- res$ld_exprfs
  z_snp <- res$z_snp
  
  save(z_gene, file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
  save(z_snp, file = paste0(outputdir, "/", outname.e, "_z_snp.Rd"))
}

####################
#set the type for each gene (e.g. tissue)
z_gene$type <- sapply(z_gene$id, function(x){paste(unlist(strsplit(unlist(strsplit(x, "[|]"))[2],"_"))[-1], collapse="_") })

####################


regionfile <- system.file("extdata", "ldetect", paste0(ld_regions, "." , ld_regions_version, ".bed"), package="ctwas")

temp_reg <- as.data.frame(data.table::fread(regionfile))
temp_reg <- temp_reg[1:20,]

write.table(temp_reg, 
            file= paste0(outputdir,"/","temp_reg.txt"),
            row.names=F, col.names=T, sep="\t", quote = F)

regionfile <- paste0(outputdir,"/","temp_reg.txt")

#####

ncore <- 4

####################

#run ctwas_rss
#ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = NULL, ld_R_dir = ld_R_dir, ld_regions = ld_regions, ld_regions_version = ld_regions_version, thin = thin, max_snp_region = max_snp_region, outputdir = outputdir, outname = outname, ncore = ncore, ncore.rerun = ncore.rerun, prob_single = prob_single)
ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = NULL, ld_R_dir = ld_R_dir, ld_regions = ld_regions, ld_regions_version = ld_regions_version, ld_regions_custom = regionfile, thin = thin, max_snp_region = max_snp_region, outputdir = outputdir, outname = outname, ncore = ncore, ncore.rerun = ncore.rerun, prob_single = prob_single)


sessionInfo()
