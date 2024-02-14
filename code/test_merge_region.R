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

ld_R_dir <- args[2]
weight <- args[3]
outname.e <- args[5]
outname <- args[6]

source(args[4]) # config

ld_exprfs <- paste0(outputdir, "/", outname, "_chr", 1:22, ".expr.gz")
ld_exprvarfs <- sapply(ld_exprfs, ctwas:::prep_exprvar)
exprvarfs <- ld_exprvarfs


z_gene <- list()
for (i in 1:22){
  load(paste0(outputdir, "/", outname, "_chr", i, ".exprqc.Rd"))
  z_gene[[i]] <- z_gene_chr
}
z_gene <- do.call(rbind, z_gene)
rownames(z_gene) <- NULL
save(z_gene, file = paste0(outputdir, "/", outname, "_z_gene.Rd"))


load(paste0(outputdir, "/", outname, "_z_snp.Rd"))
load(paste0(outputdir, "/", outname, "_z_gene.Rd"))
z_snp$type <- "SNP"
z_gene$type <- "gene"

ld_regions = "EUR"
ld_regions_version = "b38"
regionfile <- system.file("extdata", "ldetect",
                          paste0(ld_regions, "." , ld_regions_version, ".bed"), package="ctwas")

zdf <- rbind(z_snp[, c("id", "z", "type")], z_gene[, c("id", "z", "type")])

ld_Rfs <- ctwas:::write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)

regionlist <- ctwas:::index_regions(regionfile = regionfile,
                              exprvarfs = ld_exprvarfs,
                              pvarfs = NULL,
                              ld_Rfs = ld_Rfs,
                              select = zdf$id,
                              thin = 0.1, minvar = 2,
                              outname = outname,
                              outputdir = outputdir,
                              merge = FALSE,
                              ncore = 5) # susie_rss can't take 1 var.
  
saveRDS(regionlist, file=paste0(outputdir, "/", outname, ".regionlist.RDS"))