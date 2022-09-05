thin <- 0.1
ncore <- 10
ncore.rerun <- 1
prob_single <- 0.8
max_snp_region <- 20000
ld_regions <- "EUR"
ld_regions_version <- "b37"


locus_plot <- function(region_tag, rerun_ctwas = F, plot_eqtl = T, label="cTWAS"){
  region_tag1 <- unlist(strsplit(region_tag, "_"))[1]
  region_tag2 <- unlist(strsplit(region_tag, "_"))[2]
  
  a <- ctwas_res[ctwas_res$region_tag==region_tag,]
  
  regionlist <- readRDS(paste0(results_dir, "/", analysis_id, "_ctwas.regionlist.RDS"))
  region <- regionlist[[as.numeric(region_tag1)]][[region_tag2]]
  
  R_snp_info <- do.call(rbind, lapply(region$regRDS, function(x){data.table::fread(paste0(tools::file_path_sans_ext(x), ".Rvar"))}))
  
  if (isTRUE(rerun_ctwas)){
    ld_exprfs <- paste0(results_dir, "/", analysis_id, "_expr_chr", 1:22, ".expr.gz")
    temp_reg <- data.frame("chr" = paste0("chr",region_tag1), "start" = region$start, "stop" = region$stop)
    
    write.table(temp_reg, 
                #file= paste0(results_dir, "/", analysis_id, "_ctwas.temp.reg.txt") , 
                file= "temp_reg.txt",
                row.names=F, col.names=T, sep="\t", quote = F)
    
    load(paste0(results_dir, "/", analysis_id, "_expr_z_snp.Rd"))
    
    z_gene_temp <-  z_gene[z_gene$id %in% a$id[a$type=="gene"],]
    z_snp_temp <-  z_snp[z_snp$id %in% R_snp_info$id,]
    
    ctwas_rss(z_gene_temp, z_snp_temp, ld_exprfs, ld_pgenfs = NULL, 
              ld_R_dir = dirname(region$regRDS)[1],
              ld_regions_custom = "temp_reg.txt", thin = 1, 
              outputdir = ".", outname = "temp", ncore = 1, ncore.rerun = 1, prob_single = 0,
              group_prior = estimated_group_prior, group_prior_var = estimated_group_prior_var,
              estimate_group_prior = F, estimate_group_prior_var = F)
    
    
    a <- data.table::fread("temp.susieIrss.txt", header = T)
    
    rownames(z_snp_temp) <- z_snp_temp$id
    z_snp_temp <- z_snp_temp[a$id[a$type=="SNP"],]
    z_gene_temp <- z_gene_temp[a$id[a$type=="gene"],]
    
    a$z <- NA
    a$z[a$type=="SNP"] <- z_snp_temp$z
    a$z[a$type=="gene"] <- z_gene_temp$z
  }
  
  a$ifcausal <- 0
  focus <- a$id[a$type=="gene"][which.max(abs(a$z[a$type=="gene"]))]
  a$ifcausal <- as.numeric(a$id==focus)
  
  a$PVALUE <- (-log(2) - pnorm(abs(a$z), lower.tail=F, log.p=T))/log(10)
  
  R_gene <- readRDS(region$R_g_file)
  R_snp_gene <- readRDS(region$R_sg_file)
  R_snp <- as.matrix(Matrix::bdiag(lapply(region$regRDS, readRDS)))
  
  rownames(R_gene) <- region$gid
  colnames(R_gene) <- region$gid
  rownames(R_snp_gene) <- R_snp_info$id
  colnames(R_snp_gene) <- region$gid
  rownames(R_snp) <- R_snp_info$id
  colnames(R_snp) <- R_snp_info$id
  
  a$r2max <- NA
  
  a$r2max[a$type=="gene"] <- R_gene[focus,a$id[a$type=="gene"]]
  a$r2max[a$type=="SNP"] <- R_snp_gene[a$id[a$type=="SNP"],focus]
  
  r2cut <- 0.4
  colorsall <- c("#7fc97f", "#beaed4", "#fdc086")
  
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(1.5,1.5), respect = FALSE)
  par(mar = c(0, 4.1, 4.1, 2.1))
  plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 19, xlab=paste0("Chromosome ", region_tag1, " Position"),frame.plot=FALSE, col = "white", ylim= c(-0.1,1.1), ylab = "cTWAS PIP", xaxt = 'n')
  
  grid()
  points(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 21, xlab="Genomic position", bg = colorsall[1])
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$susie_pip[a$type == "SNP"  & a$r2max >r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$ifcausal == 1], a$susie_pip[a$type == "SNP" & a$ifcausal == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene"], a$susie_pip[a$type == "gene"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut], a$susie_pip[a$type == "gene"  & a$r2max > r2cut], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$ifcausal == 1], a$susie_pip[a$type == "gene" & a$ifcausal == 1], pch = 22, bg = "salmon", cex = 2)
  
  if (isTRUE(plot_eqtl)){
    for (cgene in a[a$type=="gene" & a$ifcausal == 1, ]$id){
      load(paste0(results_dir, "/",analysis_id, "_expr_chr", region_tag1, ".exprqc.Rd"))
      eqtls <- rownames(wgtlist[[cgene]])
      points(a[a$id %in% eqtls,]$pos, rep( -0.15, nrow(a[a$id %in% eqtls,])), pch = "|", col = "salmon", cex = 1.5)
    }
  }
  
  legend(min(a$pos), y= 1.1 ,c("Gene", "SNP"), pch = c(22,21), title="Shape Legend", bty ='n', cex=0.6, title.adj = 0)
  legend(min(a$pos), y= 0.7 ,c("Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = 19, col = c("salmon", "purple", colorsall[1]), title="Color Legend", bty ='n', cex=0.6, title.adj = 0)
  
  if (label=="cTWAS"){
    text(a$pos[a$id==focus], a$susie_pip[a$id==focus], labels=ctwas_gene_res$genename[ctwas_gene_res$id==focus], pos=3, cex=0.6)
  }
  
  par(mar = c(4.1, 4.1, 0.5, 2.1))
  plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 21, xlab=paste0("Chromosome ", region_tag1, " Position"), frame.plot=FALSE, bg = colorsall[1], ylab = "TWAS -log10(p value)", panel.first = grid(), ylim =c(0, max(a$PVALUE)*1.2))
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$PVALUE[a$type == "SNP"  & a$r2max > r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$ifcausal == 1], a$PVALUE[a$type == "SNP" & a$ifcausal == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene"], a$PVALUE[a$type == "gene"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut], a$PVALUE[a$type == "gene"  & a$r2max > r2cut], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$ifcausal == 1], a$PVALUE[a$type == "gene" & a$ifcausal == 1], pch = 22, bg = "salmon", cex = 2)
  abline(h=-log10(alpha/nrow(ctwas_gene_res)), col ="red", lty = 2)
  
  if (label=="TWAS"){
    text(a$pos[a$id==focus], a$PVALUE[a$id==focus], labels=ctwas_gene_res$genename[ctwas_gene_res$id==focus], pos=3, cex=0.6)
  }
}

locus_plot5 <- function(region_tag, rerun_ctwas = F, plot_eqtl = T, label="cTWAS", focus){
  region_tag1 <- unlist(strsplit(region_tag, "_"))[1]
  region_tag2 <- unlist(strsplit(region_tag, "_"))[2]
  
  a <- ctwas_res[ctwas_res$region_tag==region_tag,]
  
  regionlist <- readRDS(paste0(results_dir, "/", analysis_id, "_ctwas.regionlist.RDS"))
  region <- regionlist[[as.numeric(region_tag1)]][[region_tag2]]
  
  R_snp_info <- do.call(rbind, lapply(region$regRDS, function(x){data.table::fread(paste0(tools::file_path_sans_ext(x), ".Rvar"))}))
  
  if (isTRUE(rerun_ctwas)){
    ld_exprfs <- paste0(results_dir, "/", analysis_id, "_expr_chr", 1:22, ".expr.gz")
    temp_reg <- data.frame("chr" = paste0("chr",region_tag1), "start" = region$start, "stop" = region$stop)
    
    write.table(temp_reg, 
                #file= paste0(results_dir, "/", analysis_id, "_ctwas.temp.reg.txt") , 
                file= "temp_reg.txt",
                row.names=F, col.names=T, sep="\t", quote = F)
    
    load(paste0(results_dir, "/", analysis_id, "_expr_z_snp.Rd"))
    
    z_gene_temp <-  z_gene[z_gene$id %in% a$id[a$type=="gene"],]
    z_snp_temp <-  z_snp[z_snp$id %in% R_snp_info$id,]
    
    ctwas_rss(z_gene_temp, z_snp_temp, ld_exprfs, ld_pgenfs = NULL, 
              ld_R_dir = dirname(region$regRDS)[1],
              ld_regions_custom = "temp_reg.txt", thin = 1, 
              outputdir = ".", outname = "temp", ncore = 1, ncore.rerun = 1, prob_single = 0,
              group_prior = estimated_group_prior, group_prior_var = estimated_group_prior_var,
              estimate_group_prior = F, estimate_group_prior_var = F)
    
    
    a <- data.table::fread("temp.susieIrss.txt", header = T)
    
    rownames(z_snp_temp) <- z_snp_temp$id
    z_snp_temp <- z_snp_temp[a$id[a$type=="SNP"],]
    z_gene_temp <- z_gene_temp[a$id[a$type=="gene"],]
    
    a$z <- NA
    a$z[a$type=="SNP"] <- z_snp_temp$z
    a$z[a$type=="gene"] <- z_gene_temp$z
  }
  
  a$ifcausal <- 0
  focus <- a$id[which(a$genename==focus)]
  a$ifcausal <- as.numeric(a$id==focus)
  
  a$PVALUE <- (-log(2) - pnorm(abs(a$z), lower.tail=F, log.p=T))/log(10)
  
  R_gene <- readRDS(region$R_g_file)
  R_snp_gene <- readRDS(region$R_sg_file)
  R_snp <- as.matrix(Matrix::bdiag(lapply(region$regRDS, readRDS)))
  
  rownames(R_gene) <- region$gid
  colnames(R_gene) <- region$gid
  rownames(R_snp_gene) <- R_snp_info$id
  colnames(R_snp_gene) <- region$gid
  rownames(R_snp) <- R_snp_info$id
  colnames(R_snp) <- R_snp_info$id
  
  a$r2max <- NA
  
  a$r2max[a$type=="gene"] <- R_gene[focus,a$id[a$type=="gene"]]
  a$r2max[a$type=="SNP"] <- R_snp_gene[a$id[a$type=="SNP"],focus]
  
  r2cut <- 0.4
  colorsall <- c("#7fc97f", "#beaed4", "#fdc086")
  
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(1.5,1.5), respect = FALSE)
  par(mar = c(0, 4.1, 4.1, 2.1))
  plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 19, xlab=paste0("Chromosome ", region_tag1, " Position"),frame.plot=FALSE, col = "white", ylim= c(-0.1,1.1), ylab = "cTWAS PIP", xaxt = 'n')
  
  grid()
  points(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 21, xlab="Genomic position", bg = colorsall[1])
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$susie_pip[a$type == "SNP"  & a$r2max >r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$ifcausal == 1], a$susie_pip[a$type == "SNP" & a$ifcausal == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene"], a$susie_pip[a$type == "gene"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut], a$susie_pip[a$type == "gene"  & a$r2max > r2cut], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$ifcausal == 1], a$susie_pip[a$type == "gene" & a$ifcausal == 1], pch = 22, bg = "salmon", cex = 2)
  
  if (isTRUE(plot_eqtl)){
    for (cgene in a[a$type=="gene" & a$ifcausal == 1, ]$id){
      load(paste0(results_dir, "/",analysis_id, "_expr_chr", region_tag1, ".exprqc.Rd"))
      eqtls <- rownames(wgtlist[[cgene]])
      points(a[a$id %in% eqtls,]$pos, rep( -0.15, nrow(a[a$id %in% eqtls,])), pch = "|", col = "salmon", cex = 1.5)
    }
  }
  
  legend(max(a$pos)-0.2*(max(a$pos)-min(a$pos)), y= 1.1 ,c("Gene", "SNP"), pch = c(22,21), title="Shape Legend", bty ='n', cex=0.6, title.adj = 0)
  legend(max(a$pos)-0.2*(max(a$pos)-min(a$pos)), y= 0.7 ,c("Focal Gene", "R2 > 0.4", "R2 <= 0.4"), pch = 19, col = c("salmon", "purple", colorsall[1]), title="Color Legend", bty ='n', cex=0.6, title.adj = 0)
  
  if (label=="cTWAS"){
    text(a$pos[a$id==focus], a$susie_pip[a$id==focus], labels=ctwas_gene_res$genename[ctwas_gene_res$id==focus], pos=3, cex=0.6)
  }
  
  par(mar = c(4.1, 4.1, 0.5, 2.1))
  plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 21, xlab=paste0("Chromosome ", region_tag1, " Position"), frame.plot=FALSE, bg = colorsall[1], ylab = "TWAS -log10(p value)", panel.first = grid(), ylim =c(0, max(a$PVALUE)*1.2))
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$PVALUE[a$type == "SNP"  & a$r2max > r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$ifcausal == 1], a$PVALUE[a$type == "SNP" & a$ifcausal == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene"], a$PVALUE[a$type == "gene"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut], a$PVALUE[a$type == "gene"  & a$r2max > r2cut], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$ifcausal == 1], a$PVALUE[a$type == "gene" & a$ifcausal == 1], pch = 22, bg = "salmon", cex = 2)
  abline(h=-log10(alpha/nrow(ctwas_gene_res)), col ="red", lty = 2)
  
  if (label=="TWAS"){
    text(a$pos[a$id==focus], a$PVALUE[a$id==focus], labels=ctwas_gene_res$genename[ctwas_gene_res$id==focus], pos=3, cex=0.6)
  }
}

locus_plot4 <- function(region_tag, rerun_ctwas = F, plot_eqtl = T, label="cTWAS"){
  region_tag1 <- unlist(strsplit(region_tag, "_"))[1]
  region_tag2 <- unlist(strsplit(region_tag, "_"))[2]
  
  a <- ctwas_res[ctwas_res$region_tag==region_tag,]
  
  regionlist <- readRDS(paste0(results_dir, "/", analysis_id, "_ctwas.regionlist.RDS"))
  region <- regionlist[[as.numeric(region_tag1)]][[region_tag2]]
  
  R_snp_info <- do.call(rbind, lapply(region$regRDS, function(x){data.table::fread(paste0(tools::file_path_sans_ext(x), ".Rvar"))}))
  
  if (isTRUE(rerun_ctwas)){
    ld_exprfs <- paste0(results_dir, "/", analysis_id, "_expr_chr", 1:22, ".expr.gz")
    temp_reg <- data.frame("chr" = paste0("chr",region_tag1), "start" = region$start, "stop" = region$stop)
    
    write.table(temp_reg, 
                #file= paste0(results_dir, "/", analysis_id, "_ctwas.temp.reg.txt") , 
                file= "temp_reg.txt",
                row.names=F, col.names=T, sep="\t", quote = F)
    
    load(paste0(results_dir, "/", analysis_id, "_expr_z_snp.Rd"))
    
    z_gene_temp <-  z_gene[z_gene$id %in% a$id[a$type=="gene"],]
    z_snp_temp <-  z_snp[z_snp$id %in% R_snp_info$id,]
    
    ctwas_rss(z_gene_temp, z_snp_temp, ld_exprfs, ld_pgenfs = NULL, 
              ld_R_dir = dirname(region$regRDS)[1],
              ld_regions_custom = "temp_reg.txt", thin = 1, 
              outputdir = ".", outname = "temp", ncore = 1, ncore.rerun = 1, prob_single = 0,
              group_prior = estimated_group_prior, group_prior_var = estimated_group_prior_var,
              estimate_group_prior = F, estimate_group_prior_var = F)
    
    
    a <- data.table::fread("temp.susieIrss.txt", header = T)
    
    rownames(z_snp_temp) <- z_snp_temp$id
    z_snp_temp <- z_snp_temp[a$id[a$type=="SNP"],]
    z_gene_temp <- z_gene_temp[a$id[a$type=="gene"],]
    
    a$z <- NA
    a$z[a$type=="SNP"] <- z_snp_temp$z
    a$z[a$type=="gene"] <- z_gene_temp$z
  }
  
  a$ifcausal <- 0
  focus <- a$id[a$type=="gene"][which.max(abs(a$z[a$type=="gene"]))]
  a$ifcausal <- as.numeric(a$id==focus)
  
  a$PVALUE <- (-log(2) - pnorm(abs(a$z), lower.tail=F, log.p=T))/log(10)
  
  R_gene <- readRDS(region$R_g_file)
  R_snp_gene <- readRDS(region$R_sg_file)
  R_snp <- as.matrix(Matrix::bdiag(lapply(region$regRDS, readRDS)))
  
  rownames(R_gene) <- region$gid
  colnames(R_gene) <- region$gid
  rownames(R_snp_gene) <- R_snp_info$id
  colnames(R_snp_gene) <- region$gid
  rownames(R_snp) <- R_snp_info$id
  colnames(R_snp) <- R_snp_info$id
  
  a$r2max <- NA
  
  a$r2max[a$type=="gene"] <- R_gene[focus,a$id[a$type=="gene"]]
  a$r2max[a$type=="SNP"] <- R_snp_gene[a$id[a$type=="SNP"],focus]
  
  r2cut <- 0.4
  colorsall <- c("#7fc97f", "#beaed4", "#fdc086")
  
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(1.5,1.5), respect = FALSE)
  par(mar = c(0, 4.1, 4.1, 2.1))
  plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 19, xlab=paste0("Chromosome ", region_tag1, " Position"),frame.plot=FALSE, col = "white", ylim= c(-0.1,1.1), ylab = "cTWAS PIP", xaxt = 'n')
  
  grid()
  points(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 21, xlab="Genomic position", bg = colorsall[1])
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$susie_pip[a$type == "SNP"  & a$r2max >r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$ifcausal == 1], a$susie_pip[a$type == "SNP" & a$ifcausal == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene"], a$susie_pip[a$type == "gene"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut], a$susie_pip[a$type == "gene"  & a$r2max > r2cut], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$ifcausal == 1], a$susie_pip[a$type == "gene" & a$ifcausal == 1], pch = 22, bg = "salmon", cex = 2)
  
  if (isTRUE(plot_eqtl)){
    for (cgene in a[a$type=="gene" & a$ifcausal == 1, ]$id){
      load(paste0(results_dir, "/",analysis_id, "_expr_chr", region_tag1, ".exprqc.Rd"))
      eqtls <- rownames(wgtlist[[cgene]])
      points(a[a$id %in% eqtls,]$pos, rep( -0.15, nrow(a[a$id %in% eqtls,])), pch = "|", col = "salmon", cex = 1.5)
    }
  }
  
  #legend(min(a$pos), y= 1.1 ,c("Gene", "SNP"), pch = c(22,21), title="Shape Legend", bty ='n', cex=0.6, title.adj = 0)
  #legend(min(a$pos), y= 0.7 ,c("Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = 19, col = c("salmon", "purple", colorsall[1]), title="Color Legend", bty ='n', cex=0.6, title.adj = 0)
  
  legend(max(a$pos)-0.2*(max(a$pos)-min(a$pos)), y= 1.1 ,c("Gene", "SNP"), pch = c(22,21), title="Shape Legend", bty ='n', cex=0.6, title.adj = 0)
  legend(max(a$pos)-0.2*(max(a$pos)-min(a$pos)), y= 0.7 ,c("Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = 19, col = c("salmon", "purple", colorsall[1]), title="Color Legend", bty ='n', cex=0.6, title.adj = 0)
  
  if (label=="cTWAS"){
    text(a$pos[a$id==focus], a$susie_pip[a$id==focus], labels=ctwas_gene_res$genename[ctwas_gene_res$id==focus], pos=3, cex=0.6)
  }
  
  par(mar = c(4.1, 4.1, 0.5, 2.1))
  plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 21, xlab=paste0("Chromosome ", region_tag1, " Position"), frame.plot=FALSE, bg = colorsall[1], ylab = "TWAS -log10(p value)", panel.first = grid(), ylim =c(0, max(a$PVALUE)*1.2))
  points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$PVALUE[a$type == "SNP"  & a$r2max > r2cut], pch = 21, bg = "purple")
  points(a$pos[a$type=="SNP" & a$ifcausal == 1], a$PVALUE[a$type == "SNP" & a$ifcausal == 1], pch = 21, bg = "salmon")
  points(a$pos[a$type=="gene"], a$PVALUE[a$type == "gene"], pch = 22, bg = colorsall[1], cex = 2)
  points(a$pos[a$type=="gene" & a$r2max > r2cut], a$PVALUE[a$type == "gene"  & a$r2max > r2cut], pch = 22, bg = "purple", cex = 2)
  points(a$pos[a$type=="gene" & a$ifcausal == 1], a$PVALUE[a$type == "gene" & a$ifcausal == 1], pch = 22, bg = "salmon", cex = 2)
  abline(h=-log10(alpha/nrow(ctwas_gene_res)), col ="red", lty = 2)
  
  if (label=="TWAS"){
    text(a$pos[a$id==focus], a$PVALUE[a$id==focus], labels=ctwas_gene_res$genename[ctwas_gene_res$id==focus], pos=3, cex=0.6)
  }
}
