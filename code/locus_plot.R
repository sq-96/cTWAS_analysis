locus_plot <- function(region_tag, xlim=NULL, return_table=F, focus=NULL, label_panel="TWAS", label_genes=NULL, label_pos=NULL, plot_eqtl=NULL, rerun_ctwas=F, rerun_load_only=F, legend_side="right", legend_panel="cTWAS", twas_ymax=NULL){
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
    
        if (!rerun_load_only){
            ctwas::ctwas_rss(z_gene_temp, z_snp_temp, ld_exprfs, ld_pgenfs = NULL, 
                      ld_R_dir = dirname(region$regRDS)[1],
                      ld_regions_custom = "temp_reg.txt", thin = 1, 
                      outputdir = ".", outname = "temp", ncore = 1, ncore.rerun = 1, prob_single = 0,
                      group_prior = estimated_group_prior, group_prior_var = estimated_group_prior_var,
                      estimate_group_prior = F, estimate_group_prior_var = F)
        }
    
        a_bkup <- a         
        a <- as.data.frame(data.table::fread("temp.susieIrss.txt", header = T))
    
        rownames(z_snp_temp) <- z_snp_temp$id
        z_snp_temp <- z_snp_temp[a$id[a$type=="SNP"],]
        z_gene_temp <- z_gene_temp[a$id[a$type=="gene"],]
    
        a$genename <- NA
        a$gene_type <- NA

        a[a$type=="gene",c("genename", "gene_type")] <- a_bkup[match(a$id[a$type=="gene"], a_bkup$id),c("genename","gene_type")]
    
        a$z <- NA
        a$z[a$type=="SNP"] <- z_snp_temp$z
        a$z[a$type=="gene"] <- z_gene_temp$z
    }
  
    a_pos_bkup <- a$pos
    a$pos[a$type=="gene"] <- G_list$tss[match(sapply(a$genename[a$type=="gene"], function(x){unlist(strsplit(x, "[.]"))[1]}) ,G_list$hgnc_symbol)]
    a$pos[is.na(a$pos)] <- a_pos_bkup[is.na(a$pos)]
    a$pos <- a$pos/1000000
  
    if (!is.null(xlim)){
    
        if (is.na(xlim[1])){
            xlim[1] <- min(a$pos)
        }
    
        if (is.na(xlim[2])){
            xlim[2] <- max(a$pos)
        }
    
        a <- a[a$pos>=xlim[1] & a$pos<=xlim[2],,drop=F]
    }
  
    if (is.null(focus)){
        focus <- a$id[which.max(abs(a$z)[a$type=="gene"])]
    }
  
    if (is.null(label_genes)){
        label_genes <- focus
    }
  
    if (is.null(label_pos)){
        label_pos <- rep(3, length(label_genes))
    }
  
    if (is.null(plot_eqtl)){
        plot_eqtl <- focus
    }
  
    focus <- a$id[which(a$id==focus)]
    a$focus <- 0
    a$focus <- as.numeric(a$id==focus)
    
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
  
    start <- min(a$pos)
    end <- max(a$pos)
  
    layout(matrix(1:3, ncol = 1), widths = 1, heights = c(1.5,1.75,0.75), respect = FALSE)
  
    par(mar = c(0, 4.1, 0, 2.1))
  
    if (is.null(twas_ymax)){
        twas_ymax <- max(a$PVALUE)*1.1
    }
  
    plot(a$pos[a$type=="SNP"], a$PVALUE[a$type == "SNP"], pch = 21, xlab=paste0("Chromosome ", region_tag1, " position (Mb)"), frame.plot=FALSE, bg = colorsall[1], ylab = "-log10(p value)", panel.first = grid(), ylim =c(0, twas_ymax), xaxt = 'n', xlim=c(start, end))
  
    abline(h=-log10(alpha/nrow(ctwas_gene_res)), col ="red", lty = 2)
    points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$PVALUE[a$type == "SNP"  & a$r2max > r2cut], pch = 21, bg = "purple")
    points(a$pos[a$type=="SNP" & a$focus == 1], a$PVALUE[a$type == "SNP" & a$focus == 1], pch = 21, bg = "salmon")
    points(a$pos[a$type=="gene" & a$group=="Expression"], a$PVALUE[a$type == "gene" & a$group=="Expression"], pch = 22, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Splicing"], a$PVALUE[a$type == "gene" & a$group=="Splicing"], pch = 23, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Methylation"], a$PVALUE[a$type == "gene" & a$group=="Methylation"], pch = 24, bg = colorsall[1], cex = 2)

    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Expression"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Expression"], pch = 22, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Splicing"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Splicing"], pch = 23, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Methylation"], a$PVALUE[a$type == "gene"  & a$r2max > r2cut & a$group=="Methylation"], pch = 24, bg = "purple", cex = 2)

    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Expression"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Expression"], pch = 22, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Splicing"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Splicing"], pch = 23, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Methylation"], a$PVALUE[a$type == "gene" & a$focus == 1 & a$group=="Methylation"], pch = 24, bg = "salmon", cex = 2)
  
    if (legend_panel=="TWAS"){
        x_pos <- ifelse(legend_side=="right", max(a$pos)-0.2*(max(a$pos)-min(a$pos)), min(a$pos))
        legend(x_pos, y= twas_ymax*0.95, c("Expression","Splicing","Methylation","SNP","Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = c(22,23,24,21,19,19,19), col = c("black","black","black","black", "salmon", "purple", colorsall[1]), cex=0.7, title.adj = 0)
    }
  
    if (label_panel=="TWAS" | label_panel=="both"){
        for (i in 1:length(label_genes)){
            text(a$pos[a$id==label_genes[i]], a$PVALUE[a$id==label_genes[i]], labels=a$genename[a$id==label_genes[i]], pos=label_pos[i], cex=0.7)
        }
    }
  
    #par(mar = c(0.25, 4.1, 0.25, 2.1))
  
    #plot(NA, xlim = c(start, end), ylim = c(0, length(plot_eqtl)), frame.plot = F, axes = F, xlab = NA, ylab = NA)
  
    #for (i in 1:length(plot_eqtl)){
    #    cgene <- a$id[which(a$id==plot_eqtl[i])]
    #    load(paste0(results_dir, "/",analysis_id, "_expr_chr", region_tag1, ".exprqc.Rd"))
    #    eqtls <- rownames(wgtlist[[cgene]])
    #    eqtl_pos <- a$pos[a$id %in% eqtls]
    
        #col="grey"
    #    col="#c6e8f0"
    
    #    rect(start, length(plot_eqtl)+1-i-0.8, end, length(plot_eqtl)+1-i-0.2, col = col, border = T, lwd = 1)
  
    #    if (length(eqtl_pos)>0){
    #        for (j in 1:length(eqtl_pos)){
    #            segments(x0=eqtl_pos[j], x1=eqtl_pos[j], y0=length(plot_eqtl)+1-i-0.2, length(plot_eqtl)+1-i-0.8, lwd=1.5)  
    #        }
    #    }
    #}
  
    #text(start, length(plot_eqtl)-(1:length(plot_eqtl))+0.5,  
    #    labels = paste0(plot_eqtl, " eQTL"), srt = 0, pos = 2, xpd = TRUE, cex=0.7)
  
    par(mar = c(4.1, 4.1, 0, 2.1))
  
    plot(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 19, xlab=paste0("Chromosome ", region_tag1, " position (Mb)"),frame.plot=FALSE, col = "white", ylim= c(0,1.1), ylab = "cTWAS PIP", xlim = c(start, end))
  
    grid()
    points(a$pos[a$type=="SNP"], a$susie_pip[a$type == "SNP"], pch = 21, xlab="Genomic position", bg = colorsall[1])
    points(a$pos[a$type=="SNP" & a$r2max > r2cut], a$susie_pip[a$type == "SNP"  & a$r2max >r2cut], pch = 21, bg = "purple")
    points(a$pos[a$type=="SNP" & a$focus == 1], a$susie_pip[a$type == "SNP" & a$focus == 1], pch = 21, bg = "salmon")
    points(a$pos[a$type=="gene" & a$group=="Expression"], a$susie_pip[a$type == "gene" & a$group=="Expression"], pch = 22, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Splicing"], a$susie_pip[a$type == "gene" & a$group=="Splicing"], pch = 23, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$group=="Methylation"], a$susie_pip[a$type == "gene" & a$group=="Methylation"], pch = 24, bg = colorsall[1], cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Expression"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Expression"], pch = 22, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Splicing"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Splicing"], pch = 23, bg = "purple", cex = 2)
    points(a$pos[a$type=="gene" & a$r2max > r2cut & a$group=="Methylation"], a$susie_pip[a$type == "gene"  & a$r2max > r2cut & a$group=="Methylation"], pch = 24, bg = "purple", cex = 2)

    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Expression"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Expression"], pch = 22, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Splicing"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Splicing"], pch = 23, bg = "salmon", cex = 2)
    points(a$pos[a$type=="gene" & a$focus == 1 & a$group=="Methylation"], a$susie_pip[a$type == "gene" & a$focus == 1 & a$group=="Methylation"], pch = 24, bg = "salmon", cex = 2)
  
    if (legend_panel=="cTWAS"){
        x_pos <- ifelse(legend_side=="right", max(a$pos)-0.2*(max(a$pos)-min(a$pos)), min(a$pos))
        legend(x_pos, y= 1 ,c("Expression","Splicing","Methylation","SNP","Lead TWAS Gene", "R2 > 0.4", "R2 <= 0.4"), pch = c(22,23,24,21,19,19,19), col = c("black","black","black", "black", "salmon", "purple", colorsall[1]), cex=0.7, title.adj = 0)
    }


  
    if (label_panel=="cTWAS" | label_panel=="both"){
        for (i in 1:length(label_genes)){
        text(a$pos[a$id==label_genes[i]], a$susie_pip[a$id==label_genes[i]], labels=a$genename[a$id==label_genes[i]], pos=label_pos[i], cex=0.7)
        }
    }
  
    if (return_table){
        return(a)
    }
}


library(Gviz)

gene_track <- function(a, label_pos=NULL){
  chr <- unique(a$chrom)
  start <- min(a$pos)*1000000
  end <- max(a$pos)*1000000
  
  biomTrack <- BiomartGeneRegionTrack(chromosome = chr,
                                      start = start,
                                      end = end,
                                      name = "ENSEMBL",
                                      biomart = ensembl,
                                      filters=list(biotype="protein_coding"))
  
  
  biomTrack <- as(biomTrack, "GeneRegionTrack")
  biomTrack <- biomTrack[biomTrack@range@elementMetadata@listData$feature %in% c("protein_coding", "utr3", "utr5")]
  
  if (isTRUE(label_pos=="above")){
    displayPars(biomTrack)$just.group <- "above"
  }
  
  grid.newpage()
  
  plotTracks(biomTrack, collapseTranscripts = "meta", transcriptAnnotation = "symbol", from=start, to=end, panel.only=T, add=F)
}














read_pvar <- function(pvarf){
  pvardt <- data.table::fread(pvarf, skip = "#CHROM")
  pvardt <- dplyr::rename(pvardt, chrom = "#CHROM", pos = "POS", 
                          alt = "ALT", ref = "REF", id = "ID")
  pvardt <- pvardt[, c("chrom", "id", "pos", "alt", "ref")]
  pvardt
}

prep_pvar <- function (pgenf, outputdir = getwd()){
  if (file_ext(pgenf) == "pgen") {
    pvarf <- paste0(file_path_sans_ext(pgenf), ".pvar")
    pvarf2 <- paste0(outputdir, basename(file_path_sans_ext(pgenf)), 
                     ".hpvar")
    firstl <- read.table(file = pvarf, header = F, comment.char = "", 
                         nrows = 1, stringsAsFactors = F)
    if (substr(firstl[1, 1], 1, 1) == "#") {
      pvarfout <- pvarf
    }
    else {
      pvarfout <- pvarf2
      if (!file.exists(pvarf2)) {
        pvar <- data.table::fread(pvarf, header = F)
        if (ncol(pvar) == 6) {
          colnames(pvar) <- c("#CHROM", "ID", "CM", "POS", 
                              "ALT", "REF")
        }
        else if (ncol(pvar) == 5) {
          colnames(pvar) <- c("#CHROM", "ID", "POS", 
                              "ALT", "REF")
        }
        else {
          stop(".pvar file has incorrect format")
        }
        data.table::fwrite(pvar, file = pvarf2, sep = "\t", 
                           quote = F)
      }
    }
  }
  else if (file_ext(pgenf) == "bed") {
    pvarf <- paste0(file_path_sans_ext(pgenf), ".bim")
    pvarf2 <- file.path(outputdir, paste0(basename(file_path_sans_ext(pgenf)), 
                                          ".hbim"))
    if (!file.exists(pvarf2)) {
      pvar <- data.table::fread(pvarf, header = F)
      colnames(pvar) <- c("#CHROM", "ID", "CM", "POS", 
                          "ALT", "REF")
      data.table::fwrite(pvar, file = pvarf2, sep = "\t", 
                         quote = F)
    }
    pvarfout <- pvarf2
  }
  else {
    stop("Unrecognized genotype input format")
  }
  pvarfout
}

prep_pgen <- function(pgenf, pvarf){
  pvar <- pgenlibr::NewPvar(pvarf)
  if (file_ext(pgenf) == "pgen") {
    pgen <- pgenlibr::NewPgen(pgenf, pvar = pvar)
  }
  else if (file_ext(pgenf) == "bed") {
    famf <- paste0(file_path_sans_ext(pgenf), ".fam")
    fam <- data.table::fread(famf, header = F)
    raw_s_ct <- nrow(fam)
    pgen <- pgenlibr::NewPgen(pgenf, pvar = pvar, raw_sample_ct = raw_s_ct)
  }
  else {
    stop("unrecognized input")
  }
  pgen
}

read_pgen <- function(pgen, variantidx = NULL, meanimpute = F){
  if (is.null(variantidx)) {
    variantidx <- 1:pgenlibr::GetVariantCt(pgen)
  }
  pgenlibr::ReadList(pgen, variant_subset = variantidx, meanimpute = meanimpute)
}
