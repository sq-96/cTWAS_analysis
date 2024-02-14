library(ggplot2)
library(ggpubr)

get_sim_joint_res <- function(results_dir,runtag,configtag,simutags,thin,sample_size,PIP_threshold){
  results_df <- data.frame()
  
  for (i in 1:length(simutags)){
    simutag <- simutags[i]
    
    #load genes with true simulated effect
    load(paste0(results_dir, runtag, "_simu", simutag, "-pheno.Rd"))
    true_genes <- unlist(sapply(1:22, function(x){phenores$batch[[x]]$id.cgene}))
    true_genes_combined <- unique(sapply(true_genes, function(x){unlist(strsplit(x, "[|]"))[1]}))
    
    #load cTWAS results
    ctwas_res <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.susieIrss.txt"))
    ctwas_gene_res <- ctwas_res[ctwas_res$type!="SNP",]
    
    #number of causal genes
    n_causal <- length(true_genes)
    n_causal_combined <- length(true_genes_combined)
    
    #number of gene+tissue combinations with cTWAS PIP > threshold
    n_ctwas_genes <- sum(ctwas_gene_res$susie_pip > PIP_threshold)
    
    #number of cTWAS genes that are causal
    n_causal_detected <- sum(ctwas_gene_res$id[ctwas_gene_res$susie_pip > PIP_threshold] %in% true_genes)
    
    #collapse gene+tissues to genes and compute combined PIP
    ctwas_gene_res$gene <- sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x,"[|]"))[1]})
    ctwas_gene_res_combined <- aggregate(ctwas_gene_res$susie_pip, by=list(ctwas_gene_res$gene), FUN=sum)
    colnames(ctwas_gene_res_combined) <- c("gene", "pip_combined")
    
    #number of genes with combined PIP > threshold
    n_ctwas_genes_combined <- sum(ctwas_gene_res_combined$pip_combined > PIP_threshold)
    
    #number of cTWAS genes using combined PIP that are causal
    n_causal_detected_combined <- sum(ctwas_gene_res_combined$gene[ctwas_gene_res_combined$pip_combined > PIP_threshold] %in% true_genes_combined)
    
    #collect number of SNPs analyzed by cTWAS
    ctwas_res_s1 <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.s1.susieIrss.txt"))
    n_snps <- sum(ctwas_res_s1$type=="SNP")/thin
    rm(ctwas_res_s1)
    
    
    #load estimated parameters
    load(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.s2.susieIrssres.Rd"))
    
    #estimated group prior (all iterations)
    estimated_group_prior_all <- group_prior_rec
    estimated_group_prior_all["SNP",] <- estimated_group_prior_all["SNP",]*thin #adjust parameter to account for thin argument
    
    #estimated group prior variance (all iterations)
    estimated_group_prior_var_all <- group_prior_var_rec
    
    #set group size
    group_size <- c(table(ctwas_gene_res$type), structure(n_snps, names="SNP"))
    group_size <- group_size[rownames(estimated_group_prior_all)]
    
    #estimated group PVE (all iterations)
    estimated_group_pve_all <- estimated_group_prior_var_all*estimated_group_prior_all*group_size/sample_size #check PVE calculation
    
    #multitissue TWAS analysis with bonferroni adjusted threshold for z scores
    load(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR_z_gene.Rd"))
    
    
    alpha <- 0.05
    sig_thresh <- qnorm(1-(alpha/nrow(z_gene)/2), lower=T)
    twas_genes <- z_gene$id[abs(z_gene$z)>sig_thresh]
    twas_genes_combined <- unique(sapply(twas_genes, function(x){unlist(strsplit(x, "[|]"))[1]}))
    
    n_twas_genes <- length(twas_genes)
    n_twas_genes_combined <- length(twas_genes_combined)
    
    n_twas_genes_in_causal <- sum(twas_genes %in% true_genes)
    n_twas_genes_in_causal_combined <- sum(twas_genes_combined %in% true_genes_combined)
    
    results_current <- data.frame(simutag=as.character(simutag),
                                  n_causal=as.integer(n_causal),
                                  n_causal_combined=as.integer(n_causal_combined),
                                  n_detected_pip=as.integer(n_ctwas_genes),
                                  n_detected_pip_in_causal=as.integer(n_causal_detected),
                                  n_detected_comb_pip=as.integer(n_ctwas_genes_combined),
                                  n_detected_comb_pip_in_causal=as.integer(n_causal_detected_combined),
                                  n_detected_twas=as.integer(n_twas_genes),
                                  n_detected_twas_in_causal=as.integer(n_twas_genes_in_causal),
                                  n_detected_comb_twas=as.integer(n_twas_genes_combined),
                                  n_detected_comb_twas_in_causal=as.integer(n_twas_genes_in_causal_combined))
    
    results_current[,paste0("pve_",names(estimated_group_pve_all[,30]))] <- estimated_group_pve_all[,30]
    results_current[,paste0("prior_",names(estimated_group_prior_all[,30]))] <- estimated_group_prior_all[,30]
    results_current[,paste0("prior_var_",names(estimated_group_prior_var_all[,30]))] <- estimated_group_prior_var_all[,30]
    results_df <- rbind(results_df, results_current)
  }
  return(list(results_df,group_size))
}

get_sim_joint_res_corr <- function(results_dir,runtag,configtag,simutags,thin,sample_size,PIP_threshold){
  results_df <- data.frame(simutag=as.character(),
                           prior_weight1=as.numeric(),
                           prior_weight2=as.numeric(),
                           prior_weight3=as.numeric(),
                           prior_var_snp=as.numeric(),
                           prior_var_weight1=as.numeric(),
                           prior_var_weight2=as.numeric(),
                           prior_var_weight3=as.numeric())
  
  for (i in 1:length(simutags)){
    simutag <- simutags[i]
    
    #load genes with true simulated effect
    load(paste0(results_dir, runtag, "_simu", simutag, "-pheno.Rd"))
    true_genes <- unlist(sapply(1:22, function(x){phenores$batch[[x]]$id.cgene}))
    true_genes_combined <- unique(sapply(true_genes, function(x){unlist(strsplit(x, "[|]"))[1]}))
    
    #load cTWAS results
    ctwas_res <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.drop.merge.susieIrss.txt"))
    ctwas_gene_res <- ctwas_res[ctwas_res$type!="SNP",]
    
    #number of causal genes
    n_causal <- length(true_genes)
    n_causal_combined <- length(true_genes_combined)
    
    #number of gene+tissue combinations with cTWAS PIP > threshold
    n_ctwas_genes <- sum(ctwas_gene_res$susie_pip > PIP_threshold)
    
    #number of cTWAS genes that are causal
    n_causal_detected <- sum(ctwas_gene_res$id[ctwas_gene_res$susie_pip > PIP_threshold] %in% true_genes)
    
    #collapse gene+tissues to genes and compute combined PIP
    ctwas_gene_res$gene <- sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x,"[|]"))[1]})
    ctwas_gene_res_combined <- aggregate(ctwas_gene_res$susie_pip, by=list(ctwas_gene_res$gene), FUN=sum)
    colnames(ctwas_gene_res_combined) <- c("gene", "pip_combined")
    
    #number of genes with combined PIP > threshold
    n_ctwas_genes_combined <- sum(ctwas_gene_res_combined$pip_combined > PIP_threshold)
    
    #number of cTWAS genes using combined PIP that are causal
    n_causal_detected_combined <- sum(ctwas_gene_res_combined$gene[ctwas_gene_res_combined$pip_combined > PIP_threshold] %in% true_genes_combined)
    
    #collect number of SNPs analyzed by cTWAS
    ctwas_res_s1 <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.drop.merge.s1.susieIrss.txt"))
    n_snps <- sum(ctwas_res_s1$type=="SNP")/thin
    rm(ctwas_res_s1)
    
    
    #load estimated parameters
    load(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.drop.merge.s2.susieIrssres.Rd"))
    
    #estimated group prior (all iterations)
    estimated_group_prior_all <- group_prior_rec
    estimated_group_prior_all["SNP",] <- estimated_group_prior_all["SNP",]*thin #adjust parameter to account for thin argument
    
    #estimated group prior variance (all iterations)
    estimated_group_prior_var_all <- group_prior_var_rec
    
    #set group size
    group_size <- c(table(ctwas_gene_res$type), structure(n_snps, names="SNP"))
    group_size <- group_size[rownames(estimated_group_prior_all)]
    
    #estimated group PVE (all iterations)
    estimated_group_pve_all <- estimated_group_prior_var_all*estimated_group_prior_all*group_size/sample_size #check PVE calculation
    
    #multitissue TWAS analysis with bonferroni adjusted threshold for z scores
    load(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.drop.merge_z_gene.Rd"))
    
    
    alpha <- 0.05
    sig_thresh <- qnorm(1-(alpha/nrow(z_gene)/2), lower=T)
    twas_genes <- z_gene$id[abs(z_gene$z)>sig_thresh]
    twas_genes_combined <- unique(sapply(twas_genes, function(x){unlist(strsplit(x, "[|]"))[1]}))
    
    n_twas_genes <- length(twas_genes)
    n_twas_genes_combined <- length(twas_genes_combined)
    
    n_twas_genes_in_causal <- sum(twas_genes %in% true_genes)
    n_twas_genes_in_causal_combined <- sum(twas_genes_combined %in% true_genes_combined)
    
    results_current <- data.frame(simutag=as.character(simutag),
                                  n_causal=as.integer(n_causal),
                                  n_causal_combined=as.integer(n_causal_combined),
                                  n_detected_pip=as.integer(n_ctwas_genes),
                                  n_detected_pip_in_causal=as.integer(n_causal_detected),
                                  n_detected_comb_pip=as.integer(n_ctwas_genes_combined),
                                  n_detected_comb_pip_in_causal=as.integer(n_causal_detected_combined),
                                  pve_snp=as.numeric(rev(estimated_group_pve_all["SNP",])[1]),
                                  pve_weight1=as.numeric(rev(estimated_group_pve_all["Brain_Cerebellum_harmonized",])[1]),
                                  pve_weight2=as.numeric(rev(estimated_group_pve_all["Brain_Hippocampus_harmonized",])[1]),
                                  pve_weight3=as.numeric(rev(estimated_group_pve_all["Brain_Caudate_basal_ganglia_harmonized",])[1]),
                                  prior_snp=as.numeric(rev(estimated_group_prior_all["SNP",])[1]),
                                  prior_weight1=as.numeric(rev(estimated_group_prior_all["Brain_Cerebellum_harmonized",])[1]),
                                  prior_weight2=as.numeric(rev(estimated_group_prior_all["Brain_Hippocampus_harmonized",])[1]),
                                  prior_weight3=as.numeric(rev(estimated_group_prior_all["Brain_Caudate_basal_ganglia_harmonized",])[1]),
                                  prior_var_snp=as.numeric(rev(estimated_group_prior_var_all["SNP",])[1]),
                                  prior_var_weight1=as.numeric(rev(estimated_group_prior_var_all["Brain_Cerebellum_harmonized",])[1]),
                                  prior_var_weight2=as.numeric(rev(estimated_group_prior_var_all["Brain_Hippocampus_harmonized",])[1]),
                                  prior_var_weight3=as.numeric(rev(estimated_group_prior_var_all["Brain_Caudate_basal_ganglia_harmonized",])[1]),
                                  n_detected_twas=as.integer(n_twas_genes),
                                  n_detected_twas_in_causal=as.integer(n_twas_genes_in_causal),
                                  n_detected_comb_twas=as.integer(n_twas_genes_combined),
                                  n_detected_comb_twas_in_causal=as.integer(n_twas_genes_in_causal_combined))
    
    results_df <- rbind(results_df, results_current)
  }
  return(results_df)
}

plot_par <- function(truth, est, xlabels = c("setting 1", "setting 2", "setting 3"), ...){
  colorsall <- c("#7fc97f", "#beaed4", "#fdc086")
  col = est[,1]
  est[,1] <- jitter(est[,1])
  
  plot(est, pch = 19, xaxt = "n", xlab="" ,frame.plot=FALSE, col = colorsall[col], ...)
  axis(side=1, at=1:3, labels = xlabels, tick = F)
  
  for (t in 1:nrow(truth)){
    row <- truth[t,]
    segments(row[1]-0.2, row[2] , row[1] + 0.2, row[2],
             col = colorsall[t], lty = par("lty"), lwd = 2, xpd = FALSE)
  }
  grid()
}

get_pip_distr <- function(results_dir,runtag,configtag,simutags,thin,sample_size,PIP_threshold){
  
  Liver_attr <- c()
  Lung_attr <- c()
  Brain_attr <- c()
  for (i in 1:length(simutags)){
    simutag <- simutags[i]
    
    #load cTWAS results
    ctwas_res <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.drop.merge.susieIrss.txt"))
    ctwas_gene_res <- ctwas_res[ctwas_res$type!="SNP",]
    
    #collapse gene+tissues to genes and compute combined PIP
    ctwas_gene_res$gene <- sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x,"[|]"))[1]})
    
    df_Liver <- ctwas_gene_res[ctwas_gene_res$type=="Liver_harmonized",]
    df_Lung <- ctwas_gene_res[ctwas_gene_res$type=="Lung_harmonized",]
    df_Brain <- ctwas_gene_res[ctwas_gene_res$type=="Brain_Hippocampus_harmonized",]
    
    df_gene <- aggregate(ctwas_gene_res$susie_pip, by=list(ctwas_gene_res$gene), FUN=sum)
    colnames(df_gene) <- c("gene", "combined_pip")
    df_gene$Liver_pip <-0
    df_gene$Lung_pip <-0
    df_gene$Brain_pip <-0
    
    for(i in df_gene$gene){
      if(i %in% df_Liver$gene){
        df_gene[df_gene$gene==i,"Liver_pip"] <- round(df_Liver[df_Liver$gene==i,"susie_pip"],3)
      }
    }
    
    for(i in df_gene$gene){
      if(i %in% df_Lung$gene){
        df_gene[df_gene$gene==i,"Lung_pip"] <- round(df_Lung[df_Lung$gene==i,"susie_pip"],3)
      }
    }
    
    for(i in df_gene$gene){
      if(i %in% df_Brain$gene){
        df_gene[df_gene$gene==i,"Brain_pip"] <- round(df_Brain[df_Brain$gene==i,"susie_pip"],3)
      }
    }
    
    df_gene$combined_pip <- round(df_gene$combined_pip,3)
    #sort by combined PIP
    df_gene <- df_gene[order(-df_gene$combined_pip),]
    #genes with PIP>0.8 or 20 highest PIPs
    df_gene <- df_gene[df_gene$combined_pip>0.8,]
    
    Liver_attr <- c(Liver_attr,df_gene$Liver_pip/df_gene$combined_pip)
    Lung_attr <- c(Lung_attr,df_gene$Lung_pip/df_gene$combined_pip)
    Brain_attr <- c(Brain_attr,df_gene$Brain_pip/df_gene$combined_pip)
    
  }
  return(list(Liver_attr,Lung_attr,Brain_attr))
}

get_pip_distr_corr <- function(results_dir,runtag,configtag,simutags,thin,sample_size,PIP_threshold){
  
  Liver_attr <- c()
  Lung_attr <- c()
  Brain_attr <- c()
  for (i in 1:length(simutags)){
    simutag <- simutags[i]
    
    #load cTWAS results
    ctwas_res <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.drop.merge.susieIrss.txt"))
    ctwas_gene_res <- ctwas_res[ctwas_res$type!="SNP",]
    
    #collapse gene+tissues to genes and compute combined PIP
    ctwas_gene_res$gene <- sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x,"[|]"))[1]})
    
    df_Liver <- ctwas_gene_res[ctwas_gene_res$type=="Brain_Cerebellum_harmonized",]
    df_Lung <- ctwas_gene_res[ctwas_gene_res$type=="Brain_Hippocampus_harmonized",]
    df_Brain <- ctwas_gene_res[ctwas_gene_res$type=="Brain_Caudate_basal_ganglia_harmonized",]
    
    df_gene <- aggregate(ctwas_gene_res$susie_pip, by=list(ctwas_gene_res$gene), FUN=sum)
    colnames(df_gene) <- c("gene", "combined_pip")
    df_gene$Liver_pip <-0
    df_gene$Lung_pip <-0
    df_gene$Brain_pip <-0
    
    for(i in df_gene$gene){
      if(i %in% df_Liver$gene){
        df_gene[df_gene$gene==i,"Liver_pip"] <- round(df_Liver[df_Liver$gene==i,"susie_pip"],3)
      }
    }
    
    for(i in df_gene$gene){
      if(i %in% df_Lung$gene){
        df_gene[df_gene$gene==i,"Lung_pip"] <- round(df_Lung[df_Lung$gene==i,"susie_pip"],3)
      }
    }
    
    for(i in df_gene$gene){
      if(i %in% df_Brain$gene){
        df_gene[df_gene$gene==i,"Brain_pip"] <- round(df_Brain[df_Brain$gene==i,"susie_pip"],3)
      }
    }
    
    df_gene$combined_pip <- round(df_gene$combined_pip,3)
    #sort by combined PIP
    df_gene <- df_gene[order(-df_gene$combined_pip),]
    #genes with PIP>0.8 or 20 highest PIPs
    df_gene <- df_gene[df_gene$combined_pip>0.8,]
    
    Liver_attr <- c(Liver_attr,df_gene$Liver_pip/df_gene$combined_pip)
    Lung_attr <- c(Lung_attr,df_gene$Lung_pip/df_gene$combined_pip)
    Brain_attr <- c(Brain_attr,df_gene$Brain_pip/df_gene$combined_pip)
    
  }
  return(list(Liver_attr,Lung_attr,Brain_attr))
}

plot_par <- function(truth, est, xlabels = c("setting 1", "setting 2", "setting 3"), ...){
  colorsall <- c("#7fc97f", "#beaed4", "#fdc086")
  col = est[,1]
  est[,1] <- jitter(est[,1])
  
  plot(est, pch = 19, xaxt = "n", xlab="" ,frame.plot=FALSE, col = colorsall[col], ...)
  axis(side=1, at=1:3, labels = xlabels, tick = F)
  
  for (t in 1:nrow(truth)){
    row <- truth[t,]
    segments(row[1]-0.2, row[2] , row[1] + 0.2, row[2],
             col = colorsall[t], lty = par("lty"), lwd = 2, xpd = FALSE)
  }
  grid()
}

plot_par_7 <- function(truth, est, xlabels = c("setting 1", "setting 2", "setting 3", "setting 4", "setting 5", "setting 6", "setting 7"), ...){
  colorsall <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#f781bf", "#a65628")
  col = est[,1]
  est[,1] <- jitter(est[,1])
  
  plot(est, pch = 19, xaxt = "n", xlab="" ,frame.plot=FALSE, col = colorsall[col], ...)
  axis(side=1, at=1:7, labels = xlabels, tick = F)
  
  for (t in 1:nrow(truth)){
    row <- truth[t,]
    segments(row[1]-0.2, row[2] , row[1] + 0.2, row[2],
             col = colorsall[t], lty = par("lty"), lwd = 2, xpd = FALSE)
  }
  grid()
}

plot_par_6 <- function(truth, est, xlabels = c("setting 1", "setting 2", "setting 3", "setting 4", "setting 5", "setting 6"), ...){
  colorsall <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#f781bf")
  col = est[,1]
  est[,1] <- jitter(est[,1])
  
  plot(est, pch = 19, xaxt = "n", xlab="" ,frame.plot=FALSE, col = colorsall[col], ...)
  axis(side=1, at=1:6, labels = xlabels, tick = F)
  
  for (t in 1:nrow(truth)){
    row <- truth[t,]
    segments(row[1]-0.2, row[2] , row[1] + 0.2, row[2],
             col = colorsall[t], lty = par("lty"), lwd = 2, xpd = FALSE)
  }
  grid()
}

plot_par_4 <- function(truth, est, xlabels = c("setting 1", "setting 2", "setting 3", "setting 4"), ...){
  colorsall <- c("#7fc97f", "#beaed4", "#fdc086", "#92c5de")
  col = est[,1]
  est[,1] <- jitter(est[,1])
  
  plot(est, pch = 19, xaxt = "n", xlab="" ,frame.plot=FALSE, col = colorsall[col], ...)
  axis(side=1, at=1:4, labels = xlabels, tick = F)
  
  for (t in 1:nrow(truth)){
    row <- truth[t,]
    segments(row[1]-0.2, row[2] , row[1] + 0.2, row[2],
             col = colorsall[t], lty = par("lty"), lwd = 2, xpd = FALSE)
  }
  grid()
}

plot_par_2 <- function(truth, est, xlabels = c("setting 1", "setting 2"), ...){
  colorsall <- c("#7fc97f", "#beaed4")
  col = est[,1]
  est[,1] <- jitter(est[,1])
  
  plot(est, pch = 19, xaxt = "n", xlab="" ,frame.plot=FALSE, col = colorsall[col], ...)
  axis(side=1, at=1:2, labels = xlabels, tick = F)
  
  for (t in 1:nrow(truth)){
    row <- truth[t,]
    segments(row[1]-0.2, row[2] , row[1] + 0.2, row[2],
             col = colorsall[t], lty = par("lty"), lwd = 2, xpd = FALSE)
  }
  grid()
}

get_pip_distr <- function(results_dir,runtag,configtag,simutags,thin,sample_size,PIP_threshold){
  results_df <- data.frame()
  for (i in 1:length(simutags)){
    simutag <- simutags[i]
    #load cTWAS results
    ctwas_res <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.susieIrss.txt"))
    ctwas_gene_res <- ctwas_res[ctwas_res$type!="SNP",]
    #collapse gene+tissues to genes and compute combined PIP
    ctwas_gene_res$gene <- sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x,"[|]"))[1]})
    ctwas_gene_res <- ctwas_gene_res[,c("gene","type","susie_pip")]
    # Step 1: Aggregate values by gene
    total_pips <- ctwas_gene_res %>%
      group_by(gene) %>%
      summarise(total_pips = sum(susie_pip))
    
    # Step 2: Spread the values of each tissue into separate columns for each gene
    prop_df <- ctwas_gene_res %>%
      pivot_wider(names_from = type, values_from = susie_pip, values_fill = list(value = 0)) %>%
      left_join(total_pips, by = "gene")
    
    prop_df <- prop_df[prop_df$total_pips>0.8,]
    prop_df[is.na(prop_df)] <- 0
    results_df <- rbind(results_df,prop_df)
  }
  return(results_df)
}

get_pip_distr_corr <- function(results_dir,runtag,configtag,simutags,thin,sample_size,PIP_threshold){
  
  Liver_attr <- c()
  Lung_attr <- c()
  Brain_attr <- c()
  for (i in 1:length(simutags)){
    simutag <- simutags[i]
    
    #load cTWAS results
    ctwas_res <- data.table::fread(paste0(results_dir, runtag, "_simu", simutag, "_config", configtag, "_LDR.drop.merge.susieIrss.txt"))
    ctwas_gene_res <- ctwas_res[ctwas_res$type!="SNP",]
    
    #collapse gene+tissues to genes and compute combined PIP
    ctwas_gene_res$gene <- sapply(ctwas_gene_res$id, function(x){unlist(strsplit(x,"[|]"))[1]})
    
    df_Liver <- ctwas_gene_res[ctwas_gene_res$type=="Brain_Cerebellum_harmonized",]
    df_Lung <- ctwas_gene_res[ctwas_gene_res$type=="Brain_Hippocampus_harmonized",]
    df_Brain <- ctwas_gene_res[ctwas_gene_res$type=="Brain_Caudate_basal_ganglia_harmonized",]
    
    df_gene <- aggregate(ctwas_gene_res$susie_pip, by=list(ctwas_gene_res$gene), FUN=sum)
    colnames(df_gene) <- c("gene", "combined_pip")
    df_gene$Liver_pip <-0
    df_gene$Lung_pip <-0
    df_gene$Brain_pip <-0
    
    for(i in df_gene$gene){
      if(i %in% df_Liver$gene){
        df_gene[df_gene$gene==i,"Liver_pip"] <- round(df_Liver[df_Liver$gene==i,"susie_pip"],3)
      }
    }
    
    for(i in df_gene$gene){
      if(i %in% df_Lung$gene){
        df_gene[df_gene$gene==i,"Lung_pip"] <- round(df_Lung[df_Lung$gene==i,"susie_pip"],3)
      }
    }
    
    for(i in df_gene$gene){
      if(i %in% df_Brain$gene){
        df_gene[df_gene$gene==i,"Brain_pip"] <- round(df_Brain[df_Brain$gene==i,"susie_pip"],3)
      }
    }
    
    df_gene$combined_pip <- round(df_gene$combined_pip,3)
    #sort by combined PIP
    df_gene <- df_gene[order(-df_gene$combined_pip),]
    #genes with PIP>0.8 or 20 highest PIPs
    df_gene <- df_gene[df_gene$combined_pip>0.8,]
    
    Liver_attr <- c(Liver_attr,df_gene$Liver_pip/df_gene$combined_pip)
    Lung_attr <- c(Lung_attr,df_gene$Lung_pip/df_gene$combined_pip)
    Brain_attr <- c(Brain_attr,df_gene$Brain_pip/df_gene$combined_pip)
    
  }
  return(list(Liver_attr,Lung_attr,Brain_attr))
}