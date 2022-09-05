z_snp <- readRDS("/project2/xinhe/shared_data/GWAS_summary_statistics/aFib/ebi-a-GCST006414_aFib.df.rds")
z_snp$index <- paste(as.character(z_snp$chr),":",as.character(z_snp$pos),sep='')
UKBB_snp <- fread('./data/UKBB_SNPs_Info.text')
UKBB_snp$index <- paste(as.character(UKBB_snp$chrom),":",as.character(UKBB_snp$pos),sep='')



z_snp <- z_snp[,c('index','zscore','snp','a0','a1')]
UKBB_snp <- UKBB_snp[,c('index','id')]
total <- merge(z_snp,UKBB_snp,by="index")
total <- total[total$snp==total$id,]

