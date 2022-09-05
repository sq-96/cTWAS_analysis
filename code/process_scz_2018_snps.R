library(bigsnpr)
library(data.table)
scz_2018 <- fread("/project2/compbio/gwas_summary_statistics/scz_2018/clozuk_pgc2.meta.sumstats.txt.gz")
scz_2018 <- scz_2018[,c("CHR","BP","A1","A2","OR","SE")]
scz_2018$OR <- log(scz_2018$OR)
colnames(scz_2018) <- c("chr","pos","a0","a1","beta","se")

ukbb_snps <- fread("/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/UKBB_SNPs_Info.text")
ukbb_snps <- ukbb_snps[,c("chrom","pos","alt","ref","id")]
colnames(ukbb_snps) <- c("chr","pos","a0","a1","id")

snp_matched <- snp_match(scz_2018,ukbb_snps,strand_flip = TRUE,join_by_pos = TRUE,
                         remove_dups = TRUE, match.min.prop = 0.2, return_flip_and_rev = FALSE)


snp_matched <- snp_matched[,c("id","a0","a1","beta","se")]
colnames(snp_matched) <- c("SNP", "A1", "A2", "logOR", "SE")
saveRDS(snp_matched,file = "/project2/xinhe/shengqian/cTWAS/cTWAS_analysis/data/scz_2018.RDS")


