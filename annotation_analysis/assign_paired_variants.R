#!/usr/bin/env Rscript
#Quick script to randomly pair snps by pvalue from 2 separate lists
#arguments are: 1 and 2 are the lists (lower, upper), 3 is the summary stats, 4 is the output handle
args = commandArgs(trailingOnly=TRUE)
lower = args[1]
upper = args[2]
output = args[4]
sum_stats = args[3]
library(readr)
library(dplyr)
#create 100 bins
seqn <- seq(0,100) * 0.01
all_high <- scan(upper, what = "character()")
all_low <- scan(lower, what = "character()")
sum_stats <- read_tsv(sum_stats, col_names = F)
colnames(sum_stats) <- c("ID", "REF", "ALT", "Beta", "pval")
gwas_low_snps <- filter(sum_stats, ID %in% all_low)
gwas_high_snps <- filter(sum_stats, ID %in% all_high)
paired_high <- gwas_high_snps[0,]
paired_low <- gwas_low_snps[0,]
for (n in seqn){
  count_low_snp <- dim(gwas_low_snps %>% filter(pval >= n, pval < n + 0.01))[1]
  count_high_snp <- dim(gwas_high_snps %>% filter(pval >= n, pval < n + 0.01))[1]
  min_count <- min(count_high_snp,count_low_snp)
  t_h <- gwas_high_snps %>% filter(pval >= n, pval < n + 0.01) %>% sample_n(., min_count)
  t_l <- gwas_low_snps %>% filter(pval >= n, pval < n + 0.01) %>% sample_n(., min_count)
  paired_high <- rbind(paired_high,t_h)
  paired_low <- rbind(paired_low, t_l)

}
write_delim(paired_high %>% select(ID), paste0(output, ".paired_high_annots.tsv"), col_names = F, delim = '\t')
write_delim(paired_low %>% select(ID), paste0(output,".paired_low_annots.tsv"), col_names = F, delim = '\t')


