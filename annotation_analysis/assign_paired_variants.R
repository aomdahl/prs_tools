#!/usr/bin/env Rscript
#Quick script to randomly pair snps by pvalue from 2 separate lists
#arguments are: 1 and 2 are the lists (lower, upper), 3 is the summary stats, 4 is the output handle
suppressMessages(library(Xmisc))
library(readr)
library(dplyr)
library(data.table)
library(magrittr)


#args = commandArgs(trailingOnly=TRUE)
#lower = args[1]
#upper = args[2]
#output = args[4]
#sum_stats = args[3]

parser <- ArgumentParser$new()
parser$add_description("Pair snps by p-value based on GWAS")
parser$add_argument("--debug", action = "store_true", help = "If you want debug mode.", type = 'logical')
parser$add_argument("--lower", type = "character", help = "File containing depleted variants.")
parser$add_argument("--upper", type = 'character', help = "File containing enriched variants")
parser$add_argument("--output", type = "character", help = "Output file path and name. end in .png")
parser$add_argument("--sum_stats", type = "character", help = "Specify path to summary stats with the values for pairing")
parser$add_argument("--null", type = "character", help = "List containing a mixture of unselected samples/null distribution, for comparison. Optional.", default = "")
parser$add_argument("--pair_metric", type = "character", help = "Specify what to split on;give a column name", default = "pval")
parser$add_argument("--null_rep", type = "numeric", help = "Number of null subsets to create. requires --null.", default = 1)
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
lower <- args$lower
upper <- args$upper
output <- args$output
sum_stats <- args$sum_stats

if(args$pair_metric == "pval")
{
    seqn_n <- c(5e-8,1e-6, 1e-5, 1e-4, 1e-3, seq(1,100)*0.01)
    seqn <- c(0,5e-8, 1e-6, 1e-5, 1e-4, 1e-3,seq(1,100)*0.01)
} else
{
    seqn_n <- c(seq(1,100) * 0.01)
    seqn <- c(seq(0,100)*0.01)
}
all_high <- scan(upper, what = "character()")

all_low <- scan(lower, what = "character()")
#I mistakenly made this for my sum stats, not the GWAS ones broadly. Need to fix this...
#sum_stats <- read_tsv(sum_stats, col_names = F)
#colnames(sum_stats) <- c("ID", "REF", "ALT", "Beta", "pval")
sum_stats <- fread(sum_stats) %>% select(variant, minor_AF, pval)
colnames(sum_stats) <- c("ID", "MAF", "pval")
gwas_low_snps <- filter(sum_stats, ID %in% all_low)
gwas_high_snps <- filter(sum_stats, ID %in% all_high)
paired_high <- gwas_high_snps[0,]
paired_low <- gwas_low_snps[0,]
print("splitting beginning....")
split_on <- args$pair_metric 
print(split_on)
if(args$null == "")
{
    for (i in 1:length(seqn_n)){
      n = seqn[i]
      next_n = seqn_n[i]
      count_low_snp <- dim(gwas_low_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n))[1]
      count_high_snp <- dim(gwas_high_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n))[1]
      min_count <- min(count_high_snp,count_low_snp)
      t_h <- gwas_high_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n) %>% sample_n(., min_count)
      t_l <- gwas_low_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n) %>% sample_n(., min_count)
      paired_high <- rbind(paired_high,t_h)
      paired_low <- rbind(paired_low, t_l)

    }

} else{
all_null <- scan(args$null, what = "character")
gwas_null_snps <- filter(sum_stats, ID %in% all_null)
paired_null <- list()
for (i in 1: args$null_rep)
{
    paired_null[[i]] <-  data.frame(gwas_high_snps[0,]) 

}
for (i in 1:length(seqn_n)){
      n = seqn[i]
      next_n = seqn_n[i]
    count_low_snp <- dim(gwas_low_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n))[1]
      count_high_snp <- dim(gwas_high_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n))[1]
      count_null_snp <- dim(gwas_null_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n))[1]
      min_count <- min(count_high_snp,count_low_snp, count_null_snp)
      t_h <- gwas_high_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n) %>% sample_n(., min_count)
      t_l <- gwas_low_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n) %>% sample_n(., min_count)
      paired_high <- rbind(paired_high,t_h)
      paired_low <- rbind(paired_low, t_l)
      
    for (i in 1:args$null_rep)
    {
        t_n <- gwas_null_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n) %>% sample_n(., min_count)
        paired_null[[i]] <- rbind(paired_null[[i]], t_n)

    }
}
for (i in 1:args$null_rep)
{
    write_delim(paired_null[[i]] %>% select(ID), paste0(output, i,".paired_null_annots.tsv"), col_names = F, delim = '\t')
}
}
if(args$debug)
{
write_delim(paired_high, paste0(output, ".paired_high_annots.tsv"), col_names = F, delim = '\t')
write_delim(paired_low , paste0(output,".paired_low_annots.tsv"), col_names = F, delim = '\t')

}else{
write_delim(paired_high %>% select(ID), paste0(output, ".paired_high_annots.tsv"), col_names = F, delim = '\t')
write_delim(paired_low %>% select(ID), paste0(output,".paired_low_annots.tsv"), col_names = F, delim = '\t')
}

