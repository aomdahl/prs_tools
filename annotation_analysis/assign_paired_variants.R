#!/usr/bin/env Rscript
#Quick script to randomly pair snps by pvalue from 2 separate lists
#arguments are: 1 and 2 are the lists (lower, upper), 3 is the summary stats, 4 is the output handle
suppressMessages(library(Xmisc))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
library(microbenchmark)


#args = commandArgs(trailingOnly=TRUE)
#lower = args[1]
#upper = args[2]
#output = args[4]
#sum_stats = args[3]

parser <- ArgumentParser$new()
parser$add_description("Pair snps by p-value based on GWAS")
parser$add_argument("--debug", action = "store_true", help = "If you want debug mode.", type = 'logical', default = FALSE)
parser$add_argument("--lower", type = "character", help = "File containing depleted variants.")
parser$add_argument("--upper", type = 'character', help = "File containing enriched variants")
parser$add_argument("--output", type = "character", help = "Output file path and name. end in .png")
parser$add_argument("--sum_stats", type = "character", help = "Specify path to summary stats with the values for pairing")
parser$add_argument("--null", type = "character", help = "List containing a mixture of unselected samples/null distribution, for comparison. Optional.", default = "")
parser$add_argument("--pair_metric", type = "character", help = "Specify what to split on;give a column name", default = "MAF")
parser$add_argument("--null_rep", type = "numeric", help = "Number of null subsets to create. requires --null.", default = 1)
parser$add_argument("--background", type = "character", help = "Background list to pair with. Goes with test-sample.", default = "")
parser$add_argument("--test_sample", type = "character", help = "Sample list to pair with. Goes with background", default = "")
parser$add_argument("--lock_min", type = "logical", action = "store_true", help = "If you want to ensure the # of SNPs you pick is always based on the upper subset", default = FALSE)
parser$add_argument("--legit", type = "logical", action = "store_true", help = "The hard core way", default = FALSE)
parser$add_argument("--range", type = "numeric", help = "How big to make window size.", default = 0.0001)
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
if(args$background != "" & args$test_sample != "")
{
    lower <- args$background
    upper <- args$test_sample
}else
{
    lower <- args$lower
    upper <- args$upper
}
output <- args$output
sum_stats <- args$sum_stats

all_high <- scan(upper, what = "character()")

all_low <- scan(lower, what = "character()")
sum_stats <- fread(sum_stats) %>% select(variant, minor_AF, pval) %>% drop_na()
colnames(sum_stats) <- c("ID", "MAF", "pval")
gwas_low_snps <- filter(sum_stats, ID %in% all_low)
gwas_high_snps <- filter(sum_stats, ID %in% all_high)
paired_high <- gwas_high_snps[0,]
paired_low <- gwas_low_snps[0,]
print("Splitting beginning....")
split_on <- args$pair_metric 
print(split_on)

if(args$background != "" & args$test_sample != "")
{
    getClosest <- function(val, col)
    { which(abs(val-col)==min(abs(val-col))) } #gives index of which one matches
    
    
    print("You've opted for the simple match option")
    #Note that if you just match directly, it will give you many of the same exact snps.
    #To increase variety, we will round
    t <- round(gwas_high_snps$MAF, digits = 4)
    test_samp <- gwas_high_snps %>% select(-MAF) %>% mutate("MAF" = t)
     #%>% transmute("MAF" = round(as.numeric(MAF, digits = 5)))
     t <- round(gwas_low_snps$MAF, digits = 4)
     rand_list <- sample(seq(1, dim(gwas_low_snps)[1]), replace = F)
    back_samp <- sample(gwas_low_snps %>% select(-MAF) %>% mutate("MAF" = t), replace = F) %>% mutate("tiebreaker" = rand_list) %>% arrange(MAF, tiebreaker)
    #This step is reqlly slow.
    #match_indices <- apply(test_samp, 1, function(x) getClosest(as.numeric(x[2]), back_samp$MAF)) #gives index in back_samp that matches
    match_indices <- findInterval(sample(test_samp$MAF, replace = F), back_samp$MAF)
    dups <- duplicated(match_indices)
    keeps <- test_samp[!dups,] #this is the sampels that had a successful match
    retrys <- test_samp[dups,] #these ones had a duplicate match.
    background_keep <- back_samp[match_indices[!dups],] #These are the ones we are keeping
    background_retry <- back_samp[-(match_indices[!dups]),] #Everything NOT in this list of indices
    while(dim(retrys)[1] > 0)
    {
        #match_indices <- apply(retrys, 1, function(x) getClosest(x$MAF, background_retry$MAF))
        match_indices <- findInterval(sample(retrys$MAF, replace = F), background_retry$MAF)
        dups <- duplicated(match_indices)
        keeps <- retrys[!dups,] #samples that had a successful match
        retrys <- retrys[dups,] #these ones had a duplicate match.
        background_keep <- rbind(background_keep, background_retry[match_indices[!dups],]) #These are the ones we are keeping
        
        bst <- background_retry[-(match_indices[!dups]),]
        rand_list <- sample(seq(1, dim(bst)[1]), replace = F)
        background_retry <- sample(bst, replace = F) %>% select(-tiebreaker) %>% mutate("tiebreaker" = rand_list) %>% arrange(MAF, tiebreaker) #Everything NOT in this list of indices
    }
    #Calculate the percent overlap
    overlap = sum(background_keep$ID %in% test_samp$ID)/dim(test_samp)[1]
    print(paste("Note there is some overlap between the sets:", overlap, "%"))
    #write_delim(background_keep %>% select(ID), paste0(output, ".background_MAFpaired.tsv"), col_names = F, delim = '\t')
     write_delim(background_keep, paste0(output, ".background_MAFpaired.tsv"), col_names = F, delim = '\t')
    #write_delim(test_samp, paste0(output,".fullset_MAFpaired.tsv"), col_names = F, delim = '\t')
    quit()

}

if(args$pair_metric == "pval")
{
    seqn_n <- c(5e-8,1e-6, 1e-5, 1e-4, 1e-3, seq(1,100)*0.01)
    seqn <- c(0,5e-8, 1e-6, 1e-5, 1e-4, 1e-3,seq(1,100)*0.01)
} else
{
    seqn_n <- c(seq(1,5000) * 0.0001)
    seqn <- c(seq(0,5000)*0.0001)
}

if(args$legit)
{
    gwas_bg_snps <- as.data.table(gwas_low_snps %>% arrange(MAF))
    gwas_high_snps  <- as.data.table(gwas_high_snps %>% arrange(MAF))
    #dt_gwas_bg_snps <- as.data.table(gwas_bg_snps %>% arrange(MAF))
    #dt_gwas_high_snps  <- as.data.table(gwas_high_snps %>% arrange(MAF))

    setkey(gwas_bg_snps, MAF)
        for (i in 1:dim(gwas_high_snps)[1]){
        curr_maf <- gwas_high_snps[i,]$MAF
        range <- args$range
        #print("benchmark test")
        #n <- microbenchmark(dt_gwas_bg_snps %>% filter(MAF >= curr_maf - range, MAF < curr_maf + range) %>% filter(!(ID %in% paired_low$ID)), times = 50)
        #o <- microbenchmark(gwas_bg_snps %>% filter(MAF >= curr_maf - range, MAF < curr_maf + range) %>% filter(!(ID %in% paired_low$ID)), times = 50)
        #print("original")        
        #print(o)
        #print('new')
        #print(n)
        #readline(prompt="Press [enter] to continue") 
        #gwas_bg_samps <- gwas_bg_snps %>% filter(MAF >= curr_maf - range, MAF < curr_maf + range) %>% filter(!(ID %in% paired_low$ID))
        gwas_bg_samps <- gwas_bg_snps[(MAF >= curr_maf - range) & (MAF < curr_maf + range), ] %>% filter(!(ID %in% paired_low$ID))

        while(dim(gwas_bg_samps)[1] == 0)
        {
            range = range + args$range
            gwas_bg_samps <- gwas_bg_snps[(MAF >= curr_maf - range) & (MAF < curr_maf + range), ] %>% filter(!(ID %in% paired_low$ID))
            #gwas_bg_samps <- gwas_bg_snps %>% filter(MAF >= curr_maf - range, MAF < curr_maf + range) %>% filter(!(ID %in% paired_low$ID))
        }
        paired_low <- rbind(paired_low, sample_n(gwas_bg_samps,1))
        if(i %% 10000 == 0)
        {
            print("Passing 10000")
        }
   }
    paired_high <- gwas_high_snps 
    
} else if (args$null == "") {
    for (i in 1:length(seqn_n)){
      n = seqn[i]
      next_n = seqn_n[i]
      count_low_snp <- dim(gwas_low_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n))[1]
      count_high_snp <- dim(gwas_high_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n))[1]
      
      min_count <- min(count_high_snp,count_low_snp)
      if (args$lock_min) 
        {
            min_count <- count_high_snp
            if(min_count > count_low_snp)
            {
                print("You have more in high than in low. There may be an error")
                quit()
            }
        }
      t_h <- gwas_high_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n) %>% sample_n(., min_count)
      t_l <- gwas_low_snps %>% filter(!!sym(split_on) >= n, !!sym(split_on) < next_n) %>% sample_n(., min_count)
      paired_high <- rbind(paired_high,t_h)
      paired_low <- rbind(paired_low, t_l)

    }

} else {
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
    write_delim(paired_high, paste0(output, ".paired_high_annots.tsv"), col_names = T, delim = '\t')
    write_delim(paired_low , paste0(output,".paired_low_annots.tsv"), col_names = T, delim = '\t')
}else{
    write_delim(paired_high %>% select(ID), paste0(output, ".paired_high_annots.tsv"), col_names = F, delim = '\t')
    write_delim(paired_low %>% select(ID), paste0(output,".paired_low_annots.tsv"), col_names = F, delim = '\t')
    overlap = sum(paired_high$ID %in% paired_low$ID)/dim(paired_high)[1]
    print(paste0("Overlap: ", overlap))
}



