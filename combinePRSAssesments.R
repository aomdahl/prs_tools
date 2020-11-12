
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(cluster))
suppressMessages(library(Xmisc))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

parser <- ArgumentParser$new()
parser$add_description("Quickly combine histograms")
parser$add_argument("--label_1", type = 'character', default  = "Group 1", help = "Name for group 1")
parser$add_argument("--label_2", type = "character", help = " Name for group 2.", default = "Group 2")
parser$add_argument("--group_1", type = 'character', help = "File containing PRS data for group 1/enriched set")
parser$add_argument("--group_2", type = "character", help = " File containing PRS data for group 2/ depleted set.")
parser$add_argument("--label_mixed", type = "character", default = "Group 3", help = "Name for the mixed /3rd group")
parser$add_argument("--mixed_results", type = "character", help = " File containing PRS data for third set (optional).", default = "")
parser$add_argument("--per_snp_r2", type = "logical", action = "store_true", default = FALSE, help ="Specify if you are using files with the per-snp-r2 or not")
parser$add_argument("--group_heritability", type = "logical", action = "store_true", default = FALSE, help = "Specify if you wish to compare the relative R2 to a full group. Note that this assumes mixeed_results contains the background set to compare against")
parser$add_argument("--multi_mixed", type = "character", help = "If there are many many multi-mixed,for use as a null. Treats mixed_results as an extension.", default = "")
parser$add_argument("--output", type = "character", help = "Output file path and name. end in .png")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--annot_name', help = "Specify the name of the annotation of interest", type = "character", default = "Annotation")
parser$helpme()
args <- parser$get_args()
if(args$group_heritability & args$mixed_results == "")
{ print("For now, we assume that mixed_results is your full estimate of heritability. Please specify this arugment")
    quit()
}
if(args$group_heritability & args$multi_mixed != "")
{ print("For now, we don't handle cases with iteration and h2 calculations")
    quit()
}


high <- args$group_1
low <- args$group_2
annot <- args$annot_name
outpath <- args$output
low_corr <- read_tsv(low) %>% mutate(source= args$label_2)  %>% mutate(pval_names = as.factor(pval_names))
high_corr <- read_tsv(high) %>% mutate(source = args$label_1)%>% mutate(pval_names = as.factor(pval_names))

r2_col <- "r2"
lab_name <- expression(paste("R" ^ "2"))
if(args$per_snp_r2)
{
    r2_col <- "r2_per_snp"
    lab_name <- expression(paste("Per-snp R" ^ "2"))
}

if(args$mixed_results == "")
{
    fin <- rbind(low_corr, high_corr) #%>% mutate(bonf_p = p.adjust(pval_beta, method ="bonferroni"))
    ggplot(data = fin, aes((pval_names), !!sym(r2_col), fill = source)) + geom_bar(stat = "identity", position = "dodge") + xlab("GWAS Pvalue thresholds") + ylab(lab_name) + scale_fill_discrete(name = "Sample sets") + theme_minimal_grid(12) # labels = c(args$label_1, args$label_2))
}else{
if (args$multi_mixed != "")
{
    lf <- list.files(args$mixed_results, pattern = args$multi_mixed)
    app <- read_tsv(paste0(args$mixed_results, "/", lf[1])) %>% mutate(source= paste0("std_", "1"))  %>% mutate(pval_names = as.factor(pval_names))
    for (i in 2:length(lf))
    {
        mixed <- read_tsv(paste0(args$mixed_results, "/",lf[i])) %>% mutate(source= paste0("std_", i))  %>% mutate(pval_names = as.factor(pval_names))
        app <- rbind(app, mixed) 
    }
    #Try it with the range of scores
     #with_bars <- app %>% group_by(pval_names) %>% summarize(avg = mean(r2), std_dev = range(r2)/2)
    with_bars <- app %>% group_by(pval_names) %>% summarize(avg = mean(r2), std_dev = sd(r2))
     
    write_tsv(with_bars, paste0(output, "combined_null_var.tsv"))
    #Add error bars to the other ones,,,
    with_bars <- with_bars %>% rename("r2" = avg) %>% mutate("source" = "std")
    low_corr <- low_corr %>% mutate("std_dev" = NA) %>% select(-pval_beta)
    high_corr <- high_corr %>% mutate("std_dev" = NA) %>% select(-pval_beta)
    fin <- rbind(low_corr, high_corr, with_bars)
    ggplot(data = fin, aes((pval_names), r2, fill = source)) + geom_bar(stat = "identity", position = "dodge") + xlab("GWAS Pvalue thresholds") + ylab(expression(paste("R"^"2"))) + scale_fill_discrete(name = "Sample sets", labels = c(args$label_1, args$label_2, args$label_mixed)) + geom_errorbar( aes(x=pval_names, ymin=r2-std_dev, ymax=r2+std_dev), position = "dodge")
    

}else {
mixed <- read_tsv(args$mixed_results) %>% mutate(source= args$label_mixed)  %>% mutate(pval_names = as.factor(pval_names))
fin <- rbind(low_corr, high_corr, mixed)
ggplot(data = fin, aes((pval_names), !!sym(r2_col), fill = source)) + geom_bar(stat = "identity", position = "dodge") + xlab("GWAS Pvalue thresholds") + ylab(lab_name) + scale_fill_discrete(name = "Sample sets") + theme_minimal_grid(12)#, labels = c(args$label_1, args$label_2, args$label_mixed))
}
}
ggsave(paste0(outpath, ".png"), height = 6, width = 9)

#To calculate the trait heritability, its a simple matter of 
#(trait_r2_per_snp) * (total_r2_per_snp) ^-1
#Usually its (trait_r2/full_r2)/(trait_snps/total_snps), but when you move things around it comes to be the same thing!

if(args$group_heritability & args$per_snp_r2 & args$mixed_results != "")
{   
    high_corr <- arrange(high_corr, as.numeric(pval_names))
    low_corr <- arrange(low_corr, as.numeric(pval_names))
    mixed <- arrange(mixed, as.numeric(pval_names))
    enrich_high <- high_corr$r2_per_snp/mixed$r2_per_snp
    enrich_low <- low_corr$r2_per_snp/mixed$r2_per_snp
    herit_high <- high_corr$r2/mixed$r2
    herit_low <- low_corr$r2/mixed$r2

    ret <- data.frame("pval" = mixed$pval_names, "full_r2" = mixed$r2, "full_per_snp_r2" = mixed$r2_per_snp, "group1_prop_h2" = herit_high, "group1_enrichment" = enrich_high, "group2_prop_h2" = herit_low, "group2_enrichment" = enrich_low) %>% mutate_if(is.numeric, round, digits = 4)
    write_tsv(ret, paste0(output, "_herit-enrich.tsv"))
}   

