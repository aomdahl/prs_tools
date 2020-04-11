suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))

parser <- ArgumentParser$new()
parser$add_description("Regress annotations against snp effect sizes.")
parser$add_argument("--sum_stats", type = 'character', help = "Path to the desired summary statistics file")
parser$add_argument("--annots", type = "character", help = "File containing genomtype snp annotations")
parser$add_argument("--topn", type = "numeric", help = "Specify the top number of annotations to track", default = 25)
parser$add_argument("--iterate", type = "logical", help = "Specify if you want the iterated version", default = F, action = "store_true")
parser$add_argument("--output", type = "character", help = "Specify the output location")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')

parser$helpme()
args <- parser$get_args()

ss_in = args$sum_stats
ldsc_annots_in = args$annots
print(ldsc_annots_in)
topn <- args$topn

ldsc_annots <- fread(ldsc_annots_in) %>% drop_na()
print(paste("Reading in", dim(ldsc_annots)[0], "SNPs from annotations data."))
summary_stats <- fread(args$sum_stats) %>% drop_na()
print(paste("Reading in", dim(summary_stats)[0], "SNPs from summary stats."))


topn <- args$topn
#ldsc_annots <- fread("~/Documents/JHU/Research/LocalData/annotation_scratch/second_attempt_pruned/ldsc_db_combined.tsv")
#ldsc_annots <- ldsc_annots %>% drop_na()
#summary_stats <- fread("~/Documents/JHU/Research/LocalData/annotation_scratch/second_attempt_pruned/neal_irnt_dbp.tsv")



effect_sizes <- summary_stats %>% filter(variant %in% ldsc_annots$SNP) %>% select(variant, minor_AF, beta, pval, tstat) %>% rename(SNP = variant) %>% mutate(abs_tstat  = abs(tstat))



dat <- inner_join(ldsc_annots, effect_sizes)
annots_numeric <- dat %>% select(-CHR, -SNP, -BP, -beta, -minor_AF, -pval, -tstat, -abs_tstat) %>% select(which(!colSums(., na.rm = T) == 0))
annots_info <- dat %>% select(CHR, SNP, BP, beta, minor_AF, pval,abs_tstat,tstat)



annots_numeric_scaled <- scale(annots_numeric)
print("proceeding with regression")
#Do this at many different p-value thresholds.
pvals <- c(1)
if(args$iterate) {
pvals <- c(5e-8, 1e-6,1e-4,0.001,0.01,0.05,0.1,0.3,0.5,1)
}
ranks <- data.frame(annot = colnames(annots_numeric)) %>% arrange(annot)
in_top <- data.frame(annot = (colnames(annots_numeric))) %>% arrange(annot)
for(p in pvals)
{
  num <- annots_numeric_scaled[annots_info$pval <= p,]
  info <- annots_info[annots_info$pval <= p,]
  #all_regressions <- data.table(Fits = apply(num,2,function(x) summary(lm(info$beta ~ x + 1))))
  all_regressions <- data.table(Fits = apply(num,2,function(x) summary(lm(info$abs_tstat ~ x + 1)))) #updated on stuff.
  rsquared <-unlist(lapply(all_regressions$Fits, function(x) x$r.squared))
  annot_beta <-unlist(lapply(all_regressions$Fits, function(x) x$coefficients[2]))
  annot_pval <- unlist(lapply(all_regressions$Fits, function(x) x$coefficients[8]))
  
    annot_t <- unlist(lapply(all_regressions$Fits, function(x) x$coefficients[6]))
    annot_error <- unlist(lapply(all_regressions$Fits, function(x) x$coefficients[4]))
  annot_reg_results <- data.frame(annot = colnames(annots_numeric), R2 = rsquared, beta = annot_beta, t_stat = annot_t, pval = annot_pval, std_error = annot_error, abs_beta = abs(annot_beta),abs_t = abs(annot_t), bonf_p = p.adjust(annot_pval, method = "bonferroni")) %>% arrange(-abs_t)

#  annot_reg_results <- data.frame(annot = colnames(annots_numeric), R2 = rsquared, beta = annot_beta, t_stat = annot_t, pval = annot_pval, std_error = annot_error, abs_beta = abs(annot_beta), bonf_p = p.adjust(annot_pval, method = "bonferroni")) %>% arrange(-abs_beta)
  
  #plot
  top_annots <- annot_reg_results[1:topn,]
  plot <- ggplot(top_annots, aes(x = annot, y = beta)) + geom_pointrange(aes(ymin = beta-std_error, ymax = beta+std_error))
  maxp = top_annots[top_annots$abs_bet == max(top_annots$abs_beta),]$bonf_p
  maxb = top_annots[top_annots$abs_bet == max(top_annots$abs_beta),]$beta
  plot + scale_x_discrete(labels = abbreviate(top_annots$annot, minlength = 15)) + theme(axis.text.x = element_text(angle = 90)) + annotate(geom = "text", x =13, y =  maxb + 0.00002, label =paste("Top p=",signif(maxp,digits = 3)))
  #ggsave(width = 5, height = 7, dpi = 300, filename = paste0("top_annots_", as.character(p),".png"))
  
  #Track the ranks
  rtemp <- annot_reg_results %>% mutate(ordering = row_number()) %>% arrange(annot)
  #ranks <- ranks %>% mutate(p = rtemp$ordering)
  ranks[,as.character(p)] <-rtemp$ordering
  toppers <- rtemp %>% mutate(topTier = ifelse(ordering <= topn,1,0))
  #in_top <- in_top %>% mutate(p = toppers$topTier)
  in_top[,as.character(p)] <- toppers$topTier
  
write.table(annot_reg_results, file = paste0(args$output,p, "_effect_sizes.tsv"), quote=FALSE, sep='\t', row.names = F)

}



ranks <- ranks %>% mutate(total = rowSums(select_if(ranks, is.numeric)))
in_top <- in_top %>% mutate(total = rowSums(select_if(in_top, is.numeric)))
if(args$iterate)
{write.table(in_top, file=paste0(args$output,'top_appearances.tsv'), quote=FALSE, sep='\t', row.names = F)}
write.table(ranks, file=paste0(args$output, 'overall_ranks.tsv'), quote=FALSE, sep='\t', row.names = F)
