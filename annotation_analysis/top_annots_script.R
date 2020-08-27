suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))


regressIt <- function(x, y)
{
    c <- round(runif(1,0,1), digits = 2)
   #Sample a random number each time. 2 in every 100 runs should print a qqplot.  
    lmreg <- lm(y ~ x + 1) #Does this even add the intercept?
    if(c < 0.05)
    {
        
        png(paste0(c,"qq.png"))
        plot(lmreg, which = 2)
        dev.off()
    }
    return(lmreg)
}


parser <- ArgumentParser$new()
parser$add_description("Regress annotations against snp effect sizes.")
parser$add_argument("--sum_stats", type = 'character', help = "Path to the desired summary statistics file; set up for default Neale lab UKBB format, with columns variant, minor_AF, beta, pval, and tstat")
parser$add_argument("--annots", type = "character", help = "File containing genomtype snp annotations. Be sure this has a column labelled 'SNP' with variant IDs of the same format as those in sum_stats, as well as Chr, SNP, Pos columns")
parser$add_argument("--topn", type = "numeric", help = "Specify the top number of annotations to track", default = 25)
parser$add_argument("--iterate", type = "logical", help = "Specify if you want to regress against p-value subsets (i.e. GWAS s.s. at p < 5e-8, 0.01,...)-the  iterated version", default = F, action = "store_true")
parser$add_argument("--output", type = "character", help = "Specify the output location")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')

parser$helpme()
args <- parser$get_args()

ss_in = args$sum_stats
ldsc_annots_in = args$annots
print(ldsc_annots_in)
topn <- args$topn

print("Reading in SNPs from annotations data.")
ldsc_annots <- fread(ldsc_annots_in)
print(paste0("Total of ", nrow(ldsc_annots), " in annotation data"))
print("Reading in SNPs from summary stats.")
summary_stats <- fread(args$sum_stats)
print(paste0("Total of ", nrow(summary_stats), " read in from summary stats."))

effect_sizes <- summary_stats %>% filter(variant %in% ldsc_annots$SNP) %>% select(variant, minor_AF, beta, pval, tstat) %>% rename(SNP = variant) %>% mutate(abs_tstat  = abs(tstat))



dat <- inner_join(ldsc_annots, effect_sizes) %>% drop_na()
annots_numeric <- dat %>% select(-Chr, -SNP, -Pos, -beta, -minor_AF, -pval, -tstat, -abs_tstat) %>% select(which(!colSums(., na.rm = T) == 0))
annots_info <- dat %>% select(Chr, SNP, Pos, beta, minor_AF, pval,abs_tstat,tstat)



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
  #lmreg <- lm(info$abs_tstat ~ x + 1) #Does this even add the intercept?
  #all_regressions <- data.table(Fits = apply(num,2,function(x) summary(lm(info$abs_tstat ~ x + 1)))) #updated on stuff.
  all_regressions <- data.table(Fits = apply(num,2,function(x) summary(regressIt(x,info$abs_tstat)))) #updated on stuff.
  rsquared <-unlist(lapply(all_regressions$Fits, function(x) x$r.squared))
  annot_beta <-unlist(lapply(all_regressions$Fits, function(x) x$coefficients[2]))
  annot_pval <- unlist(lapply(all_regressions$Fits, function(x) x$coefficients[8]))
  
    annot_t <- unlist(lapply(all_regressions$Fits, function(x) x$coefficients[6]))
    annot_error <- unlist(lapply(all_regressions$Fits, function(x) x$coefficients[4]))
  annot_reg_results <- data.frame(annot = colnames(annots_numeric), R2 = rsquared, beta = annot_beta, t_stat = annot_t, pval = annot_pval, std_error = annot_error, abs_beta = abs(annot_beta),abs_t = abs(annot_t), bonf_p = p.adjust(annot_pval, method = "bonferroni")) %>% arrange(-abs_t)

#  annot_reg_results <- data.frame(annot = colnames(annots_numeric), R2 = rsquared, beta = annot_beta, t_stat = annot_t, pval = annot_pval, std_error = annot_error, abs_beta = abs(annot_beta), bonf_p = p.adjust(annot_pval, method = "bonferroni")) %>% arrange(-abs_beta)
  
  #plot
  #top_annots <- annot_reg_results[1:topn,]
  #plot <- ggplot(top_annots, aes(x = annot, y = beta)) + geom_pointrange(aes(ymin = beta-std_error, ymax = beta+std_error))
  #maxp = top_annots[top_annots$abs_bet == max(top_annots$abs_beta),]$bonf_p
  #maxb = top_annots[top_annots$abs_bet == max(top_annots$abs_beta),]$beta
  #plot + scale_x_discrete(labels = abbreviate(top_annots$annot, minlength = 15)) + theme(axis.text.x = element_text(angle = 90)) + annotate(geom = "text", x =13, y =  maxb + 0.00002, label =paste("Top p=",signif(maxp,digits = 3)))
  #ggsave(width = 5, height = 7, dpi = 300, filename = paste0("top_annots_", as.character(p),".png"))
  
  #Track the ranks
  rtemp <- annot_reg_results %>% mutate(ordering = row_number()) %>% arrange(annot)
  #ranks <- ranks %>% mutate(p = rtemp$ordering)
  ranks[,as.character(p)] <-rtemp$ordering
  toppers <- rtemp %>% mutate(topTier = ifelse(ordering <= topn,1,0))
  #in_top <- in_top %>% mutate(p = toppers$topTier)
  in_top[,as.character(p)] <- toppers$topTier
print(args$output) 
write.table(annot_reg_results, file = paste0(args$output,".",p, ".effect_sizes.tsv"), quote=FALSE, sep='\t', row.names = F)

}



ranks <- ranks %>% mutate(total = rowSums(select_if(ranks, is.numeric)))
in_top <- in_top %>% mutate(total = rowSums(select_if(in_top, is.numeric)))
if(args$iterate)
{write.table(in_top, file=paste0(args$output,'top_appearances.tsv'), quote=FALSE, sep='\t', row.names = F)
write.table(ranks, file=paste0(args$output, 'overall_ranks.tsv'), quote=FALSE, sep='\t', row.names = F)}
