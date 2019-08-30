suppressMessages(library(stringr))
suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(rms))
suppressMessages(library(ggplot2))

parser <- ArgumentParser$new()
parser$add_description("R script for generating figures associated with PRS results and their phenotypes. Currently (8/26) creates a histogram and a scatter plot, with points labelled by case/control disease condition.")
parser$add_argument("--prs_results", type = 'character', help = "Prefix to the files with PRS scores, listed in a .tsv by ID and score", default = "/work-zfs/abattle4/ashton/prs_dev/scratch/calc_prs_testing/beta_results_test.tsv")
parser$add_argument("--pheno", type = "character", help = "Path to phenotype file, as given by Marios' analysis.", default = "/work-zfs/abattle4/ashton/prs_dev/CHS_race1_toy/pheno_mendrand.txt")
parser$add_argument("--output", type = "character", help = "Output file path.")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()


'%ni%' <- Negate('%in%')

#get the files we are looking at
fl <- Sys.glob(paste0(args$prs_results,"*.tsv"))
nagel_r <-c()
pval_names <- c()
phenos <- read_tsv(args$pheno) %>% rename(ID = IID)

for (f in fl)
{
    filename <- basename(f)
    prs <- read_tsv(f)
    MI_pred <- inner_join(prs, phenos) %>% select(ID, Score, MI)
    MI_pred$MI <- as.factor(MI_pred$MI)
    #ggplot(MI_pred,aes(x=Score, color = MI)) + geom_histogram() +labs(x="PRS", y="Count", title="Distribution of PRS Scores")
    #ggsave(filename = paste0(args$output,substr(filename,1,nchar(f)-4), ".histogram.png"))
    #ggplot(MI_pred, aes(x=c(0), y = Score, color = MI)) + geom_jitter(position=position_jitter(0.1)) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    #ggsave(filename = paste0(args$output,substr(filename, 1, nchar(f)-4), ".dotplot.png"))
    
    #Calculate Nagelkerke R
    lm <- lrm(MI_pred$MI ~ MI_pred$Score)
    nagel_r <- c(nagel_r, lm$stats[10])
    pval_names <- c(pval_names, str_extract(f, "[01]\\.*[\\-e\\d]+"))
    print(paste("Completed review for", f)) 
}
#Do that for each individually. Now get the
dat <- data.frame(pval_names, nagel_r)
dat$pval_names <- as.numeric(as.character(dat$pval_names))
dat <- dat[order(dat$pval_names),]
dat$pval_names <- as.factor(dat$pval_names)
print(dat)
ggplot(dat, aes(pval_names, nagel_r)) + geom_bar(stat = "identity") + labs(x="P-value threshold", y = "R2 (Nagelkerke)")
ggsave(filename = paste0(args$ouput, "bar_plot.png"))
