suppressMessages(library(stringr))
suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(rms))
suppressMessages(library(ggplot2))

source("/work-zfs/abattle4/ashton/prs_dev/prs_tools/liability_pseudoR2.R")

parser <- ArgumentParser$new()
parser$add_description("R script for generating figures associated with PRS results and their phenotypes. Currently (8/26) creates a histogram and a scatter plot, with points labelled by case/control disease condition.")
parser$add_argument("--prs_results", type = 'character', help = "Prefix to the files with PRS scores, listed in a .tsv by ID and score", default = "/work-zfs/abattle4/ashton/prs_dev/scratch/calc_prs_testing/beta_results_test.tsv")
parser$add_argument("--pheno", type = "character", help = "Path to phenotype file, as given by Marios' analysis.", default = "/work-zfs/abattle4/ashton/prs_dev/CHS_race1_toy/pheno_mendrand.txt")
parser$add_argument("--output", type = "character", help = "Output file path.")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--risk_trait', help = "Specify the name of the condition we are scoring on, a string", type = "character", default = "MI")
parser$add_argument("--case_control", type = "logical", action = "store_true", help = "Specify this if it isn't a continuous trait. Otherwise program will guess based on data", default = FALSE)
parser$add_argument("--r2", help = "Specify which type of r2 or pseudo-r2 to use. Default is either pearson or liability scale.", default = "Pearson", type = "character")



parser$helpme()
args <- parser$get_args()

#A function to plot the correlation bar plots

plotCorr <- function(dat, output, style_name)
{
    print(dat)
    print(output)
    print(style_name)
    dat[,1] <- as.numeric(as.character(dat[,1]))
    dat <- dat[order(dat$pval_names),]
    dat$pval_names <- as.factor(dat$pval_names)
    ggplot(dat, aes(pval_names, corr_r)) + geom_bar(stat = "identity") + labs(x="P-value threshold", y = paste("R2",style_name))
    ggsave(filename = paste0(output, "bar_plot.",style_name, ".png"))    
}

'%ni%' <- Negate('%in%')

#get the files we are looking at
fl <- Sys.glob(paste0(args$prs_results,"*.tsv"))
#fl <- c(args$prs_results)
print(fl)
corr_r <- c()
corr_name <- args$r2
pval_names <- c()
trait <- args$risk_trait
print(trait)
phenos <- read_tsv(args$pheno) %>% select("IID", "gender", "age", "DBP")
print(phenos)
print("Starting soon...")
for (f in fl)
{
    filename <- basename(f)
    prs <- read_tsv(f)
    print(prs)
    print(f)
    MI_pred <- inner_join(prs, phenos, by = "IID")
    pval_names <- names(prs %>% select(-IID))
    print(pval_names)
    MI_pred <- na.omit(MI_pred)
    #write_tsv(MI_pred, "sandbox_copy.tsv")
    #If the trait is continuous, do regular R2. 
    if(length(unique(MI_pred[trait])) > 2 || !(args$case_control))
    {
        for (n in pval_names)
{
        corr_r <- c(corr_r, cor(MI_pred[trait], MI_pred[n]))
        #Some other plots
        ggplot(MI_pred, aes(x=c(0), y =n)) + geom_jitter(position=position_jitter(0.1), aes(fill = MI_pred$DBP)) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        ggsave(filename = paste0(n, "dotplot.png"))
} 
   }
    else #Its a case control trait
    {
        #This step is problematic.
        print("Checkpoint")
        #MI_pred[,4] <- as.factor(MI_pred[,4])
        #ggplot(MI_pred,aes(x=Score, color = trait)) + geom_histogram() +labs(x="PRS", y="Count", title="Distribution of PRS Scores")
        #ggsave(filename = paste0(args$output,substr(filename,1,nchar(f)-4), ".histogram.png"))
        #ggplot(MI_pred, aes(x=c(0), y = Score, color = MI)) + geom_jitter(position=position_jitter(0.1)) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        #ggsave(filename = paste0(args$output,substr(filename, 1, nchar(f)-4), ".dotplot.png"))
        
        
        #otherwise
        #Calculate Nagelkerke R
        for(p in pval_names){
            if(args$r2 == "nagelkerke")
            {
                corr_r <- c(corr_r, R2ObservedNagelkerke(trait,p, MI_pred))
            }
            else
            {
                corr_r <- c(corr_r, R2LiabilityProbit(trait, p, MI_pred))
            }
        }
    }
    #pval_names <- c(pval_names, str_extract(f, "[501]\\.*[\\-e\\d]+"))
    print(paste("Completed review for", f)) 
}
#Do that for each individually. Now get the
dat <- data.frame(pval_names, corr_r)
write_tsv(dat, "corr_counts.tsv")
plotCorr(dat, args$output, corr_name)

