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
parser$add_argument("--prs_results", type = 'character', help = "Prefix to the files with PRS scores, listed in a .tsv by ID and score")
parser$add_argument("--pheno", type = "character", help = "Path to phenotype file, as given by Marios' analysis.", default = "/work-zfs/abattle4/ashton/prs_dev/CHS_race1_toy/pheno_mendrand.txt")
parser$add_argument("--output", type = "character", help = "Output file path.")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--risk_trait', help = "Specify the name of the condition we are scoring on, a string", type = "character", default = "LDL")
parser$add_argument("--case_control", type = "logical", action = "store_true", help = "Specify this if it isn't a continuous trait. Otherwise program will guess based on data", default = FALSE)
parser$add_argument("--delim", type = "character", help = "Specify which kind of delimiter in file, could be a ' ', '\\t' or otherwise", default = "\t")
parser$add_argument("--r2", help = "Specify which type of r2 or pseudo-r2 to use. Default is either linear R2 or liability scale.", default = "linear", type = "character")
parser$add_argument("--covars", type = "character", help = "Specify column names with covariates to correct for, such as 'gender'. Default is none. Write as covar1+covar2+..", default = "")
parser$helpme()
args <- parser$get_args()

#A function to plot the correlation bar plots

plotCorr <- function(dat, output, style_name)
{
    dat[,1] <- as.numeric(as.character(dat[,1]))
    dat <- dat[order(dat$pval_names),]
    dat$pval_names <- as.factor(dat$pval_names)
    ggplot(dat, aes(pval_names, r2)) + geom_bar(stat = "identity") + labs(x="P-value threshold", y = expression(paste(R ^ 2,style_name)))
    ggsave(filename = paste0(output, "bar_plot.",style_name, ".png"))    
}
'%ni%' <- Negate('%in%')
if(length(args$prs_results) == 0)
{
    print("No prs data provided. Program will terminate.")
    quit(save = "no")
}
#get the files we are looking at

#fl <- Sys.glob(paste0(args$prs_results,"*.tsv"))
fl <- c(args$prs_results)
r2 <- c()
r2_name <- args$r2
pval_names <- c()
trait <- args$risk_trait
#NEED TO FIX THIS TO TAKE any in the arcument
select_cols <- c("IID", trait, as.list(strsplit(args$covars, '+', fixed = T)[[1]]))
phenos <- read_delim(args$pheno, delim = args$delim) %>% select(select_cols) %>% na.omit()
#phenos <- read_delim(args$pheno, delim = args$delim) %>% select("IID", "gender", "age", trait) %>% na.omit()
#phenos <- read_tsv(args$pheno) %>% select("IID", trait) %>% na.omit()
print(paste0("Phenotype data for ", nrow(phenos), " individuals"))
for (f in fl)
{
    filename <- basename(f)
    prs <- read_delim(f, delim = args$delim)
    print(paste0("PRS scores for ", nrow(prs), " individuals."))
    MI_pred <- inner_join(prs, phenos, by = "IID")
    print(paste0(nrow(MI_pred), " individuals had both PRS scores and phenotype data"))
    pval_names <- names(prs %>% select(-IID))
    MI_pred <- na.omit(MI_pred)
    #If the trait is continuous, do regular R2.
    if(nrow(unique(MI_pred[trait])) > 2 && !(args$case_control))
    {
        
        for (n in pval_names)
        {
            if(args$covars == "")
            {
            #r2 <- c(r2, cor(MI_pred[trait], MI_pred[n]))
            lmv <- lm(MI_pred[[trait]] ~ MI_pred[[n]])
            r2 <- c(r2, summary(lmv)$r.squared)
            }
            else
            {
                pname = paste0("\`",n,"\`")
                f_full <- as.formula(paste(trait, paste(pname,"+", args$covars), sep = " ~ "))
                lmv_corr <- lm(f_full, data = MI_pred)
                #lmv_corr <- lm(MI_pred[[trait]] ~ MI_pred[[n]] + MI_pred[[args$covars]])
                f_part <-  as.formula(paste(trait, args$covars,sep = " ~ "))
                lmv_t <- lm(f_part, data = MI_pred)
                #lmv_t <-  lm(MI_pred[[trait]] ~ MI_pred[[args$covars]])
                r2 <- c(r2, summary(lmv_corr)$r.squared - summary(lmv_t)$r.squared)
            }
        #Some other plots
        #ggplot(MI_pred, aes(x=c(0), y =n)) + geom_jitter(position=position_jitter(0.1), aes(fill = MI_pred[trait])) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        #ggsave(filename = paste0(args$output, n, "dotplot.png"))
        #print("Succesfully generated dotplot.")
        }        
   }
    else #Its a case control trait
    {
        #This step is problematic.
        #MI_pred[,4] <- as.factor(MI_pred[,4])
        #ggplot(MI_pred,aes(x=Score, color = trait)) + geom_histogram() +labs(x="PRS", y="Count", title="Distribution of PRS Scores")
        #ggsave(filename = paste0(args$output,substr(filename,1,nchar(f)-4), ".histogram.png"))
        #ggplot(MI_pred, aes(x=c(0), y = Score, color = MI)) + geom_jitter(position=position_jitter(0.1)) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        #ggsave(filename = paste0(args$output,substr(filename, 1, nchar(f)-4), ".dotplot.png"))
        
        
        #otherwise
        #TODO: set this up for correcting for covariates too
        #Calculate Nagelkerke R
        for(p in pval_names){
            if(args$r2 == "nagelkerke")
            {
                r2 <- c(r2, R2ObservedNagelkerke(trait,p, MI_pred))
                r2_name <- "Nagelkerke" 
            }
            else
            {
                r2 <- c(r2, R2LiabilityProbit(trait, p, MI_pred))
                r2_name <- "Liability Scale (Probit)"
            }
        }
    }
    #pval_names <- c(pval_names, str_extract(f, "[501]\\.*[\\-e\\d]+"))
    print(paste("Completed review for", f)) 
}
#Do that for each individually. Now get the
dat <- data.frame(pval_names, r2)
write_tsv(dat, paste0(args$output,"_r2counts.tsv"))
plotCorr(dat, args$output, r2_name)

