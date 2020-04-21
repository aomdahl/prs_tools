suppressMessages(library(stringr))
suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(rms))
suppressMessages(library(ggplot2))

source("./prs_tools/liability_pseudoR2.R")

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
parser$add_argument("--category_split", type = "character", help = "Specify different subgroups to score on, such as ancestry. This will score the groups separately.", default = "")
parser$add_argument("--hide_pvals", type = "logical", help = "Select this to omit pvalues on bar plots", default = FALSE, action = "store_true")
parser$helpme()
args <- parser$get_args()
options(readr.num_columns = 0)
#A function to plot the correlation bar plots

plotCorr <- function(dat, output, style_name, category_var, no_pvals)
{
    library("ggsci")
    #dat[,1] <- as.numeric(as.character(dat[,1]))
    #dat[,1] <- as.character(dat[,1])
    dat$pval_names <- as.character(dat$pval_names)
    dat <- dat[order(dat$pval_names),]
    dat$pval_names <- as.factor(dat$pval_names)
    dat$logp <- -log2(dat$pval_beta)
    if(category_var == ""){
        base <- ggplot(dat, aes(x = pval_names, y = r2, fill =logp )) + geom_bar(stat = "identity") + 
            labs(x="P-value threshold", y = expression(paste(R ^ 2))) + labs(fill = "-log2(pval)") + scale_fill_gsea()
        if(no_pvals) {base}
        else {base + geom_text(aes(label=as.character(round(pval_beta, digits =3))), position=position_dodge(width = 0.9), vjust = -0.2)}
         
    } else {
        base <- ggplot(dat, aes(x = pval_names, y = r2, fill =category )) + geom_bar(stat = "identity", position = "dodge") + 
            labs(x="P-value threshold", y = expression(paste(R ^ 2))) + scale_fill_discrete(name = category_var)
        if(no_pvals) {base}
        else {base  + geom_text(aes(label = as.character(round(pval_beta, digits = 3))), position = position_dodge(width = 0.9), vjust = -0.3, size = 3)}
    }
     
    ggsave(filename = paste0(output, "bar_plot.",style_name, ".png"), height = 7, width = 8.5) 
}

'%ni%' <- Negate('%in%')
if(length(args$prs_results) == 0)
{
    print("No prs data provided. Program will terminate.")
    quit(save = "no")
}
#get the files we are looking at

fl <- c(args$prs_results)
r2 <- c()
pval_beta <- c()
r2_name <- args$r2
pval_names <- c()
pval_list <- c()
cat_name_tracker <- c()
trait <- args$risk_trait
cat_split <- args$category_split
if(cat_split == '')
{
    select_cols <- c("IID", trait, unlist(as.list(strsplit(args$covars, '+', fixed = T)[[1]])))
}else{
    select_cols <- c("IID", trait,cat_split, unlist(as.list(strsplit(args$covars, '+', fixed = T)[[1]])))
}
phenos <- read_delim(args$pheno, delim = args$delim) %>% select(all_of(select_cols)) %>% na.omit()
print(paste0("Phenotype data for ", nrow(phenos), " individuals"))
#Set up if splitting by categories

if(cat_split == '')
{
    cat_names <- c('')
}else{
    cat_names <- unlist(unique(phenos %>% select(all_of(cat_split))))
    }
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
    for (cn in cat_names)
    {
        #event of no category splits
        if(cn =="") { 
            filt_list <- MI_pred
        } else {
            
            #filt_list <- filter(MI_pred, cat_split == cn)
            filt_list <- MI_pred[MI_pred[[cat_split]] == cn,]
        }
        
    if(nrow(unique(MI_pred[trait])) > 2 && !(args$case_control))
    {
            for (n in pval_names)
            {
                if(args$covars == "")
                {
                    #r2 <- c(r2, cor(MI_pred[trait], MI_pred[n]))
                    lmv <- lm(filt_list[[trait]] ~ filt_list[[n]])
                    r2 <- c(r2, summary(lmv)$r.squared)
                    pval_beta <- c(pval_beta, summary(lmv)$coefficients[2,4])
                }
                else{
                    
                    pname = paste0("\`",n,"\`")
                    f_full <- as.formula(paste(trait, paste(pname,"+", args$covars), sep = " ~ "))
                    lmv_complete <- lm(f_full, data = filt_list)
                    #lmv_corr <- lm(MI_pred[[trait]] ~ MI_pred[[n]] + MI_pred[[args$covars]])
                    f_part <-  as.formula(paste(trait, args$covars,sep = " ~ "))
                    lmv_null <- lm(f_part, data = filt_list)
                    #lmv_t <-  lm(MI_pred[[trait]] ~ MI_pred[[args$covars]])
                    r2 <- c(r2, summary(lmv_complete)$r.squared - summary(lmv_null)$r.squared)
                    pval_beta <- c(pval_beta, summary(lmv_complete)$coefficients[2,4])
                }
                cat_name_tracker <- c(cat_name_tracker, cn)
                pval_list <- c(pval_list, n)
            }
    } else {
        #Some other plots
        #ggplot(MI_pred, aes(x=c(0), y =n)) + geom_jitter(position=position_jitter(0.1), aes(fill = MI_pred[trait])) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        #ggsave(filename = paste0(args$output, n, "dotplot.png"))
        #print("Succesfully generated dotplot.")
        #}else #Its a case control trait
        #{
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
                r2 <- c(r2, R2ObservedNagelkerke(trait,p, filt_list))
                r2_name <- "Nagelkerke" 
            }
            else
            {
                r2 <- c(r2, R2LiabilityProbit(trait, p, filt_list))
                r2_name <- "Liability Scale (Probit)"
            }
            cat_name_tracker <- c(cat_name_tracker, cn)
            pval_list <- c(pval_list,p)
        }
    }
    #pval_names <- c(pval_names, str_extract(f, "[501]\\.*[\\-e\\d]+"))

    }
    dat <- data.frame("pval_names" = pval_list, r2, pval_beta, "category"=cat_name_tracker)
    write_tsv(dat, paste0(args$output,"_r2counts.tsv"))
    plotCorr(dat, args$output, r2_name,cat_split, args$hide_pvals)
    print(paste("Completed review for", f)) 
}


