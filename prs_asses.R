suppressMessages(library(stringr))
suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(RNOmni))

source("/work-zfs/abattle4/ashton/prs_dev/prs_tools/liability_pseudoR2.R")

parser <- ArgumentParser$new()
parser$add_description("Script to asses PRS results and their phenotypes. Can create bar plot of scores by different groups, report R2 performance for different p-values, and quantile plots.")
parser$add_argument("--prs_results", type = 'character', help = "File containing PRS scores.")
parser$add_argument("--pheno", type = "character", help = "Path to phenotype file, as given by Marios' analysis. Include covariates in this file too.", default = "/work-zfs/abattle4/ashton/prs_dev/CHS_race1_toy/pheno_mendrand.txt")
parser$add_argument("--output", type = "character", help = "Output file path.")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--risk_trait', help = "Specify the name of the condition we are scoring on, a string", type = "character", default = "LDL")
parser$add_argument("--case_control", type = "logical", action = "store_true", help = "Specify this if it isn't a continuous trait. Otherwise program will guess based on data", default = FALSE)
parser$add_argument("--delim", type = "character", help = "Specify which kind of delimiter in file, could be a ' ', '\\t' or otherwise", default = "\t")
parser$add_argument("--no_plot", help = "Specify this to say no plot to be included", action = "store_true", type = "logical", default = FALSE)
parser$add_argument("--r2", help = "Specify which type of r2 or pseudo-r2 to use. Default is either linear R2 or liability scale.", default = "linear", type = "character")
parser$add_argument("--covars", type = "character", help = "Specify column names with covariates to correct for, such as 'gender'. Default is none. Write as covar1+covar2+..", default = "")
parser$add_argument("--category_split", type = "character", help = "Specify different subgroups to score on, such as ancestry. This will score the groups separately.", default = "")
parser$add_argument("--hide_pvals", type = "logical", help = "Select this to omit pvalues on bar plots", default = FALSE, action = "store_true")
parser$add_argument("--irnt", type = "logical", help = "Select this to scale, center, and inverse-rank normalize the phenotype data.", default = FALSE, action = "store_true")
parser$add_argument("--quantile", type = "logical", help = "Select this to create a quantile plot", default = FALSE, action = "store_true")
parser$add_argument("--n_quants", type = "integer", help = "Specify the number of quantiles to put on quantile plot", default = 10)
parser$add_argument("--quant_style", type = "character", help = "Specify the type of quantile plot: simple (just PRS quantiles vs phenotype), or relative (PRS quantiles vs effect size, corrected for covariates, relative to median quantile)", default = "relative")
parser$helpme()
args <- parser$get_args()
options(readr.num_columns = 0)
#A function to plot the correlation bar plots

plotCorr <- function(dat, output, style_name, category_var, no_pvals)
{
    library("ggsci")
    dat$pval_names <- as.numeric(as.character(dat$pval_names))
    dat <- arrange(dat, pval_names)
    dat$pval_names <- as.factor(dat$pval_names)
    dat$logp <- -log10(dat$pval_beta)
    if(category_var == ""){
        base <- ggplot(dat, aes(x = pval_names, y = r2, fill =logp )) + geom_bar(stat = "identity") + 
            labs(x="P-value threshold", y = expression(paste(R ^ 2))) + labs(fill = "-log10(pval)") + scale_fill_gsea()
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
#Make a quantile plot
plotQuantile <- function(dat, trait, output, n_quants, style_name, covars_arg)
{
    covars <- str_split(covars_arg, fixed("+"))[[1]]
    for (n in pval_names) #is pval_names global here?
    {
        for_plot <-  dat %>% mutate(quantile = ntile(dat[[n]], n_quants)) 
        if(covars[1] != "")
        {for_plot <- for_plot %>% select(IID, all_of(n), quantile, all_of(trait), all_of(covars)) %>% filter(!(is.na(trait))) %>% ungroup()}
        else
        {for_plot <- for_plot %>% select(IID, all_of(n), quantile, all_of(trait)) %>% filter(!(is.na(trait))) %>% ungroup()}
        xounts <- count(for_plot, quantile)
        avgs <- for_plot %>% group_by(quantile) %>% dplyr::summarize(avg = mean(DBP), std_dev = sd(DBP)) %>% arrange(quantile) %>% mutate(num = xounts$n) %>% mutate(sem = std_dev/sqrt(num))
        #Simple one
        if (style_name == "simple")
        {
            plt <- (ggplot(avgs, aes(x=quantile, y=avg)) + geom_pointrange(aes(ymin=avg-sem, ymax=avg+sem)) + 
                        ylab(trait) + xlab("Quantile") + scale_x_continuous(breaks=c(1:n_quants), labels=c(1:n_quants)) + 
                        ggtitle(paste("Simple Quantile plot, p = ", n)) + scale_color_npg())
            
            ggsave(filename = paste0(output, ".quantile_plot.",n, ".", style_name, ".png"), plot = plt, height = 7, width = 8.5) 
        } else
        {
            #Nuanced one that allows for correcting for covariate- identify reference quantile- the median one.
            med_val <- which(for_plot[[n]] == median(for_plot[[n]]))
            med_quant <- for_plot[med_val,]$quantile
            beta <- c()
            std_error <- c()
            #$Compare each other quantile to the reference quantile
            for(i in 1:n_quants)
            {
                dat_c <- for_plot %>% filter(quantile %in% c(i, med_quant)) %>% mutate(explan  = ifelse(quantile == i, 1,0))
                #Still need to add lines for residual standard error, or stadnard error on beta, depending on what we choose.
                if(i == med_quant)
                {
                    beta <- c(beta, 0)
                    std_error <- c(std_error, 0)
                    next
                }
                if(covars[1] == "")
                {
                    lmv_prs <- lm(dat_c[[trait]] ~ dat_c$explan)
                }
                else{
                    
                    f_covars <- as.formula(paste(trait, paste(covars), sep = " ~ "))
                    lmv_covars <- lm(f_covars, data = dat_c)
                    lmv_prs <- lm(lmv_covars$residuals ~ dat_c$explan) #Could also do it all together and just pull out the betas. Would be worth looking empirically if they are the same, which sthey should be GIVEN independence
                }
                beta <- c(beta, summary(lmv_prs)$coefficients[2,1])
                std_error <- c(std_error, summary(lmv_prs)$coefficients[2,2])
            }
            quantile_covars <- data.frame(beta, std_error, quant = 1:n_quants) 
            plt <- (ggplot(quantile_covars, aes(x=quant, y=beta)) + geom_pointrange(aes(ymin=beta-std_error, ymax=beta+std_error)) + 
                        ylab("Effect size") + xlab("Quantile") +  scale_x_continuous(breaks=c(1:n_quants), labels=c(1:n_quants)) + 
                        ggtitle(paste("Relative Quantile plot, p = ", n, "Covariates:", covars_arg)))
            ggsave(filename = paste0(output, ".quantile_plot.",n, ".", style_name, ".png"), plot = plt, height = 7, width = 8.5) 
        }
       
    }
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
phenos <- read_delim(args$pheno, delim = args$delim) %>% select(unlist(select_cols)) %>% na.omit()

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
    prs$IID <- as.character(prs$IID)
    phenos$IID <- as.character(phenos$IID)
    print(paste0("PRS scores for ", nrow(prs), " individuals."))
    #full_dat contains the PRS scores at each p-value for each sample, their phenotypes, and the associated covariates to correct for.
    full_dat <- inner_join(prs, phenos, by = "IID")
    print(paste0(nrow(full_dat), " individuals had both PRS scores and phenotype data"))
    pval_names <- names(prs %>% select(-IID))
    full_dat <- na.omit(full_dat)
    #If the trait is continuous, do regular R2.
    for (cn in cat_names)
    {
         #event of no category splits
        if(cn =="") { 
            filt_list <- full_dat
        } else {
            
            #filt_list <- filter(full_dat, cat_split == cn)
            filt_list <- full_dat[full_dat[[cat_split]] == cn,]
        }
            
        if(nrow(unique(full_dat[trait])) > 2 && !(args$case_control)) #Trait is continuous
        {
                #Center and scale the data:
                if (args$irnt)
                {
                    filt_list[[paste0("unscaled_", trait)]] = filt_list[[trait]]
                    temp_ <- as.numeric(scale(filt_list[[trait]]))
                    filt_list[[trait]] <- rankNorm(temp_)
                }
                for (n in pval_names)
                {
                    if(args$covars == "")
                    {
                        #r2 <- c(r2, cor(full_dat[trait], full_dat[n]))
                        lmv <- lm(filt_list[[trait]] ~ filt_list[[n]])
                        r2 <- c(r2, summary(lmv)$r.squared)
                        pval_beta <- c(pval_beta, summary(lmv)$coefficients[2,4])
                    }
                    else{
                        
                        pname = paste0("\`",n,"\`")
                        f_full <- as.formula(paste(trait, paste(pname,"+", args$covars), sep = " ~ "))
                        lmv_complete <- lm(f_full, data = filt_list)
                        #lmv_corr <- lm(full_dat[[trait]] ~ full_dat[[n]] + full_dat[[args$covars]])
                        f_part <-  as.formula(paste(trait, args$covars,sep = " ~ "))
                        lmv_null <- lm(f_part, data = filt_list)
                        #lmv_t <-  lm(full_dat[[trait]] ~ full_dat[[args$covars]])
                        r2 <- c(r2, summary(lmv_complete)$r.squared - summary(lmv_null)$r.squared)
                        pval_beta <- c(pval_beta, summary(lmv_complete)$coefficients[2,4])
                    }
                    cat_name_tracker <- c(cat_name_tracker, cn)
                    pval_list <- c(pval_list, n)
                }
        } else {
            #Some other plots
            #ggplot(full_dat, aes(x=c(0), y =n)) + geom_jitter(position=position_jitter(0.1), aes(fill = full_dat[trait])) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
            #ggsave(filename = paste0(args$output, n, "dotplot.png"))
            #print("Succesfully generated dotplot.")
            #}else #Its a case control trait
            #{
            #This step is problematic.
            #full_dat[,4] <- as.factor(full_dat[,4])
            #ggplot(full_dat,aes(x=Score, color = trait)) + geom_histogram() +labs(x="PRS", y="Count", title="Distribution of PRS Scores")
            #ggsave(filename = paste0(args$output,substr(filename,1,nchar(f)-4), ".histogram.png"))
            #ggplot(full_dat, aes(x=c(0), y = Score, color = MI)) + geom_jitter(position=position_jitter(0.1)) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
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
    #write_tsv(dat, paste0(args$output,"_r2counts.tsv"))
    plotCorr(dat, args$output, r2_name,cat_split, args$hide_pvals)
}

#Do that for each individually. Now get the
dat <- data.frame(pval_names, r2)
write_tsv(dat, paste0(args$output,"_r2counts.tsv"))

if(args$quantile)
{
    plotQuantile(full_dat, trait, args$output, args$n_quants, args$quant_style, args$covars)
}
print("Finished plotting!")

