suppressMessages(library(stringr))
suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(rms))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
 #Data should be clumped...
parser <- ArgumentParser$new()
parser$add_description("R script for generating figures associated with PRS results and their phenotypes.")
parser$add_argument("-o", "--output", type = "character", help = "Output file name")
parser$add_argument("--gwas", type = "character", help = "GWAS summary statistics file. Assuming default format")
parser$add_argument("--include", type = "character", help = "Path to list of pre-clumped SNPs to consider. Recommend that these are independent, not in LD snps.")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument("--upper_bound", type = "numeric", help = "Specify the highest value on output plot for beta", default = 0.05)
parser$helpme()
args <- parser$get_args()


filterList <- function(filt_path, gwas)
{
  filt_list <- read_tsv(filt_path,col_names = FALSE)
  gwas %>% filter(ID %in% filt_list$X1)
}

gwas_in <- fread(args$gwas)

if(length(args$include) > 0)
{
  print("Filtering by given list...")
  gwas_in <- filterList(args$include, gwas_in)
}

print("Data has been read in...")
#Get the frequency at 

gnmwide_thresh <- gwas_in %>% drop_na() %>% filter(pval <= 5e-8)
print("first filter...")

betas_thresh_marginal <- data.frame("betas" = gnmwide_thresh$beta,"betas_abs" = abs(gnmwide_thresh$beta), thresh = as.factor(5e-8))
omit_list <- gnmwide_thresh$ID
#thresh_vals <- c(1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
thresh_vals <- c(1e-4, 0.01, 0.1, 0.3, 0.5, 1)
prev = 5e-8
for(thresh in thresh_vals)
{
  print(paste("Currently working on", as.character(thresh)))
  dat <- (gwas_in %>% filter(pval <= as.numeric(thresh)) %>% filter(pval > prev))
  tmp <- data.frame("betas" = dat$beta,"betas_abs" = abs(dat$beta), thresh = as.factor(thresh))
  betas_thresh_marginal <- rbind(betas_thresh_marginal, tmp)
  prev = as.numeric(thresh)
}

print("Finished iterations")
#Add lines for mean of each group
means <- betas_thresh_marginal %>% group_by(thresh) %>% dplyr::summarise(mu = median(betas_abs))
ggplot(betas_thresh_marginal, aes(betas_abs,fill = thresh)) + geom_density(alpha = 0.1) + scale_x_continuous(limits = c(0, args$upper_bound)) + geom_vline(data = means, aes(xintercept=mu, color = thresh)) 
ggsave(args$output, height = 6, width = 8)
