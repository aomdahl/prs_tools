#!/usr/bin/env Rscript
#Quick script to calculate scores at p = 1
#arguments are 1: matrix file in, 2: weights with ids as keys, 3: output file location
args = commandArgs(trailingOnly=TRUE)
fin = args[1]
weights = args[2]
output = args[3]
library(readr)
library(data.table)

lower <- fread(fin)
scores <- 2- lower[,7:ncol(lower)]

weights <- fread(weights)
colnames(weights) <- c("IDS", "betas")
names_low <- colnames(scores)
names_low <- data.frame("IDS" = sub("_[ACTG]", "", names_low))
combined <- merge(x = names_low, y = weights, by = "IDS", all.x = TRUE)
reorder_attempt <- combined[names_low$IDS,]
print(paste("ID errors:",sum(reorder_attempt$IDS != names_low$IDS))) #Its 0, we are good!
betas <- as.matrix(reorder_attempt$betas)
geno_m <- as.matrix(scores)
scores_low <- round((geno_m %*% betas), digits = 5)
low_out <- data.frame("IDs" = lower[,1], "p=1" = scores_low)
write_tsv(low_out, output)
