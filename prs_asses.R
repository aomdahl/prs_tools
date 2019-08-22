library(magrittr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
'%ni%' <- Negate('%in%')

prs <- read_tsv("/work-zfs/abattle4/ashton/prs_dev/scratch/calc_prs_testing/beta_results_test.tsv")
phenos <- read_tsv("/work-zfs/abattle4/ashton/prs_dev/CHS_race1_toy/pheno_mendrand.txt") %>% rename(ID = IID)
MI_pred <- inner_join(prs, phenos) %>% select(ID, Score, MI)
MI_pred$MI <- as.factor(MI_pred$MI)
ggplot(MI_pred,aes(x=Score, color = MI)) + geom_histogram() +labs(x="PRS", y="Count", title="Distribution of PRS Scores")

ggplot(MI_pred, aes(x=c(0), y = Score, color = MI)) + geom_jitter(position=position_jitter(0.1)) + xlim(-0.4, 0.4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
