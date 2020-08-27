#!/bin/bash
#template script for running the prs calculator
#Please be sure you have installed plink2: 
set -e
ml python/3.7-anaconda-2019.03
SC="/work-zfs/abattle4/ashton/prs_dev/"
GENO="/work-zfs/abattle4/marios/HF_GWAS/cohorts/CHS/merged_genotypes/CNV370v1/preimpute/race1/unrelated/unrelated_pfile"
SUMSTAT="/work-zfs/abattle4/ashton/prs_dev/dbp_ss_neal/4079_irnt.gwas.imputed_v3.both_sexes.tsv"
LDREF="/work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref"
PHENO="/work-zfs/abattle4/marios/HF_GWAS/cohorts/CHS/merged_genotypes/CNV370v1/preimpute/race1/unrelated/pheno_mendrand.txt"
#To run the tool in just one step:
python $SC/prs_tools/calc_prs_lowmem.py -snps $GENO -ss $SUMSTAT --pvals 1,0.5,0.3,0.1,0.01,0.001,1e-4,1e-6,5e-8  -o results_jan9 -pv 2 --reset_ids --ss_format NEAL --align_ref --clump --ld_ref $LDREF --debug
 

#You can separately run the preprocessing to make sure everything is working alright, but this isn't necessary
python $SC/prs_tools/calc_prs_lowmem.py -snps $GENO -ss $SUMSTAT --pvals 0.1,0.01,0.001,1e-4,1e-6,5e-8  -o results_jan9 -pv 2 --reset_ids --ss_format NEAL --align_ref --clump --ld_ref $LDREF --only_preprocess 
    #You can check the outputs from above to see if they look right. 

#Now run the actual calculation step
python $SC/prs_tools/calc_prs_lowmem.py -snps new_plink -ss filtered_NEAL.ss --pvals 1,0.5,0.1,0.01,0.001,1e-4,1e-6,5e-8  -o results_jan9 -pv 2 --reset_ids --ss_format DEFAULT --prefiltered_ss --preprocessing_done
    #This should generate your output file! Note that many of the file names changed, since they are filtered and printed out for convenience inthe step above. In the future, I may set the outputs just into a serialized file for speed, but this version is convenient for debugging.

#To asses the performance of the scoring:
PA="/work-zfs/abattle4/ashton/prs_dev/prs_tools"
module load gcc/5.5.0
module load R
Rscript $PA/prs_asses_lowmem.R  --prs_results results_jan9.tsv --pheno $PHENO --output figs_chs_race1 --risk_trait DBP --covars gender
    #this will generate a plot showing the r2 of scores with DBP, as well as a matrix, corr_count.tsv, that shows the raw r2 scores

