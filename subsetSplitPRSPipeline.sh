#!/bin/bash
#Divided PRS pipeline
#inputs:
#1: Setup source file where everything goes.
#paired PRS run for comparison:
##Define vars 
PRS="/work-zfs/abattle4/ashton/prs_dev/prs_tools/calc_prs.py"
PRS_ASSES="/work-zfs/abattle4/ashton/prs_dev/prs_tools/prs_asses.R"
COMBINE_ASSES="/work-zfs/abattle4/ashton/prs_dev/prs_tools/combinePRSAssesments.R"
source $1
#mkdir -p $enriched_dir
cd $enriched_dir
set -e
echo "Going to run a paired analysis on $ANNOT, for trait $TRAIT, with background PRS directory $background_dir"
ml python/3.7-anaconda
ml gcc/5.5.0
ml R
#Oct 12: after some discussion, modifying the pipeline from being annotation selection first to doing clumping first, and then using that as the base list for analysis

#Step 1: create a generous clumped dataset to sample from
python $PRS  -snps $TARGET  -ss $SS -o ./$ANNOT -pv 2 --ss_format NEAL --preprocess_only --clump --ld_ref /work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref --no_na --reset_ids --align_ref

#Step 2: Select from those snps the ones for which we have annotations
#Requires the annotation database name and column number
*************
python extractAnnots.py --output extract_list.tsv 

#Get the median score in the extract list

MED=`getMedian extract_list.tsv 2`

#Split into upper and lower lists with awk.
awk -v med="$MED" '($2 + 0.0  > med) {print $1}' extract_list.tsv > high_${ANNOT}.txt
awk -v med="$MED" '($2 + 0.0  < med) {print $1}' extract_list.tsv > low_${ANNOT}.txt
awk '{print $1}' extract_list.tsv > combined_snps.txt

#Pair these by p-value. Might also be worth looking into a maf-pairing approach.

Rscript /work-zfs/abattle4/ashton/prs_dev/prs_tools/annotation_analysis/assign_paired_variants.R --lower low_${ANNOT}.txt --upper high_${ANNOT}.txt --sum_stats $SS --null combined_snps.txt --output ./split_lists --inverted_snp_encoding 

#Calculate scores on each list subset. May actually be faster to use the 
#Note that here, I am runnign each separately. This allows me to get things moving a bit faster
#however, at some point you may want to run these steps combined into 1. Then you should do the following:
#python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ./extract_list.tsv --split_scores  ../split_lists.paired_high_annots.tsv, ../split_lists.paired_low_annots.tsv

mkdir -p high_annot
cd high_annot
python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ../split_lists.paired_high_annots.tsv --inverted_snp_encoding &
#Run the low_annot version
cd ../
mkdir -p low_annot
cd low_annot
python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ../split_lists.paired_low_annots.tsv --inverted_snp_encoding &

#Run the standard background version
cd ../
#aside- will there be an issue because geno_ids not in the local directory? I hope not:
mkdir -p std 
python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ../_lists_combined.tsv --inverted_snp_encoding &
cd ../
wait
#Now asses them all
Rscript $PRS_ASSES --prs_results ./low_annot/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_low --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3
Rscript $PRS_ASSES --prs_results ./std_annot/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_std --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3
Rscript $PRS_ASSES --prs_results ./high_annot/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_high --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3

Rscript $COMBINE_ASSES --high_results ./high_annot/${ANNOT}_r2counts.tsv --low_results ./low_annot/${ANNOT}_r2counts.tsv --output ${ANNOT}_multi.png --annot_name ${ANNOT} --mixed_results ./std_annot/${ANNOT}_r2counts.tsv
#^Modify this to include a third option of results if you'd like (is it worth it if you're dropping the project??)
echo "Finished!"
