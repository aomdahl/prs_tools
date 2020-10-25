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
source ~/.bashrc
set -e
echo "Going to run a paired analysis on $ANNOT, for trait $TRAIT"
ml python/3.7-anaconda
source ~/.bashrc
#Step 1: create a generous clumped dataset to sample from
python $PRS  -snps $TARGET  -ss $SS -o ./$ANNOT -pv 2 --ss_format NEAL --preprocess_only --clump --ld_ref /work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref --no_na --reset_ids --align_ref

#Step 2: Select from those snps the ones for which we have annotations
#Requires the annotation database name and column number
python /work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/annot_extract_cleanup/annotExtractionLite.py --output extract_list.tsv --search geno_ids.f --annot $ANNOT --database $DB  

#Get the median score in the extract list
##need to fix this function, not working
if [[ "$ANNOT_TYPE" == "continuous" ]]; then
    MED=`cut -f 2 extract_list.tsv | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`
    #Split into upper and lower lists with awk.

    awk -v med="$MED" '($2 + 0.0  > med) {print $0}' extract_list.tsv > high_${ANNOT}.txt
    awk -v med="$MED" '($2 + 0.0  < med) {print $0}' extract_list.tsv > low_${ANNOT}.txt
fi
if [[ "$ANNOT_TYPE" == "binary" ]]; then
    
    awk -v med="$MED" '($2 + 0.0  > 0) {print $0}' extract_list.tsv > high_${ANNOT}.txt
    awk -v med="$MED" '($2 + 0.0  == 1) {print $0}' extract_list.tsv > low_${ANNOT}.txt
fi    
awk '{print $1}' extract_list.tsv > combined_snps.txt
#
##DEBUGGING:
##Check that the distributions split correctly
cut -f 2 high_${ANNOT}.txt > upper_vals.tmp
cut -f 2 low_${ANNOT}.txt > lower_vals.tmp
ml gcc/5.5.0
ml R
Rscript /work-zfs/abattle4/ashton/gene_utils/quick_hist.R --input upper_vals.tmp --output upper_${ANNOT}.png --metric NONE
Rscript /work-zfs/abattle4/ashton/gene_utils/quick_hist.R --input lower_vals.tmp --output lower_${ANNOT}.png --metric NONE


cut -f 1 high_${ANNOT}.txt > upper.tmp && mv upper.tmp high_${ANNOT}.txt
cut -f 1 low_${ANNOT}.txt > lower.tmp && mv lower.tmp low_${ANNOT}.txt
#
##Pair these by p-value. Might also be worth looking into a maf-pairing approach.
ml gcc/5.5.0
ml R
Rscript /work-zfs/abattle4/ashton/prs_dev/prs_tools/annotation_analysis/assign_paired_variants.R --lower low_${ANNOT}.txt --upper high_${ANNOT}.txt --sum_stats $SS --null combined_snps.txt --output ./split_lists --inverted_snp_encoding 

#Calculate scores on each list subset. May actually be faster to use the 
#Note that here, I am runnign each separately. This allows me to get things moving a bit faster
#however, at some point you may want to run these steps combined into 1. Then you should do the following:
#python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ./extract_list.tsv --split_scores  ../split_lists.paired_high_annots.tsv, ../split_lists.paired_low_annots.tsv

mkdir -p high_annot
cd high_annot
pwd
python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ../split_lists.paired_high_annots.tsv --inverted_snp_encoding --reset_ids &
#Run the low_annot version
cd ../
mkdir -p low_annot
cd low_annot
python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ../split_lists.paired_low_annots.tsv --inverted_snp_encoding --reset_ids &
#Run the standard background version
cd ../
#aside- will there be an issue because geno_ids not in the local directory? I hope not:
mkdir -p std 
cd std
python $PRS  -snps ../new_plink  -ss ../filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ../split_lists.paired_null_annots.tsv --inverted_snp_encoding --reset_ids & 
cd ../ 
wait
#Now asses them all
Rscript $PRS_ASSES --prs_results ./low_annot/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_low --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3
Rscript $PRS_ASSES --prs_results ./std/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_std --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3
Rscript $PRS_ASSES --prs_results ./high_annot/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_high --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3

Rscript $COMBINE_ASSES --high_results ${ANNOT}_high_r2counts.tsv --low_results ${ANNOT}_low_r2counts.tsv --output ${ANNOT}_multi.png --annot_name ${ANNOT} --mixed_results ${ANNOT}_std_r2counts.tsv
#^Modify this to include a third option of results if you'd like (is it worth it if you're dropping the project??)
echo "Finished!"
