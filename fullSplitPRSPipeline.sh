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

python $PRS  -snps $TARGET  -ss $SS -o ./$ANNOT -pv 2 --ss_format NEAL --preprocess_only --clump --ld_ref /work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref --no_na --reset_ids --select_var $enriched  --align_ref
#We lost a TON of variants through clumping. Not really sure why or what happened, but here we are.
#Preprocess through clumping separately, get the max number of SNPS possible.
Rscript /work-zfs/abattle4/ashton/prs_dev/prs_tools/annotation_analysis/assign_paired_variants.R $background ./geno_ids.f $SS ./split_lists

python $PRS  -snps new_plink  -ss ./filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --no_overwrite --select_vars ./split_lists.paired_high_annots.tsv

Rscript $PRS_ASSES --prs_results $ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output $ANNOT --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3

#Now run a comparable background
echo "Running background PRS. This assumes you have a clumped subset already to go."
enriched_dir=`pwd`
cd $background_dir
python $PRS -snps new_plink  -ss ./filtered_NEAL.ss -o ./${ANNOT}_paired -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --no_overwrite --select_vars $enriched_dir/split_lists.paired_low_annots.tsv

Rscript $PRS_ASSES --prs_results ${ANNOT}_paired.full_prs.scores.tsv --pheno $PHENO --output paired_with_${ANNOT} --risk_trait BMI --covars age+gender+PC1+PC2+PC3
cd $enriched_dir
Rscript $COMBINE_ASSES --high_results ${ANNOT}_r2counts.tsv --low_results ${background_dir}/paired_with_${ANNOT}_r2counts.tsv --output ${ANNOT}_vs_std.png --annot_name ${ANNOT}

echo "Finished!"
