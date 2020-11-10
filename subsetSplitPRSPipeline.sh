#!/bin/bash
#Divided PRS pipeline
#inputs:
#1: Setup source file where everything goes.
#paired PRS run for comparison:
##Define vars 
source $1
source ~/.bashrc
set -e
echo "Going to run a comparative PRS analysis on $ANNOT, for trait $TRAIT"
ml python/3.7-anaconda
source ~/.bashrc
#REplacing slurm_run.sh with ${SLURM_PATH}
#Step 1: create a generous clumped dataset to sample from, unless a list has been already specified
if [ "$PRECLUMP" == "NA" ]; then
  echo "Creating a clumped target of SNPs to score on...."
  python $PRS  -snps $TARGET  -ss $SS -o ./$ANNOT -pv 2 --ss_format NEAL --preprocess_only --clump --ld_ref /work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref --no_na --reset_ids --align_ref --maf $MAF_THRESH
  CLUMPLIST="geno_ids.f"
  PRECLUMP=`pwd`
else
  CLUMPLIST="${PRECLUMP}/geno_ids.f"
  echo "SNPs for scoring have been provided in ${CLUMPLIST}"
fi

#Step 2: Select from those snps the ones for which we have annotations
#Requires the annotation database name and column number
FNH="divided_${ANNOT}"
python /work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/annot_extract_cleanup/annotExtractionLite.py --output ${FNH} --database ${DB} --annot ${ANNOT} --split_metric quantiles --search ${CLUMPLIST} --parts $SD --quantiles ${N_QUANT}


ml gcc/5.5.0
ml R
#Types of analysis approach options
if [ "$PAIRED" == "TRUE" ]; then #paired, which means # SNPs in scoring will be matched either by P-value or MAF
  if [ "$BG_ITER" != "0" ]; then #Iterate on a background set of SNPs, to get a sense for overall performance
  #1) if you want to do all paired vs a backgroudn with iteration
    Rscript /work-zfs/abattle4/ashton/prs_dev/prs_tools/annotation_analysis/assign_paired_variants.R --lower ${FNH}_low_annot.txt --upper ${FNH}_high_annot.txt --sum_stats $SS --null ${FNH}_extract_list.tsv --output ./split_lists --null_rep $BG_ITER  --pair_metric $PAIRING
  else
  #2) If you want paired with the full background (for comparison and enrichment calculations.)
    Rscript /work-zfs/abattle4/ashton/prs_dev/prs_tools/annotation_analysis/assign_paired_variants.R --lower ${FNH}_low_annot.txt --upper ${FNH}_high_annot.txt --sum_stats $SS --output ./split_lists --pair_metric $PAIRING
  fi
  #Specify where the high and low subset lists are.
  HIGH_VARS="split_lists.paired_high_annots.tsv"
  LOW_VARS="split_lists.paired_low_annots.tsv"
else #Paired = false, doesn't allow for background permutations
    #3) If you want no paired with just the full background (for enrichment calculations)
    HIGH_VARS="${FNH}_high_annot.txt"
    LOW_VARS="${FNH}_low_annot.txt"
fi

#Actually calculating the PRS
echo "Launching calculation of PRS for high annotation set..."
mkdir -p high_annot
cp ${SLURM_PATH} ./high_annot/slurm_run_high.sh
cd high_annot
echo "source ../${1}" >> slurm_run_high.sh
echo "python $PRS  -snps ${PRECLUMP}/new_plink  -ss "${PRECLUMP}/filtered_NEAL.ss" -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --prefiltered_ss --select_vars ../${HIGH_VARS} --inverted_snp_encoding --reset_ids" >> slurm_run_high.sh
sbatch slurm_run_high.sh
#Run the low_annot version
cd ../

echo "Launching calculation of PRS for low annotation set..."
mkdir -p low_annot
cp ${SLURM_PATH} ./low_annot/slurm_run_low.sh
cd low_annot
echo "source ../${1}" >> slurm_run_low.sh
echo "python $PRS -snps ${PRECLUMP}/new_plink  -ss ${PRECLUMP}/filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --prefiltered_ss --select_vars ../split_lists.paired_low_annots.tsv --inverted_snp_encoding --reset_ids" >> slurm_run_low.sh
sbatch slurm_run_low.sh
cd ../


echo "Launching calculation of PRS for standard annotation set..."
mkdir -p std
cd std
#If its a non-iterating one..... (2 or 3)
if [ "$BG_ITER" == "0" ]; then
  cp ${SLURM_PATH} ./slurm_std_run.sh #Beware, this may need more time.
  echo "source ../${1}" >> slurm_std_run.sh
  echo "python $PRS  -snps ${PRECLUMP}/new_plink  -ss ${PRECLUMP}/filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --prefiltered_ss --no_na --select_vars ../${FNH}_extract_list.tsv --inverted_snp_encoding --reset_ids"  >> slurm_std_run.sh
  sbatch ./slurm_std_run.sh
else
  #Run the standard background version
    for (( i=1; i<= $BG_ITER; i++)); do
        mkdir -p ${i}_std
        cd ${i}_std
        cp ${SLURM_PATH} ./slurm_${i}_std_run.sh
        echo "source ../../${1}" >> slurm_${i}_std_run.sh
        echo "python $PRS  -snps ${PRECLUMP}/new_plink  -ss ${PRECLUMP}/filtered_NEAL.ss -o ./$ANNOT -pv 2 --ss_format DEFAULT --preprocessing_done --no_na --select_vars ../../split_lists${i}.paired_null_annots.tsv --inverted_snp_encoding --reset_ids"  >> slurm_${i}_std_run.sh
        sbatch slurm_${i}_std_run.sh
        cd ../
    done
fi
  cd ../


#Now asses them all
until [[ -f ./low_annot/$ANNOT.full_prs.scores.tsv && -f ./low_annot/$ANNOT.full_prs.scores.tsv ]]
do
     sleep 5
done
echo "High/low PRS calculation done. Performing assesments..."
Rscript $PRS_ASSES --prs_results ./low_annot/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_low --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3
Rscript $PRS_ASSES --prs_results ./high_annot/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_high --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3


if [ "$BG_ITER" != "0" ]; then
  until [[ -f ./std/1_std/$ANNOT.full_prs.scores.tsv && -f ./std/${BG_ITER}_std/$ANNOT.full_prs.scores.tsv ]]
  do
      sleep 5
  done
  cd std
  
for (( i=1; i<= $BG_ITER; i++)); do
    cd ${i}_std
    Rscript $PRS_ASSES --prs_results $ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ../${i}_${ANNOT}_std --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3
    cd ../
  done
  cd ../
else #No iterations
  until [ -f ./std/$ANNOT.full_prs.scores.tsv ]
    do
        sleep 5
    done
  Rscript $PRS_ASSES --prs_results ./std/$ANNOT.full_prs.scores.tsv --pheno ${PHENO} --output ${ANNOT}_std --risk_trait $TRAIT --covars age+gender+PC1+PC2+PC3
fi
echo "Background assesments done."

#Final step- combining the assesments
COMBINE_ASSES="/work-zfs/abattle4/ashton/prs_dev/prs_tools/combinePRSAssesments.R"
echo "All assesments complete; combining analysis...."
#All versions- just compare the high and low
Rscript $COMBINE_ASSES --high_results ${ANNOT}_high_r2counts.tsv --low_results ${ANNOT}_low_r2counts.tsv --output ${ANNOT}_high_low.png --annot_name ${ANNOT}

if [ "$BG_ITER" != "0" ]; then
#If you did the permuted version:
Rscript $COMBINE_ASSES --high_results ${ANNOT}_high_r2counts.tsv --low_results ${ANNOT}_low_r2counts.tsv --output ${ANNOT}_mixed.png --annot_name ${ANNOT} --mixed_results ./std/ --multi_mixed "*_std_r2counts.tsv"
else
  Rscript $COMBINE_ASSES --high_results ${ANNOT}_high_r2counts.tsv --low_results ${ANNOT}_low_r2counts.tsv --output ${ANNOT}_mixed.png --annot_name ${ANNOT} --mixed_results ./${ANNOT}_std_r2counts.tsv
fi
