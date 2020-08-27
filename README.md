# prs_tools
Basic tools for performing polygenic risk score (PRS) analysis.
Developed as part of my research at JHU in the lab of Alexis Battle, starting in August of 2019.\
\
**Scripts:**
* *calc_prs.py*: Tool to perform Polygenic Risk score (PRS) analysis. Call `python calc_prs.py -h` for usage. This is the current version of the tool, which uses plink to write out plink genotype files and parse them line-by-line.
* *liability_pseudoR2.R*: tools for calculating psuedo R2, which is used for case/control studies. Based on the implementation/pseudocode given by Lee et al in 2012 Genetic Epidemiology Paper, "A Better Coefficient of Determination for Genetic Profile Analysis"
* *prs_asses.R*: tool to asses the quality of calculated PRS scores; outputs various figures (bar chart, quantile plot) to show score distribution and correlation. Requires liability_pseudoR2.R in same directory to run.
* *run_prs_template.sh*: Example template of how to run the basic pipeline.
* *compare_prs_results.R*: Script to directly compare the outputs of 2 calc_prs.py runs. Helpful for debugging.



* **annotation_analysis:**
  * *top_annots_script.R*: Regress a set of genomic annotations against GWAS summary statistics. Optimized for the NEALE lab GWAS sumstats format.
  * *assign_paired_variants.R*: Script to pair SNPS by GWAS pvalue to generate 2 sets of SNPS of ~equal size.
* **old_versions_calc_prs:**
  * *calc_prs_dask.py* (OBSOLETE): Early attempts at making the tool using dask. 
  * *calc_prs_speed.py* (OBSOLETE): PRS calculations (attempts) using dask and one-time id filtering.
  * *calc_prs_lowmem.py* (OBSOLETE): PRS calcuations that parse the file in a line-by-line manner
