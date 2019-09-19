# prs_tools
Basic tools for performing polygenic risk score (PRS) analysis.
Developed as part of my research at JHU in the lab of Alexis Battle, starting in August of 2019.
Scripts:
    -liability_pseudoR2.R: tools for calculating psuedo R2, which is used for case/control studies. Based on the implementation/pseudocode given by Lee et al in 2012 Genetic Epidemiology Paper, "A Better Coefficient of Determination for Genetic Profile Analysis"
    -prs_asses.R: tool to asses the quality of calculated PRS scores; outputs various figures (scatter plot, bar chart) to show score distribution and correlation. Requires liability_pseudoR2.R in same directory to run.
    -calc_prs.py: Tool to perform Polygenic Risk score (PRS) analysis. Call `python calc_prs.py -h` for usage. This is the original version of the tool, with a pandas-based approach
    -calc_prs_dask.py: Early attempts at making the tool using dask.
    -calc_prs_speed.py: PRS calculations (attempts) using dask and one-time id filtering.
    -calc_prs_lowmem.py: PRS calcuations that parse the file in a line-by-line manner
