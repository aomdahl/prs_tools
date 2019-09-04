"""
Ashton Omdahl
August 2019
"""
#Make sure you are using the NON-PYMOL interpreter.
import argparse
import numpy as np
import scipy as sp
import plinkio 
import pandas as pd
from datetime import date
import sys

#General reference things
CHR = "chr"
POS = "pos"
REFID = "ref_id"
EFFAL = "effect_allele"
EFFALFREQ = "effect_allele_freq"
BETA = "beta"
PVAL = "pval"
REF= "REF" #reference allele
ALT = "ALT" #Alternate allele

def getStatHeaders(setting):
    if setting == "SAIGE":
        hl = "#CHROM POS ID REF ALT ac af num_cases num_controls beta sebeta Tstat pval pval_SAIGE_NoSPA Is_Converged varT varTstar".split(" ")
        return {CHR: hl[0], POS:hl[1], REFID: hl[2], EFFAL: hl[4], EFFALFREQ:hl[6], BETA:hl[9], PVAL:hl[12], ALT:hl[4], REF:hl[3]}
        #Can certainly add more later, but this is a fine start.
        
        
        
def loadSummaryStats(sum_path, sum_format, no_ambig):
    h = getStatHeaders(sum_format)
    #Make sure the delimiter is correct!
    fline=open(sum_path).readline().rstrip()
    t = fline.split('\t')
    if len(t) <= 1:
        print("We suspect the delimiter in your summary stats file is not a tab but something else. Please check and correct.")
        sys.exit()
    
    ss = pd.read_csv(sum_path, sep = '\t', index_col = False, 
                     dtype= { h[CHR]:'category', h[POS]:'int64', h[ALT]: 'category', h[REF]:'category', h[BETA]:'float64', h[PVAL]:'float64'})
    ss[REFID] = ss[h[CHR]].astype(str) + ":" + ss[h[POS]].astype(str)
    #Remove ambiguous SNPs
    if not no_ambig:
        print("Please make sure ambiguous alleles are removed. Program will terminate")
        sys.exit()
        #ss = ss.query((ALT != 'A' and REF != 'T') & (ALT != 'T' & REF != 'A') & (ALT != 'G' & REF != 'C') & (ALT != 'C' & REF != 'G'))
        #Note: these are ambiguous because in a comparison to another dataset its unclear if A is T visa versa, so including these cases with complementary ALT and REF alleles is ambiguous.
    return ss, h

def loadGenoVars(snp_file):
    ##CHROM  POS     ID      REF     ALT     FILTER  INFO
    bim_file = pd.read_csv(snp_file, sep = '\t', header = None, names = [CHR, POS, REFID, REF, ALT], dtype={CHR:'category', REFID:'str', POS:'int32', ALT: 'category', REF:'category'}, comment="#") 
    #its not just 24
    #bim_file = pd.read_csv(snp_file, sep = '\t', header = None, names = [CHR, POS, REFID, REF, ALT], dtype={CHR:'category', REFID:'str', POS:'int32', ALT: 'category', REF:'category'}, skiprows = list(range(0,24))) 
    return bim_file

    #Speed up the read in by specifying the type
    #Omit the first 23 lines- they are header information
    #Saves TONS of memory
    return bim_file

def filterSNPs(pval, maf, ss, vars, h):
    #@param h is the header data.
    #Get list of intersecting SNPs that pass the p-value threshold. Can add the MAF threshold at a later time.
    ss_new, vars_new = getIntersectingData(ss[(ss[h[PVAL]] <= pval) & (ss[h[CHR]] != "X")], vars)
    return ss_new[REFID], ss_new, vars_new

def calculatePRS(ss, geno_mat):
    id_list = ss[REFID].values
    snp_matrix = geno_mat.loc[id_list]
    #do the necessary math: choose all the right rows.
    #Assuming the alt are the affect alleles.. that is what we assume. Was true in the small number of cases I looked at.
    #Check all the values are correct....
    assert(sameAlleleAssesment(ss[ALT].values, snp_matrix[ALT].values))
    print("Same alt alleles!")
    assert(sameAlleleAssesment(ss[REF].values, snp_matrix[REF].values))
    print("Same REF alleles!")
    #get the ALT column of the SNP matrix and make sure it matches with the alleles for which we have effects.
    dat = (2-snp_matrix.iloc[:,5:]).values
    betas = stats_filtered[BETA].values
    scores = np.matmul(betas, dat)
    return scores

def sameAlleleAssesment(s1, s2):
    return (s1 == s2).all


#TODO: This is slow, probably because of the sorting. Come up with a way to maintain order and speed.
#Checked functionality theoretically, seems to work.
def getIntersectingData(ss, vars):
    id_keep = ss[REFID]
    #select out the values in vars that are in the vars_keep list (make sure IN ORDER)
    vars_keep = vars[vars[REFID].isin(id_keep)].sort_values(by=REFID) #Remember, if its a 3 we ignore it- it means the data isn't available for them.
    return ss[ss[REFID].isin(vars_keep[REFID])].sort_values(by=REFID), vars_keep 

def plinkToMatrix(snp_keep, args):
    #Write out the list of SNPs to keep
    snp_order = list(snp_keep.values)
    list_dir = writeSNPIDs(snp_keep)
    #Make the new matrix
    from subprocess import call 
    call(["plink2", "--pfile", args, "--not-chr", "X", "--extract", list_dir, "--out", "mat_form_tmp", "--export", "A-transpose", "--max-alleles", "2"]) 
    
    ret = pd.read_csv("mat_form_tmp.traw", sep = "\t", index_col="SNP")
    ret.rename(columns={"COUNTED":REF}, inplace=True)
    return ret

def getPatientIDs(snp_matrix):
    tnames = list(snp_matrix.columns)[5:]
    return [int(i.split("_")[0]) for i in tnames]


def writeScores(scores, ids, destination):
    with open(destination, 'w') as ostream:
        ostream.write("ID\tScore\n")
        for s, id, in zip(scores, ids):
            ostream.write(str(id) + '\t' + str(s) + '\n')

    return

def writeSNPIDs(ids):
    ids.to_csv("snp_ids.tmp", header = False, index = False, sep = '\t') #(should have 2 columns of the family/infamily ids, which appear to be the same.)
    return "snp_ids.tmp"

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Basic tools for calculating rudimentary Polygenic Risk Scores (PRSs). This uses PLINK bed/bim/fam files and GWAS summary stats as standard inputs. Please use PLINK 2, and ensure that there are no multi-allelic SNPs")
    parser.add_argument("-snps", "--plink_snps", help = "Plink file handle for actual SNP data (omit the .extension portion)")
    parser.add_argument("-ss", "--sum_stats", help = "Path to summary stats file. Be sure to specify format if its not DEFAULT. Assumes tab delimited", required = True)
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE"])
    parser.add_argument("-p", "--pval", default = 5e-8, type = float, help = "Specify what pvalue threshold to use.")
    parser.add_argument("--pvals", help= "Use this if you wish to specify multiple at once, in a list separated by commas")
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--no_ambig", default = False, action = "store_true", help = "Specify this option if you have already filter for ambiguous SNPs and bi-alleleic variants. Recommended for speed.")
    args = parser.parse_args()
    pvals = [args.pval]
    if args.pvals:
        pvals = [float(x) for x in (args.pvals).split(",")]
        
    print("Loading summary statistics and relevant genotype information....")
    stats, header = loadSummaryStats(args.sum_stats, args.ss_format, args.no_ambig)
    variants = loadGenoVars(args.plink_snps + ".pvar")
    print("Filtering SNPs...")
    for pval in pvals:
        snp_list, stats_filtered, variants_filtered = filterSNPs(pval, args.maf, stats, variants, header) #TODO: makes sure the variants are in order.
        #Run plink filtering on the bed file, and then load it into memory
        print("Building a PLINK matrix for pval", str(pval), "...")
        snp_matrix = plinkToMatrix(snp_list, args.plink_snps)
        #stats was output of filterSNPs
        print("Calculating PRSs")
        #stats_qc, snp_qc = qualityControl(stats_filtered, snp_matrix)
        scores = calculatePRS(stats_filtered, snp_matrix)
        patient_ids  = getPatientIDs(snp_matrix)
        writeScores(scores, patient_ids, args.output + "_" + str(pval) + ".tsv")
        #Match SNPs by address and/or reference ID.
        print("Scores written out to", args.output + "_" + str(pval) + ".tsv")
        #add option to accomodate large size files with limited memory- chunk it.


