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
        
        
        
def loadSummaryStats(sum_path, sum_format):
    h = getStatHeaders(sum_format)
    ss = pd.read_csv(sum_path, sep = ' ', index_col = False, 
                     dtype={h[CHR]:'category', h[POS]:'int32', h[ALT]: 'category', h[REF]:'category', h[BETA]:'float64', h[PVAL]:'float64'})
    ss[REFID] = ss[h[CHR]].astype(str) + ":" + ss[h[POS]].astype(str)
    #Remove biallelic SNPs
    #Remove ambiguous SNPs
    return ss, h

def loadGenoVars(snp_file):
    bim_file = pd.read_csv(snp_file, sep = '\t', header = None, names = [CHR, POS, REFID, REF, ALT], dtype={CHR:'category', REFID:'str', POS:'int32', ALT: 'category', REF:'category'}, skiprows = list(range(0,24))) 
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
    #get the ALT column of the SNP matrix and make sure it matches with the alleles for which we have effects.
    assert(sameAlleleAssesment(ss[ALT].values, geno_mat[ALT].values))
    assert(sameAlleleAssesment(ss[REF].values, geno_mat[REF].values))
    dat = (2-snp_matrix.iloc[:,5:]).values
    betas = stats_filtered[BETA].values
    scores = np.matmul(betas, dat)
    return scores

def sameAlleleAssesment(s1, s2):
    sames = (s1 == s2)
    if sames.all():
        return True
    else:
        print("There are differences. Investigate further (line 73)")
        sys.exit()

#TODO: This is slow, probably because of the sorting. Come up with a way to maintain order and speed.
def getIntersectingData(ss, vars):
    id_keep = ss[REFID]
    vars_ids = vars[REFID].sort_values()
    #select out the values in vars that are in the vars_keep list (make sure IN ORDER)
    vars_keep = vars[vars[REFID].isin(id_keep)].sort_values(by=REFID) #Remember, if its a 3 we ignore it- it means the data isn't available for them.
    return ss[ss[REFID].isin(vars_keep[REFID])].sort_values(by=REFID), vars_keep, 
    


def plinkToMatrix(snp_keep, args):
    #Write out the list of SNPs to keep
    snp_order = list(snp_keep.values)
    list_dir = writeSNPIDs(snp_keep)
    #Make the new matrix
    from subprocess import call 
    call(["plink2", "--pfile", args, "--not-chr", "X", "--extract", list_dir, "--out", "mat_form_tmp", "--export", "A-transpose", "--biallelic-only", "strict"]) 
    return pd.read_csv("mat_form_tmp.traw", sep = "\t", index_col="SNP")

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

    parser = argparse.ArgumentParser(description = "Basic tools for calculating rudimentary Polygenic Risk Scores (PRSs). This uses PLINK bed/bim/fam files and GWAS summary stats as standard inputs. Please use PLINK 2, and ensure that there are no multi-allelic SNPs (PLINK biallelic strict))
    parser.add_argument("-snps", "--plink_snps", help = "Plink file handle for actual SNP data (omit the .extension portion)")
    parser.add_argument("-ss", "--sum_stats", help = "Path to summary stats file. Be sure to specify format if its not DEFAULT. Assumes tab delimited")
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE"])
    parser.add_argument("-p", "--pval", default = 5e-8, type = float, help = "Specify what pvalue threshold to use.")
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--train", action = "store_true", help = "Select this if training a PRS model. This will save a list of the SNPs we include in the test, and can also draw some plots (?)")

    args = parser.parse_args()
    print("Loading summary statistics and relevant genotype information....")
    stats, header = loadSummaryStats(args.sum_stats, args.ss_format)
    variants = loadGenoVars(args.plink_snps + ".pvar")
    print("Filtering SNPs...")
    snp_list, stats_filtered, variants_filtered = filterSNPs(args.pval, args.maf, stats, variants, header) #TODO: makes sure the variatns are in order.

    #Run plink filtering on the bed file, and then load it into memory
    print("Building a PLINK matrix...")
    snp_matrix = plinkToMatrix(snp_list, args.plink_snps)
    #stats was output of filterSNPs
    print("Calculating PRSs")
    scores = calculatePRS(stats_filtered, snp_matrix)
    patient_ids  = getPatientIDs(snp_matrix)
    writeScores(scores, patient_ids, args.output)

    #Match SNPs by address and/or reference ID.
    print("scores written out to", args.output)
    #add option to accomodate large size files with limited memory- chunk it.

    #Load the summary stats into memory quickly.

    #Make any checks or corrections that are needed
    #Check for ambiguous SNP changes- they should be removed


    #MAke sure 

    #Add an option to do clumping if you want?

    #Make another script for easy visualization of the outpu.
    
    


