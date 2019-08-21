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

#General reference things
CHR = "chr"
POS = "pos"
REFID = "ref_id"
EFFAL = "effect_allele"
EFFALFREQ = "effect_allele_freq"
BETA = "beta"
PVAL = "pval"
REF= "ref_allele"
ALT = "alt_allele"

def getStatHeaders(setting):
    if setting == "SAIGE":
        hl = "#CHROM POS ID REF ALT ac af num_cases num_controls beta sebeta Tstat pval pval_SAIGE_NoSPA Is_Converged varT varTstar".split(" ")
        return {CHR: hl[0], POS:hl[1], REFID: hl[2], EFFAL: hl[4], EFFALFREQ:hl[6], BETA:hl[9], PVAL:hl[12], ALT:hl[4], REF:hl[3]}
        #Can certainly add more later, but this is a fine start.
        
        
        
def loadSummaryStats(sum_path, sum_format):
    h = getStatHeaders(sum_format)
    ss = pd.read_csv(sum_path, sep = ' ', index_col = False, 
                     dtype={h[CHR]:'category', h[POS]:'int32', h[ALT]: 'category', h[REF]:'category', h[BETA]:'float64', h[PVAL]:'float64', })
    ss[REFID] = ss[h[CHR]].astype(str) + ":" + ss[h[POS]].astype(str)
    return ss, h

def loadGenoVars(snp_file):
    #snp_file = plinkio.PlinkFile(snp_file)
    #put it into a format we can do something with... We actually really only need the bim file it looks like.  
    #bim_file = pd.read_csv(sum_path, sep = '\t', index_col = 1, header = None, Names = [CHR, REFID, "centi", POS, "allele1", "allele2"]) #1 is usually the minor
    #TODO: possible speedup, by choosing index column for quick lookup? (as above)
    #bim_file = pd.read_csv(sum_path, sep = '\t', header = None, names = [CHR, REFID, "centi", POS, "allele1", "allele2"], dtype={CHR:'int32', REFID:'str', 'centi':'int32', POS:'int32', "allele1": 'category', 'allele2':'category'}, skiprows=list(range(0,23))) 
    bim_file = pd.read_csv(snp_file, sep = '\t', header = None, names = [CHR, POS, REFID, REF, ALT], dtype={CHR:'category', REFID:'str', POS:'int32', ALT: 'category', REF:'category'}, skiprows = list(range(0,24))) 

    #Speed up the read in by specifying the type
    #Omit the first 23 lines- they are header information
    #Saves TONS of memory
    return bim_file

def filterSNPs(pval, maf, ss, vars, h):
    #@param h is the header data.
    #Get list of intersecting SNPs that pass the p-value threshold. Can add the MAF threshold at a later time.
    ss_new, vars_new = getIntersectingData(ss[(ss[h[PVAL]] <= pval) & (ss[h[CHR]] != "X")], vars)
    return ss_new[REFID], ss_new, vars_new


def calculatePRS(ss, vars,h, geno_mat, args):
    
    #do the necessary math: choose all the right rows.
    #Assuming the alt are the affect alleles.. that is what we assume. Was true in the small number of cases I looked at.
    dat = 2 - geno_mat.iloc[,6:] #Columns from 6 on
    betas = (ss[BETA].values).T
    scores = np.matmul(betas, dat)
    """
    #This current way of doing it is a bit janky, but gets the job done.
    effect_alleles = ss_k[EFFAL]
    #Alternate way
    ss_k["allele_match"]  = (effect_alleles == ss_k['allele1']).astype('int32') + (effect_alleles == ss_k['allele2']).astype('int32')

    #ss_k["allele_count"] = [0 if (a != e and b!= e) else 2 if (a == e and b == e) else 1 for a,b, e in zip(ss_k['allele1'], ss_k['allele2'], s_k[EFFAL])]   
    #I don't know about the speed of this command, it seems like it will be pretty slow.
    #USING PYPLINK
    """
    #return np.dot(ss[h[BETA]].values, vars_keep)
    return scores
    #wouldn't we want something scaled- as in values that have more SNPs?

#TODO: This is slow, probably because of the sorting. Come up with a way to maintain order and speed.
def getIntersectingData(ss, vars):
    id_keep = ss[REFID]
    vars_ids = vars[REFID].sort_values()
    #select out the values in vars that are in the vars_keep list (make sure IN ORDER)
    vars_keep = vars[vars[REFID].isin(id_keep)].sort_values(by=REFID) #Remember, if its a 3 we ignore it- it means the data isn't available for them.
    return ss[ss[REFID].isin(vars_keep[REFID])].sort_values(by=REFID), vars_keep, 
    

def plinkToMatrix(snp_keep, ss, vars, args):
    #Write out the list of SNPs to keep
    snp_order = list(snp_keep.values)
    list_dir = writeSNPIDs(snp_keep)
    #Make the new matrix
    from subprocess import call 
    call(["plink2", "--pfile", args.plink_snps, "--not-chr", "X", "--extract", list_dir, "--out", "mat_form_tmp", "--export", "A-transpose"]) 
    return pd.read_csv("mat_form_tmp.traw", sep = "\t")
    #we are gonna do it all with plink2. BOwho
    #Note that the order is as given by snp_order.
    #Run plink, write out as a matrix (now much smaller hopefully)
    #Call outsite of python.
    #Read in the matrix and associated files (ah... the order will be redone so we will need to anyway. Shoot. Wasted resources. Unless we can finalize resources.


def writeScores(s, destination):
    return

def writeSNPIDs(ids):
    ids.to_csv("snp_ids.tmp", header = False, index = False, sep = '\t') #(should have 2 columns of the family/infamily ids, which appear to be the same.)
    return "snp_ids.tmp"    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Basic tools for calculating rudimentary Polygenic Risk Scores (PRSs). This uses PLINK bed/bim/fam files and GWAS summary stats as standard inputs. Please use PLINK 1.9")
    parser.add_argument("-snps", "--plink_snps", help = "Plink file handle for actual SNP data (omit the .extension portion)")
    parser.add_argument("-ss", "--sum_stats", help = "Path to summary stats file. Be sure to specify format if its not DEFAULT. Assumes tab delimited")
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE"])
    parser.add_argument("-p", "--pval", default = 5e-8, type = float, help = "Specify what pvalue threshold to use.")
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--train", action = "store_true", help = "Select this if training a PRS model. This will save a list of the SNPs we include in the test, and can also draw some plots (?)")

    args = parser.parse_args()
    stats, header = loadSummaryStats(args.sum_stats, args.ss_format)
    variants = loadGenoVars(args.plink_snps)
    snp_list, stats_filtered, variants_filtered = filterSNPs(args.pval, args.maf, stats, variants, header) #TODO: makes sure the variatns are in order.

    #Run plink filtering on the bed file, and then load it into memory
    snp_matrix = plinkToMatrix(snp_list, stats_filtered, variants_filtered)
    #stats was output of filterSNPs
    scores = calculatePRS(stats_filtered, header, snp_matrix)
    writeScores(scores, args.output)

    #Match SNPs by address and/or reference ID.
    print("Success?")
    print(stats)

    input()
    #add option to accomodate large size files with limited memory- chunk it.

    #Load the summary stats into memory quickly.

    #Make any checks or corrections that are needed
    #Check for ambiguous SNP changes- they should be removed


    #MAke sure 

    #Add an option to do clumping if you want?

    #Make another script for easy visualization of the outpu.
    
    


