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
import dask.dataframe
from subprocess import call 
import os

#General reference things
CHR = "chr"
POS = "pos"
REFID = "ref_id"
EFFAL = "effect_allele"
EFFALFREQ = "effect_allele_freq"
BETA = "beta"
SEBETA="sebeta"
TSTAT = "tstat"
PVAL = "pval"
REF= "REF" #reference allele
ALT = "ALT" #Alternate allele

def getStatHeaders(setting):
    if setting == "SAIGE":
        hl = "#CHROM POS ID REF ALT ac af num_cases num_controls beta sebeta Tstat pval pval_SAIGE_NoSPA Is_Converged varT varTstar".split(" ")
        return {CHR: hl[0], POS:hl[1], REFID: hl[2], EFFAL: hl[4], EFFALFREQ:hl[6], BETA:hl[9], PVAL:hl[12], ALT:hl[4], REF:hl[3]}
        #Can certainly add more later, but this is a fine start.
        
def readInCol(fin):
    ret_list = list()
    with open(fin, 'r') as istream:
        for line in istream:
            line = line.strip()
            ret_list.append(line)
    return ret_list

def readInSummaryStats(s_path, ids_path):
    ss = dask.dataframe.read_csv(s_path, sep = ' ', header = None,  names = [REFID, REF,ALT, BETA,SEBETA,TSTAT, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64', SEBETA:'float64', TSTAT:'float64'})
    print("Data in dask...")
    return ss, readInCol(ids_path)

        
def prepSummaryStats(geno_ids, sum_path, pval_thresh):
    """
    This function uses awk to filter the summary stats by pvalue, and then reads them into a dask data structure.
    """
    #If the id is in the geno_ids file, and it passes the p-value threshold, keep it
    command = "awk '(FNR == NR) {a[$1];next} ($1 in a && $7 <= " + str(pval_thresh) + ") {print $0}' " + geno_ids + " " + sum_path + " > ss.tmp"
    call(command, shell = True) 
    #call(["awk", "'FNR==NR{a[$1];next} ($1 in a && $7 <=", str(pval_thresh), "{ print $0}'", geno_ids, sum_path, ">", "ss.tmp"])
    command = "cut -f 1 -d ' ' ss.tmp > ss_ids.tmp"
    call(command, shell = True)
   # call(["cut", "-f", "1", "-d","' '", "ss.tmp", ">", "ss_ids.tmp"]) #This gives the order of all the ids.
    return "ss.tmp", "ss_ids.tmp"

#This will extract the ids we want to be working with.
def prepSNPIDs(snp_file, ss_file, ss_type, ambig):
    """
    This extracts all the IDs from the genotype data (pvar file) that we wish to be using and selects just the data from the summary stats data we want
    @return path to the genotype ids
    @return path to the summary stats.
    """
    if not ambig:
        print("Please remove ambiguous snps. program will terminate")
        sys.exit()

    #Remove all the header lines, just get what we need.
    if not os.path.isfile('geno_ids.f'):
        command = "awk '(!/##|ID/ && substr($3,1,2) !~ /X:/) {print $3 }' " + snp_file + ".pvar > geno_ids.f"
        call(command, shell = True)
    #call(["awk", "'!/##|ID/ {print $3 } '", (snp_file + ".pvar"), ">", "geno_ids.f"])
    if not os.path.isfile("ss_filt.f"):
        if ss_type == "SAIGE":
            command = '''awk '(FNR == 1) {next;} {print $1":"$2, $4,$5,$10,$11,$12,$13}' ''' + ss_file + " > ss_filt.f"
            call(command, shell = True)
            #call(["awk", """'(FNR == 1) {next;} {print $1":"$2, $4,$5,$10,$11,$12,$13}'""",ss_file, ">", "ss_filt.f"]) #TODO: find ways to speed up this step and downstream ones, maybe filter out X chrom?
        else:
            print("The type of ss file hasn't been specified. Please specify this.")

    return "geno_ids.f", "ss_filt.f"


def filterSNPs(pval, maf, ss, vars, h):
    #@param h is the header data.
    #Get list of intersecting SNPs that pass the p-value threshold. Can add the MAF threshold at a later time.
    ss_new, vars_new = getIntersectingData(ss[(ss[h[PVAL]] <= pval) & (ss[h[CHR]] != "X")], vars)
    return ss_new.index, ss_new, vars_new
    #return ss_new[REFID], ss_new, vars_new

def calculatePRS(ss, snp_matrix, snp_list):
    #check = ss.set_index(REFID)
    #The above step isn't really necessary either, everything should be in order?
    #assert(sameAlleleAssesment(check[ALT], snp_matrix[ALT]))
    #print("Same alt alleles!")
    #assert(sameAlleleAssesment(check[REF], snp_matrix[REF]))
    #print("Same REF alleles!")
    #get the ALT column of the SNP matrix and make sure it matches with the alleles for which we have effects.
    #okay, I need to look at this interactively >_<
    dat = (2-snp_matrix.iloc[:,6:])
    #dat = (2-snp_matrix.loc[:,snp_list]) #this should yield an n x m array, n patients by m snps
    betas = ss[BETA] #Select the    Betas, should be a list of length m
    scores = np.matmul(betas, dat)
    return scores

def sameAlleleAssesment(s1, s2):
    t = (s1 != s2).sum()
    if t == 0:
        return True
    else:
        print(s1, s2)
        print("^Unusual case")
        return False


#TODO: This is slow, probably because of the sorting. Come up with a way to maintain order and speed.
#Checked functionality theoretically, seems to work.
def getIntersectingData(ss, vars):
    """
    @param ss: the summary statistics file, in a filtered dask
    @param vars: the read in pvar file, all the SNPs for which we have information.
    """
    #Because ss is a dask structure, we should select IDs from vars instead.
    #Also note that the dask is currently indexed by ID, so this should be FAST
    #id_keep = ss[REFID].values
    #select out the values in vars that are in the vars_keep list (make sure IN ORDER)
    #Dask cannot sort values
    ss_filt = ss.loc[vars.index]
    #vars_keep = vars[vars[REFID].isin(id_keep)].sort_values(by=REFID) #Remember, if its a 3 we ignore it- it means the data isn't available for them.
    #return ss[ss[REFID].isin(vars_keep[REFID])].sort_values(by=REFID), vars_keep 
    return ss_filt, vars_keep
    
def plinkToMatrix(snp_keep, args):
    """
    @param snp_keep is the path to the file with the SNP list. Currently just a single column, TODO check if this will work for plink filtering.
    """
    #Write out the list of SNPs to keep
    snp_order = snp_keep
    #list_dir = writeSNPIDs(snp_keep)
    #Make the new matrix
    #call(["plink2", "--pfile", args, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "A","'include-alt'", "--max-alleles", "2"]) 
    
    call(["plink2", "--pfile", args, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "A-transpose", "--max-alleles", "2"]) 
    #print("Cleaning up header")
    #clean_header= "sed -i '1 s/\([1-22]:[0-9]*\)_[ATCGN](\/[ATGCN])/\1/g' mat_form_tmp.raw"
    #call(clean_header, shell = True)
    ret = dask.dataframe.read_csv("mat_form_tmp.traw", sep = "\t") #Changed this to a dask, let's see if its faster...
    #No need to set an index- everything is in order!
    #ret = ret.set_index("SNP")
    ret = ret.rename(columns={"COUNTED":REF})
    return ret

def getPatientIDs(snp_matrix):
    tnames = list(snp_matrix.columns)[6:]
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
        
    print("Selecting IDs for analysis...")
    geno_ids, ss_parse = prepSNPIDs(args.plink_snps, args.sum_stats,args.ss_format, args.no_ambig)
    for pval in pvals:
        print("For", str(pval), "threshold, selecting relevant SNP IDs") #The best thing to do up front would be to assign each SNP to a category based on the P-value thresholds it passes.... this owuld be much faster overall
        stats_file, ids_file = prepSummaryStats(geno_ids, ss_parse, pval) #Gives a list of ids in order and the summary stats
        print("Reading filtered data into memory")
        stats_filtered, snp_list = readInSummaryStats(stats_file,ids_file)
        #stats_filtered is our ss dask, snp_list is our ids in memory
        #snp_list, stats_filtered, variants_filtered = filterSNPs(pval, args.maf, stats, variants, header)
        #Run plink filtering on the bed file, and then load it into memory
        print("Building a PLINK matrix for pval", str(pval), "...")
        snp_matrix = plinkToMatrix(ids_file, args.plink_snps) #Already written to file, so can use this...
        #stats was output of filterSNPs
        print("Calculating PRSs")
        #stats_qc, snp_qc = qualityControl(stats_filtered, snp_matrix)
        scores = calculatePRS(stats_filtered, snp_matrix, snp_list)
        print(scores)
        input()
        patient_ids  = getPatientIDs(snp_matrix)
        writeScores(scores, patient_ids, args.output + "_" + str(pval) + ".tsv")
        #Match SNPs by address and/or reference ID.
        print("Scores written out to", args.output + "_" + str(pval) + ".tsv")
        #add option to accomodate large size files with limited memory- chunk it.


