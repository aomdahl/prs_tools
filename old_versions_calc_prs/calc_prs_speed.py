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

def readInSummaryStats(s_path):
    ss = pd.read_csv(s_path, sep = ' ', header = None,  names = [REFID, REF,ALT, BETA,SEBETA,TSTAT, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64', SEBETA:'float64', TSTAT:'float64'})
    print("Data in memory...")
    ss.info(memory_usage='deep') 
    
    return ss#, readInCol(ids_path)

def filterSumStats(pval,ss,):
    #@param h is the header data.
    #Get list of intersecting SNPs that pass the p-value threshold. Can add the MAF threshold at a later time.
    ss_new = ss[(ss[PVAL] <= pval)]
    return ss_new[REFID], ss_new
        
def prepSummaryStats(geno_ids, sum_path, pval_thresh):
    """
    This function uses awk to filter the summary stats by pvalue, and then reads them into a dask data structure.
    """
    #If the id is in the geno_ids file, and it passes the p-value threshold, keep it
    #With the most recent 9/13 update, this is superfluousi
    #In fact, we should just read this into a dask and call it a day
    command = "awk '($7 <= " + str(pval_thresh) + ") {print $0}' " + sum_path + " > ss.tmp"
    #command = '''awk '(FNR == NR) {a[$1];next} ($1":"$3 in a && $7 <= ''' + str(pval_thresh) + ") {print $0}' " + geno_ids + " " + sum_path + " > ss.tmp"
    #command = "awk '(FNR == NR) {a[$1];next} ($1 in a && $7 <= " + str(pval_thresh) + ") {print $0}' " + geno_ids + " " + sum_path + " > ss.tmp"
    call(command, shell = True) 
    #call(["awk", "'FNR==NR{a[$1];next} ($1 in a && $7 <=", str(pval_thresh), "{ print $0}'", geno_ids, sum_path, ">", "ss.tmp"])
    command = "cut -f 1 -d ' ' ss.tmp > ss_ids.tmp"
    call(command, shell = True)
    #TODO: determine if the above awk/cut commands are faster than dask filtering + computing the IDs. Or maybe just doing that with pandas.
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
        #command = "awk '(!/##|ID/ && substr($3,1,2) !~ /X:/) {print $3 }' " + snp_file + ".pvar > geno_ids.f"
        #Generate your own pvar file and use that, since the IDs are non-unique >_<
        #Still working on this (below_)
        
        command = '''awk '(!/##/ && $1 !~ /X/) {print $1"\t"$2"\t"$3":"$4":"$5"\t"$4"\t"$5}' ''' + snp_file + ".pvar " + " | sed '1 s/ID:REF:ALT/ID/' > local_geno.pvar"
        #command = '''awk '(!/##|ID/ && $1 !~ /X/) {print $3":"$4":"$5}' ''' + snp_file + ".pvar > geno_ids.f"
        call(command, shell = True)
        command = "tail -n +2 local_geno.pvar | cut -f 3 > geno_ids.f"
        call(command, shell = True)
    if not os.path.isfile("ss_filt.f"):
        if ss_type == "SAIGE":
            #ID, REF, ALT, BETA, SEBETA, Tstat,, pval
            command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1":"$2":"$4":"$5 in a) {print $1":"$2":"$4":"$5, $4,$5,$10,$11,$12,$13}' geno_ids.f ''' + ss_file + " > ss_filt.f"
            #command = '''awk '(FNR == 1) {next;} {print $1":"$2, $4,$5,$10,$11,$12,$13}' ''' + ss_file + " > ss_filt.f"
            call(command, shell = True)
            #call(["awk", """'(FNR == 1) {next;} {print $1":"$2, $4,$5,$10,$11,$12,$13}'""",ss_file, ">", "ss_filt.f"]) #TODO: find ways to speed up this step and downstream ones, maybe filter out X chrom?
        else:
            print("The type of ss file hasn't been specified. Please specify this.")

    return "local_geno.pvar","geno_ids.f", "ss_filt.f"



def calculatePRS(ss, snp_matrix):
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
    
def plinkToMatrix(snp_keep, args, local_pvar):
    """
    @param snp_keep is the path to the file with the SNP list. Currently just a single column, TODO check if this will work for plink filtering.
    """
    #Write out the list of SNPs to keep
    if os.path.isfile("mat_form_tmp.traw"):
        print("Using currently existing genotype matrix...")
    else:
        call(["plink2", "--pfile", args, "--pvar", local_pvar, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "A-transpose", "--max-alleles", "2"]) 
    #ret = dask.dataframe.read_csv("mat_form_tmp.traw", sep = "\t") #Changed this to a dask, let's see if its faster...
    #No need to set an index- everything is in order!
    import time
    start = time.time()
    ret = pd.read_csv("mat_form_tmp.traw", sep = '\t', index_col = "SNP")   
    end = time.time()
    print("read in time:", str(end-start))
    print("In memory size")
    ret.info(memory_usage='deep') 
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

def datasetsAlign(snps, ss, snp_list):
    """
    Goals is to make sure that everything checks out before calculating PRS. Specifically:
    Make sure the REF and ALT alleles are the same
    Make sure that we have the right dimensions.
    This shouldn't need to be run, but we can include it as needed...
    """
    if not sameAlleleAssesment(ss[ALT], snps[ALT]):
        print("The alternate alleles don't match as we'd expect, please review your data!")
        return False
    if not sameAlleleAssesment(ss[REF], snps[REF]):
        print("The reference alleles don't match as we'd expect, please review your data!")
        return False
    return True
    #print("Same alt alleles!")
    #assert(sameAlleleAssesment(check[REF], snp_matrix[REF]))
    #print("Same REF alleles!")
    

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
    local_pvar, geno_ids, ss_parse = prepSNPIDs(args.plink_snps, args.sum_stats,args.ss_format, args.no_ambig)
    print("Reading data into memory (only done on first pass)")
    stats_complete = readInSummaryStats(ss_parse)
    snp_matrix = plinkToMatrix(geno_ids, args.plink_snps, local_pvar) #Would a line-by-line approach be faster? #Or would a dask be better with the index set in?
    
    for pval in pvals:
        print("For", str(pval), "threshold, selecting relevant SNP IDs") #The best thing to do up front would be to assign each SNP to a category based on the P-value thresholds it passes.... this owuld be much faster overall
        #stats_file, ids_file = prepSummaryStats(geno_ids, ss_parse, pval) #Gives a list of ids in order and the summary stats
        fss_ids, filtered_ss = filterSumStats(pval,stats_complete)
        print("Subsetting the PLINK matrix for pval", str(pval), "...")
        #stats was output of filterSNPs
        geno_relevant = snp_matrix.loc[fss_ids] 
        print("Checking dimensions and allele references...")
        #assert(datasetsAlign(geno_relevant, filtered_ss, fss_ids))

        print("Calculating PRSs")
        #stats_qc, snp_qc = qualityControl(stats_filtered, snp_matrix)
        scores = calculatePRS(filtered_ss, geno_relevant)
        patient_ids  = getPatientIDs(geno_relevant)
        writeScores(scores, patient_ids, args.output + "_" + str(pval) + ".tsv")
        #Match SNPs by address and/or reference ID.
        print("Scores written out to", args.output + "_" + str(pval) + ".tsv")
        #add option to accomodate large size files with limited memory- chunk it.


