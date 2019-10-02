"""
Ashton Omdahl
August 2019
"""
#Make sure you are using the NON-PYMOL interpreter.
import argparse
import numpy as np
import scipy as sp
import pandas as pd
from datetime import date
import sys
from subprocess import call 
import os
import time
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
INDEX = 0
B = 1

        
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
    
    #return ss.sort_values(by=[REFID])#, readInCol(ids_path)
    return ss

def filterSumStats(pvals,ss):
    pvals_ref = dict() #dictionary containing a pointer for each p value
    for p in pvals:
        print("Filtering sum stats, see if we get the right indices")
        samp = ss[(ss[PVAL] <= float(p))]
        #print(samp)
        #print(samp.index.tolist())
        pvals_ref[p] = [samp.index.tolist(), samp[BETA].values] #Get just the pvalues and the indices.
        #print("Scores as extracted from sum stats....")
        #print(pvals_ref[p][0], pvals_ref[p][1])
    return pvals_ref
        
def prepSummaryStats(geno_ids, sum_path, pval_thresh):
    """
    This function uses awk to filter the summary stats by pvalue, and then reads them into a dask data structure.
    """
    #If the id is in the geno_ids file, and it passes the p-value threshold, keep it
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
def prepSNPIDs(snp_file, ss_file, ss_type, ambig, id_type):
    """
    This extracts all the IDs from the genotype data (pvar file) that we wish to be using and selects just the data from the summary stats data we want
    @return path to the genotype ids
    @return path to the summary stats.
    """
    local_geno = snp_file + ".pvar"
    if id_type:
        local_geno = "local_geno.pvar"
    if not ambig:
        print("Please remove ambiguous snps. program will terminate")
        sys.exit()

    #Remove all the header lines, just get what we need.
    if not os.path.isfile('geno_ids.f'):
        
        #command = '''awk '(!/##/ && $1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' ''' + snp_file + ".pvar  > local_geno.pvar"
        #Simplifying things, we aren't going to re-write IDs....
        if id_type:
            command = '''awk '(!/##/ && $1) {print $1"\t"$2"\t"$3":"$4":"$5"\t"$4"\t"$5}' ''' + snp_file + ".pvar " + " | sed '1 s/ID:REF:ALT/ID/' > local_geno.pvar"
        #command = '''awk '(!/##|ID/ && $1 !~ /X/) {print $3":"$4":"$5}' ''' + snp_file + ".pvar > geno_ids.f"
            call(command, shell = True)
        
        command_n = "tail -n +2 " + local_geno + " | cut -f 3 > geno_ids.f"
        call(command_n, shell = True)
    
    if not os.path.isfile("ss_filt.f"):
        if ss_type == "SAIGE":
            #ID, REF, ALT, BETA, SEBETA, Tstat,, pval
            #command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1":"$2":"$4":"$5 in a) {print $1":"$2":"$4":"$5, $4,$5,$10,$11,$12,$13}' geno_ids.f ''' + ss_file + " > ss_filt.f"
            command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1":"$2 in a) {print $1":"$2, $4,$5,$10,$11,$12,$13}' geno_ids.f ''' + ss_file + " > ss_filt.f"
        elif ss_type == "NEAL":
            print("Nealing it out....")
            command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a) {print $1, "N",$2,$8,$9,$10,$11}' geno_ids.f ''' + ss_file + " > ss_filt.f"   

        else:
            print("The type of ss file hasn't been specified. Please specify this.")
    

        #command = '''awk '(FNR == 1) {next;} {print $1":"$2, $4,$5,$10,$11,$12,$13}' ''' + ss_file + " > ss_filt.f"
        call(command, shell = True)
             
        #reset the ids we use downatream
        command = "cut -f 1 -d ' ' ss_filt.f > geno_ids.f"
        call(command, shell = True)
    
    #print("Did this work out right? Sample the files")
    #input()
    return local_geno,"geno_ids.f", "ss_filt.f"

def scoreCalculation(geno, betas):
    return np.dot(2- geno, betas)   

def hashMapBuild(ss_snps, plink_snps, ret_tap):
    pref = dict()
    for i in range(0, len(plink_snps)):
        pref[plink_snps[i][:-2]] = i
    for j in range(0, len(ss_snps)):
        curr_id = ss_snps.iloc[j]
        try:
            ret_tap[j] = int(pref[curr_id])
        except KeyError:
            print("Unable to find a match for ", curr_id)
            ret_tap[j] = -1
            #input()
    return ret_tap
        
def buildVarMap(plink_order, ss_order):
    """
        Return a numpy array that maps an index in the ss file to a number in plink
    """
    ret = np.zeros(len(ss_order), dtype = np.int32)
    search_count = 0
    for i in range(0, len(ss_order)):
        curr_ss = ss_order.iloc[i]
        map_index = i #assume they are ordered the same
        if curr_ss != plink_order[i][:-2]: #But if for some reason it isn't a match....
            #search_count += 1
            #map_index = plink_order.index(curr_ss)
            ret = hashMapBuild(ss_order, plink_order, ret)
            break
        ret[i] = int(map_index)
    return ret
    

def linearGenoParse(snp_matrix, snp_indices,pvals, scores, patient_ids, ss_ids):
    with open(snp_matrix, 'r') as istream:
        counter = 0
        line = istream.readline().strip()
        while line:
        #for line in istream:
            line = line.strip()
            if counter == 0:
                
                #line = line.replace("_", ":")
                var_map = buildVarMap(line.split()[6:], ss_ids)
                line = istream.readline().strip()
                counter += 1
                continue
            dat = np.array(line.split(), dtype = np.float64)
            patient_ids.append(dat[0])
            dat = dat[6:] #Just the numberic values at the end
            for p in pvals:
                snp_index = snp_indices[p][INDEX]
                sel = var_map[snp_index]
                scores[p].append(scoreCalculation(dat[sel],snp_indices[p][B]))
            counter += 1
            if counter % 100 == 0:
                print("Currently at patient", counter)
            line = istream.readline()
            #if counter % 10 == 0 :
            #    break
    return scores, patient_ids

def calculatePRS(pvals, snp_matrix, snp_indices, snp_list):

    #scores = dict() #where the scores get stored
    #for p in pvals:
    #    scores[p] = list()
    #Seraj shoutout
    scores = { p : [] for p in pvals}
    patient_ids = list() #store the patient's ids...
    scores, patient_ids = linearGenoParse(snp_matrix, snp_indices,pvals, scores, patient_ids, snp_list)
    return scores, patient_ids

def sameAlleleAssesment(s1, s2):
    t = (s1 != s2).sum()
    if t == 0:
        return True
    else:
        print(s1, s2)
        print("^Unusual case")
        return False

def plinkToMatrix(snp_keep, args, local_pvar):
    """
    @param snp_keep is the path to the file with the SNP list. Currently just a single column, TODO check if this will work for plink filtering.
    """
    #Write out the list of SNPs to keep
    if os.path.isfile("mat_form_tmp.raw"):
        print("Using currently existing genotype matrix...")
    else:
        #call(["plink2", "--pfile", args, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "A", "--max-alleles", "2"]) 
        
        call(["plink2", "--pfile", args, "--pvar", local_pvar, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "A", "--max-alleles", "2"]) 
        print("New plink matrix made!")
        #input()
        if not os.path.isfile("mat_form_tmp.raw"):
            print("PLINK file was not correctly created. Program will terminate.")
            sys.exit()
    return "mat_form_tmp.raw"

def getPatientIDs(snp_matrix):
    tnames = list(snp_matrix.columns)[6:]
    return [int(i.split("_")[0]) for i in tnames]


def writeScores(scores, ids, destination):
    with open(destination, 'w') as ostream:
        p_string = "IID\t"
        num = 0
        for p in scores:
            p_string = p_string + str(p) + '\t'
            num = len(scores[p]) #number of patients
        ostream.write(p_string[:-1] + '\n')
        for i in range(0, num):
            write_s = str(int(ids[i])) + '\t'
            for p in scores:
                write_s = write_s + str(scores[p][i]) + '\t'
            write_s = write_s[:-1] + '\n'
            ostream.write(write_s)

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
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE", "NEAL"])
    parser.add_argument("-p", "--pval", default = 5e-8, type = float, help = "Specify what pvalue threshold to use.")
    parser.add_argument("--pvals", help= "Use this if you wish to specify multiple at once, in a list separated by commas")
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--no_ambig", default = False, action = "store_true", help = "Specify this option if you have already filter for ambiguous SNPs and bi-alleleic variants. Recommended for speed.")
    parser.add_argument("--var_in_id", action = "store_true", help = "Specify this if you wish to use more specific id of the format chr:loc:ref:alt. For default SNP in file, just leave this.")
    args = parser.parse_args()
    pvals = [args.pval]

    start = time.time()
    if args.pvals:
        pvals = [float(x) for x in (args.pvals).split(",")]
        
    print("Selecting IDs for analysis...")
    local_pvar, geno_ids, ss_parse = prepSNPIDs(args.plink_snps, args.sum_stats,args.ss_format, args.no_ambig, args.var_in_id)
    print("Reading data into memory (only done on first pass)")
    stats_complete = readInSummaryStats(ss_parse)
    snp_matrix = plinkToMatrix(geno_ids, args.plink_snps, local_pvar) 
    snp_indices = filterSumStats(pvals,stats_complete)
    #Now do it 
    print("Parsing the genotype file...")
    scores, patient_ids = calculatePRS(pvals, snp_matrix, snp_indices,stats_complete[REFID])
    writeScores(scores, patient_ids, args.output + ".tsv")
    #Match SNPs by address and/or reference ID.
    print("Scores written out to", args.output + ".tsv")
    stop = time.time()
    print("Total runtime:", str(stop - start))



