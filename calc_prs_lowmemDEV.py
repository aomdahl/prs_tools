"""
Ashton Omdahl
August 2019
"""
import argparse
import numpy as np
import scipy as sp
import pandas as pd
from datetime import date
import sys
from subprocess import call 
import os
import time
from plinkio import plinkfile
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
        #print("Filtering sum stats, see if we get the right indices")
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
        
        command = '''awk '(!/##/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' ''' + snp_file + ".pvar  > local_geno.pvar"
        #Simplifying things, we aren't going to re-write IDs....
        if id_type:
            print("ID type")
            command = '''awk '(!/##/) {print $1"\t"$2"\t"$3":"$4":"$5"\t"$4"\t"$5}' ''' + snp_file + ".pvar " + " | sed '1 s/ID:REF:ALT/ID/' > local_geno.pvar"
        #command = '''awk '(!/##|ID/ && $1 !~ /X/) {print $3":"$4":"$5}' ''' + snp_file + ".pvar > geno_ids.f"
        call(command, shell = True)
        #command_n = "tail -n +2 local_geno.pvar | awk ' ($1 !~ /X/) {print $1}' |   cut -f 3 > geno_ids.f"
        command_n = "awk ' (NR> 1 && $1 !~ /X/) {print $3}'  local_geno.pvar > geno_ids.f" #can do this in one awk command too you goon.
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

    return local_geno,"geno_ids.f", "ss_filt.f"

def scoreCalculation(geno, betas):
    return np.dot(2 - geno, betas)   

def hashMapBuild(ss_snps, plink_snps, ret_tap): #we need to deal with this issue of multiallelic snps
    pref = dict()
    for i in range(0, len(plink_snps)):
        id_name = str(plink_snps[i].name) + ":" + str(plink_snps[i].allele2) + ":" + str(plink_snps[i].allele1) #why is this order so? nervous making
        pref[id_name] = i
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
        try:
            id_name = str(plink_order[i].name) + ":" + str(plink_order[i].allele2) + ":" + str(plink_order[i].allele1)
            if curr_ss != plink_order[i].name and curr_ss != id_name: #But if for some reason it isn't a match....
                #TODO: check the alt and ref here potentially.
                ret = hashMapBuild(ss_order, plink_order, ret)
                break
        except IndexError:
            print("Something went wrong...")
            print(i, curr_ss)
            print(curr_ss, plink_order[i])
            input()
        
        ret[i] = int(map_index)
    return ret
    

def linearGenoParse(snp_matrix, snp_indices,pvals, scores, ss_ids):
    #C
    from operator import itemgetter
    plinkf = plinkfile.PlinkFile(snp_matrix)
    print("Genotype binary loaded.")
    patient_ids = plinkf.get_samples()
    num_patients = len(patient_ids)
    report_index = int(num_patients * 0.1)
    #Check that our mapping of loci is right....
    loci = plinkf.get_loci()
    var_map = buildVarMap(loci, ss_ids) #Need to check that this goes the way we expect.
    for i in range(0, num_patients):
        line = next(plinkf)
        #for line in istream:
        #time comparison
        """
        start_np = time.time()
        dat = np.array(line,dtype = np.float64)
        for p in pvals:
            snp_index = snp_indices[p][INDEX]
            sel = var_map[snp_index]
            scores[p].append(scoreCalculation(dat[sel],snp_indices[p][B])) 
        end_np = time.time()
        print("time for np method:", str(end_np - start_np))
        #__________ 
        #now doing it the pythonic way
        start_ig = time.time()
        """
        for p in pvals:
            snp_index = snp_indices[p][INDEX]
            sel = var_map[snp_index]
            res_list = np.array((itemgetter(*sel)(line)))  #sadly we still have this...
            scores[p].append(scoreCalculation(res_list,snp_indices[p][B]))
        #end_ig = time.time()
        #print("time for itemgettermethod:", str(end_ig - start_ig))



        if i % report_index == 0:
            print("Currently at patient", i + 1, "of", num_patients)

    return scores, patient_ids

def calculatePRS(pvals, snp_matrix, snp_indices, snp_list):

    #scores = dict() #where the scores get stored
    #for p in pvals:
    #    scores[p] = list()
    #Seraj shoutout
    scores = { p : [] for p in pvals}
    scores, patient_ids = linearGenoParse(snp_matrix, snp_indices,pvals, scores, snp_list)
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
    if os.path.isfile("mat_form_tmp.bed"):
        print("Using currently existing genotype matrix...")
    else:        
        #call(["plink2", "--pfile", args, "--pvar", local_pvar, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "A", "--max-alleles", "2"]) 
        #Okay, this new version changes it to plink 1, doesn't write out a massive matrix, and makes it patient first indexed
        call(["plink2", "--pfile", args, "--pvar", local_pvar, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "ind-major-bed", "--max-alleles", "2"])
        print("New plink matrix made!")
        if not os.path.isfile("mat_form_tmp.bed"): #C
            print("PLINK file was not correctly created. Program will terminate.")
            sys.exit()
    return "mat_form_tmp"

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
            write_s = str(ids[i].iid) + '\t'
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

    #So what we need to do differently:
        #Let's keep using plink to filter down the SNP list, but save it out as a bed/bim format file. 
        #assuming that works....
        #Load it into plinkio\
        #Load the transposed one into memory (or just transpose?) Not sure if this will do it.
        #Then return the plinkf object to iterate.
        #Iterate through it like the lines of the matrix file, except you aren't reading those in line by line. Hopefully that proves to be faster.
        #This method has slowdowns in transposing the data. If I could avoid doing that then it wouldn't be an issue.
        #I guess the boost of going this other directino- not loading everything into memory- is no longer an issue, as plinkio is efficient.
        #So, make one pass through the plinkio, 
        #Note- the key xommands to know are just that plinkfile.PlinkFile has an iterator
        #i.e. line = next(plinkf)
        #repeat until youve gone through all samples.


        #plinkf = plinkfile.PlinkFile("./ref") 
        # #okay, the transpose call is slower than I'd like, certainly slower than reading in. I think that's because it does a write out, I wonder if this can be overridden. 
        #Sad news- the transpose function doens't store anything, just writes it out to memory.
        #Can plink transpose the data for me? hmmmmm
        #Maybe. It looks like maybe. In that case, what if we use plink to subset our SNPs (although, this isn't the smartest way to do it), then transpose the data, then load it in and iterate it?
        #This is an option to pursue.
        #Alternatively, just go through only to the relevant SNPs andg et the info, but that seems wasteful.
        #So there is something rotated of the file

        #So it doesn't look like plink2 will do the transposition for me. That's weird
        #So is it better to transpose it use libplinkio or just marse it the way they have?
        #Or should I rebuild my whole parser?


    #snp_matrix_plinkio = plinkioReadIn(args.plink_snps)
    snp_indices = filterSumStats(pvals,stats_complete)
    #Now do it 
    print("Parsing the genotype file...")
    scores, patient_ids = calculatePRS(pvals, snp_matrix, snp_indices,stats_complete[REFID])
    writeScores(scores, patient_ids, args.output + ".tsv")
    #Match SNPs by address and/or reference ID.
    print("Scores written out to", args.output + ".tsv")
    stop = time.time()
    print("Total runtime:", str(stop - start))



