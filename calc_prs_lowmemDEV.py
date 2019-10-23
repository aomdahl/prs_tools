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
VAR_IN_ID = False
        
def readInCol(fin):
    ret_list = list()
    with open(fin, 'r') as istream:
        for line in istream:
            line = line.strip()
            ret_list.append(line)
    return ret_list

def readInSummaryStats(s_path):
    #ID, ref, alt, beta, pvalue. Th
    ss = pd.read_csv(s_path, sep = ' ', header = None,  names = [REFID, REF,ALT, BETA, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64', SEBETA:'float64', TSTAT:'float64'})
    #ss = pd.read_csv(s_path, sep = ' ', header = None,  names = [REFID, REF,ALT, BETA,SEBETA,TSTAT, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64', SEBETA:'float64', TSTAT:'float64'})
    print("Data in memory...")
    #ss.info(memory_usage='deep') 
    
    #return ss.sort_values(by=[REFID])#, readInCol(ids_path)
    return ss

def filterSumStats(pvals,ss):
    pvals_ref = dict() #dictionary containing a pointer for each p value
    for p in pvals:
        #print("Filtering sum stats, see if we get the right indices")
        samp = ss[(ss[PVAL] <= float(p))]
        #print(len(samp),"<- at threshold of", p)
        pvals_ref[p] = [samp.index.tolist(), samp[BETA].values] #Get just the pvalues and the indices.
        #print(samp[BETA].values[-5:-1])
        #print(samp.index.tolist()[-5:-1])
    #input()
    return pvals_ref
        
def prepSummaryStats(geno_ids, sum_path, pval_thresh):
    """
    This function uses awk to filter the summary stats by pvalue, and then reads them into a dask data structure.
    """
    #If the id is in the geno_ids file, and it passes the p-value threshold, keep it
    command = "awk '($7 <= " + str(pval_thresh) + ") {print $0}' " + sum_path + " > ss.tmp"
    call(command, shell = True) 
    command = "cut -f 1 -d ' ' ss.tmp > ss_ids.tmp"
    call(command, shell = True)
    #TODO: determine if the above awk/cut commands are faster than dask filtering + computing the IDs. Or maybe just doing that with pandas.
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
            command = '''awk '(!/##/) {print $1"\t"$2"\t"$3":"$4":"$5"\t"$4"\t"$5}' ''' + snp_file + ".pvar " + " | sed '1 s/ID:REF:ALT/ID/' > local_geno.pvar"
        #command = '''awk '(!/##|ID/ && $1 !~ /X/) {print $3":"$4":"$5}' ''' + snp_file + ".pvar > geno_ids.f"
        call(command, shell = True)
        #command_n = "tail -n +2 local_geno.pvar | awk ' ($1 !~ /X/) {print $1}' |   cut -f 3 > geno_ids.f"
        command_n = "awk ' (NR> 1 && $1 !~ /X/) {print $3}'  local_geno.pvar > geno_ids.f"
        call(command_n, shell = True)
        #input("blah") 
    if not os.path.isfile("ss_filt.f"):
        #ss = pd.read_csv(s_path, sep = ' ', header = None,  names = [REFID, REF,ALT, BETA,SEBETA,TSTAT, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64', SEBETA:'float64', TSTAT:'float64'})
        #Note that for all types, we expect IDs to be of the form chr:#:ref:alt
        #CHANGE on 10-22: only outputing the pare minimum, ID, ref, alt, beta, pvalue. The rest isn't needed
        if ss_type == "SAIGE":
            #ID, REF, ALT, BETA, SEBETA, Tstat,, pval
             command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($3 in a) {print $3, $4,$5,$10, $13}' geno_ids.f ''' + ss_file + " > ss_filt.f"  
            #command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($3 in a) {print $3, $4,$5,$10,$11,$12,$13}' geno_ids.f ''' + ss_file + " > ss_filt.f"  
        elif ss_type == "NEAL":
            command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a) {print $1, "N",$2,$8,$11}' geno_ids.f ''' + ss_file + " > ss_filt.f"    
            #command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a) {print $1, "N",$2,$8,$9,$10,$11}' geno_ids.f ''' + ss_file + " > ss_filt.f"   

        else:
            print("The type of ss file hasn't been specified. Please specify this.")
    
        #command = '''awk '(FNR == 1) {next;} {print $1":"$2, $4,$5,$10,$11,$12,$13}' ''' + ss_file + " > ss_filt.f"
        call(command, shell = True)
        #reset the ids we use downstream
        command = "cut -f 1 -d ' ' ss_filt.f > geno_ids.f"
        call(command, shell = True)
    return local_geno,"geno_ids.f", "ss_filt.f"

def scoreCalculation(geno, betas):
    return np.dot(2 - geno, betas)   

#There must be a better way than this...
def hashMapBuild(ss_snps, plink_snps, ret_tap): #we need to deal with this issue of multiallelic snps
    pref = dict()
    
    #Put all of the plink genotype SNPs into a dictionary for fast lookup
    for i in range(0, len(plink_snps)):
        if VAR_IN_ID:
            id_name= str(plink_snps[i].name)
            #print(id_name)
        else:
            id_name = str(plink_snps[i].name) + ":" + str(plink_snps[i].allele2) + ":" + str(plink_snps[i].allele1) #why is this order so? nervous making
        pref[id_name] = i
    ss_snp_names = ss_snps.values
    for j in range(0, len(ss_snps)):
        curr_id = ss_snp_names[j]
        try:
            ret_tap[j] = int(pref[curr_id])
        except KeyError:
            print("Unable to find a match for ", curr_id)
            ret_tap[j] = -1
    print("Re-indexing complete")
    return ret_tap
        
#Method for determine if IDs are in the right order
#Doesn't check every id, just samples a few to improve our confidence that we have the same order
#By checking >= 5 locations at random and ensuring they are the same. 
#If each position is equally likely for each entry, then the probability of having different entries by chance
#For this many checks is really small (i.e. (1/number of snps)^5).
def quickIDCheck(ss_order, plink_order):
    #Try a number of samples to you get a high certainty
    PROB_ACC = 0.0000001
    import random
    min_run = (np.log(PROB_ACC) / np.log(1/float(len(ss_order))))
    for i in range(0, int(min_run)+5): #I want to do at least 5. Actually a better way may be to sum up the ids across an interval and see if those are the same
        rand = random.randint(0, len(ss_order))
        if plink_order[rand].name != ss_order.iloc[rand]:
            #print("It didn't work as we wanted.....", plink_order[rand].name, ss_order.iloc[rand])
            print("SNP order not aligned, realigning now...")
            print(plink_order[rand].name)
            print(ss_order.iloc[rand])
            
            return False
        
    return True

def buildVarMap(plink_order, ss_order):
    """
        Return a numpy array that maps an index in the ss file to a number in plink
        TODO: Speed up this step, its too slow
    """

    if quickIDCheck(ss_order, plink_order):
        return np.array(range(0, len(ss_order) + 1)) 
     
    ret = np.zeros(len(ss_order), dtype = np.int32)
    ret = hashMapBuild(ss_order, plink_order, ret) 
    search_count = 0
    """
    #There is a better way:
    for i in range(0, len(ss_order)):
        curr_ss = ss_order.iloc[i]
        map_index = i #assume they are ordered the same
        try:
            id_name = str(plink_order[i].name) + ":" + str(plink_order[i].allele2) + ":" + str(plink_order[i].allele1)
            if curr_ss != plink_order[i].name and curr_ss != id_name: #But if for some reason it isn't a match....
                ret = hashMapBuild(ss_order, plink_order, ret)
                break
        except IndexError:
            print("Something went wrong...")
            print(i, curr_ss)
            print(curr_ss, plink_order[i])
            input()
        
        ret[i] = int(map_index)
    """
    return ret
    

def linearGenoParse(snp_matrix, snp_indices,pvals, scores, ss_ids):
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
        for p in pvals:
            snp_index = snp_indices[p][INDEX]
            sel = var_map[snp_index]  
            res_list = np.array((itemgetter(*sel)(line)))
            scores[p].append(scoreCalculation(res_list,snp_indices[p][B]))

        if i % report_index == 0:
            print("Currently at patient", i + 1, "of", num_patients)

    return scores, patient_ids

def calculatePRS(pvals, snp_matrix, snp_indices, snp_list):
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

#TODO- you can do this once with the 1000 genomes cohort and then just repeat with that list of "best" genes.
#Doesn't need to be done denovo every time!
def plinkClump(reference_ld, clump_ref, geno_ids, ss):
    """
    Run plink and clump, get the new list of ids.
    Select just these from the summary statistics.
    """
    if clump_ref == "NA":
        #Build the file for plink to run on 
        command = "echo SNP P > clump_ids.tmp"
        call(command, shell = True) 
        command = "awk ' {print $1,$NF} ' " + ss + ">> clump_ids.tmp"
        call(command, shell = True)
        #Run plink clumping
        plink_command = "plink --bfile " + reference_ld + "--clump clump_ids.tmp --clump-best"
        call(plink_command, shell = True)
        clump_file = "plink.clumped.best"
    else:
        print("Using provided reference clump...")
        clump_file = clump_ref
    #Refilter summary stats....
    #Actually, I should be able to just filter the list based on the plink input... buttt then we still have an order issue. Which I guess doesn't matter but.
    command = "awk '(FNR == NR) {a[$1];next} ($1 in a) {print $0}' " + clump_ref + " " + ss + " > t && mv t " + ss
    call(command, shell = True)
    #command = "awk '(FNR == NR) {a[$1];next} ($1 in a) {print $1}' " + clump_ref + " " + geno_ids + " > t && mv t " + geno_ids
    command = "cut -f 1 -d ' ' " + ss " > " geno_ids
    call(command, shell = True)
    return ss, geno_ids
    

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
    parser.add_argument("--clump", help = "Specify if you wish to do clumping. DEfault threshold is 0.5", action = "store_true", required = "--ld_ref" in sys.argv or "--clump_ref" in sys.argv)
    parser.add_argument("--ld_ref", help = "Specify an LD reference panel, required if you want to do clumping. Pass in plink1.9 format please.")
    parser.add_argument("--clump_ref", help = "Specify this option if you wish to perform clumping and already have a plink clump file produced (the output of plink1.9 with the --clump-best argument). Provide the path to this file.", default = "NA")
    args = parser.parse_args()
    pvals = [args.pval]

    start = time.time()
    if args.pvals:
        pvals = [float(x) for x in (args.pvals).split(",")]
    if args.var_in_id:
        VAR_IN_ID = True    
    
    print("Selecting IDs for analysis...")
    local_pvar, geno_ids, ss_parse = prepSNPIDs(args.plink_snps, args.sum_stats,args.ss_format, args.no_ambig, args.var_in_id)
    print("Reading data into memory (only done on first pass)")
    #Selecting out SNPs that are likely in LD
    #at this point, the geno ids and ss parse list  has what we want. We will reduce it more though through clumping via plink.
    if args.clump:
        ss_parse, geno_ids = plinkClump(args.ld_ref, args.clump_ref, geno_ids, ss_parse)
    stats_complete = readInSummaryStats(ss_parse)
    

    snp_matrix = plinkToMatrix(geno_ids, args.plink_snps, local_pvar) 

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



