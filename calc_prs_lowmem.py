"""
Ashton Omdahl
August 2019
"""
import os
import argparse
import numpy as np
import scipy as sp
import pandas as pd
from datetime import date
import sys
import subprocess
from subprocess import check_call 
import os
import time
import multiprocessing as mp
from operator import itemgetter
#General reference things
CHR = "chr"
POS = "pos"
REFID = "ref_id"
SNP = "snp"
BETA = "beta"
SEBETA="sebeta"
PVAL = "pval"
REF= "REF" #reference allele
ALT = "ALT" #Alternate allele
INDEX = 0
B = 1
global VAR_IN_ID
VAR_IN_ID = True #Variants are included in ID, assumed by default
DEBUG = True        
DEBUG_DAT = dict()
LOG="log_file.txt"

def updateLog(*prints):
    with open(LOG, 'a') as ostream:
        for a in prints:
            ostream.write(a + " ")
        ostream.write("\n")

def readInSummaryStats(s_path):
    #ID, ref, alt, beta, pvalue. Th
    #BEWARE: may have issue with sep, depending..
    ss = pd.read_csv(s_path, sep = '\t', header = None,  names = [REFID, REF,ALT, BETA, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64'})
    print("Data in memory...")
    return ss

def indexSSByPval(pvals,ss):
    """
    Generates an index that lists the location of pvalues between a given range
    i.e. betas < p= 0.001 [3,2,5,11]
    i.e. betas from 0.001 < p < 0.01 are at [12,14,15,161]
    i.e. betas from 0.01 < p < 1 are at [1,4,17,18...]
    listed relative to the next smallest SNP, for the way we calculate downstream 
    """
    pvals_ref = dict() #dictionary containing a pointer for each p value
    pvals.sort() #smallest one first
    """
    for p in pvals:
        samp = ss[(ss[PVAL] <= float(p))]
        pvals_ref[p] = [samp.index.values.tolist(), samp[BETA].values]      
    """
    prev = -1
    for p in pvals:
        samp = ss[(ss[PVAL] <= float(p)) & (ss[PVAL] > prev)]
        pvals_ref[p] = [samp.index.tolist(), samp[BETA].values]
        prev = float(p) 
    
    return pvals_ref
        
def prepSummaryStats(geno_ids, sum_path, pval_thresh):
    """
    This function uses awk to filter the summary stats by pvalue, and then reads them into a dask data structure.
    """
    #If the id is in the geno_ids file, and it passes the p-value threshold, keep it
    command = "awk '($7 + 0.0 <= " + str(pval_thresh) + ") {print $0}' " + sum_path + " > ss.tmp"
    check_call(command, shell = True) 
    #Pro tip- in benchmarks, cut performs better here than awk
    command = "cut -f 1 -d ' ' ss.tmp > ss_ids.tmp"
    check_call(command, shell = True)
    #TODO: determine if the above awk/cut commands are faster than dask filtering + computing the IDs. Or maybe just doing that with pandas.
    return "ss.tmp", "ss_ids.tmp"

#This will extract the ids we want to be working with.
#(args.plink_snps, args.sum_stats,args.ss_format, vplink = args.plink_version, max_p = max(pvals))
def prepSNPIDs(snp_file, ss_file, ss_type, vplink = 2, max_p = 1):
    """
    This extracts all the IDs from the genotype data (pvar file) that we wish to be using and selects just the data from the summary stats data we want
    @return path to the genotype ids
    @return path to the summary stats.
    """
    global VAR_IN_ID
    max_p = str(max_p)
    #Report the number of snps in the genotype file for log file
    if vplink == 2:
        local_geno = snp_file + ".pvar"
        updateLog("Number of variants detected in original genotype file", str(getEntryCount(snp_file + ".pvar", True)))
    else: #vplink == 1:
        local_geno = snp_file + ".bim"
        updateLog("Number of variants detected in original genotype file", str(getEntryCount(local_geno, False)))
    geno_id_list = "geno_ids.f"
    inter_sum_stats = "ss_filt.f" #Intersected summary stats with genotype ids
    #If we have already generated our geno_id_list,skip this tep 
    if not os.path.isfile(geno_id_list) or not VAR_IN_ID: #If we've already done this step, no point in doing more work....   
        if not VAR_IN_ID: #need to rework the IDs so they are chr:pos:ref:alt
            local_geno = "local_geno.pvar" #We will be making our own pvar
            command = ''' awk '(!/##/ && /#/) {print $1"\t"$2"\tID\t"$4"\t"$5} (!/#/) {print $1"\t"$2"\t"$1":"$2":"$4":"$5"\t"$4"\t"$5}' ''' + snp_file + ".pvar > " + local_geno
            if vplink == 1:
                local_geno = "local_geno.bim"
                command = '''awk '{print $1"\t"$1":"$4":"$5":"$6"\t"$3"\t"$4"\t"$5"\t"$6}' ''' + snp_file + ".bim > " + local_geno
            try:
                check_call(command, shell = True)
            except subprocess.CalledProcessError:
                print("Unable to generate modified pvar/bim file. Are you sure you specified the correct plink version?")
                print("Failed on command", command)
                sys.exit()
        if vplink == 2:
            command_n = "awk ' (NR> 1 && $1 !~ /X/) {print $3}' " + local_geno + " > " + geno_id_list
        else:
            command_n = "awk ' ($1 !~ /X/) {print $2}' " + local_geno + " > " + geno_id_list
        check_call(command_n, shell = True)
    
    VAR_IN_ID = True #we have imposed the variant information into the IDs
    if not os.path.isfile(inter_sum_stats):
        write_out = geno_id_list + " " + ss_file + " > " + inter_sum_stats
        if ss_type == "SAIGE":
            #ID, REF, ALT, BETA, SEBETA, Tstat, pval
             command = '''awk '(FNR==NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($3 in a && $13 + 0.0 <= ''' + max_p + " ) {print $3, $4,$5,$10, $13}' " + write_out 
        elif ss_type == "NEAL":
            command = '''awk '(FNR==NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a && $11 + 0.0 <= ''' + max_p + ''' ) {print $1, "N",$2,$8,$11}' ''' + write_out 
        elif ss_type == "HELIX": #SNP chr BP, A1 A2 P Beta --> ID ref alt beta pvalue
            command = '''awk '(FNR==NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1":"$5":"$4 in a && $6 + 0.0 <= ''' + max_p + ''' ) {print $1":"$5":"$4, $5,$4,$7,$6 }' ''' + write_out
        elif ss_type == "DEFAULT": #The kind of output we like, with a header ID ALT REF BETA SNP
            command = '''awk '(FNR==NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a && $5 + 0.0 <= ''' + max_p + ") {print $0}' " + write_out
        else:
            print("The type of ss file hasn't been specified. Please specify this.")
            sys.exit()
        check_call(command, shell = True)
        #reset the ids we use downstream
        command = " awk '{print $1}' " + inter_sum_stats + " > " + geno_id_list
        check_call(command, shell = True)

    return local_geno, geno_id_list, inter_sum_stats

def getEntryCount(fin, header):
    if header:
        return sum(1 for line in open(fin)) -1
    else:
        return sum(1 for line in open(fin))

def getPatientCount(snp_file, vers):
    try:
        if vers == 1:
            return getEntryCount(snp_file + ".fam", False)
        return getEntryCount(snp_file + ".psam", True)
    except:
        print("It appears the plink family file is missing or has the incorrect file extension. Please check it out!")
        sys.exit()

def scoreCalculation(geno, betas):
    #Note that what is done depends if the genotype data codes for the major or minor allele
    try: 
        m = 2-geno
        if max(m) > 2 or min(m) < 0: #bad genotype data we missed.
            print("Missing data detected. Attempting to correct appropriately")
            removal_indices = np.where(m >2)
            removal_indices.extend(np.where(m < 0))
            geno = np.delete(geno, removal_indices)
            betas_n = np.delete(betas, removal_indices)
            return np.dot(m, betas_n)
             
        return np.dot(m, betas)
    except ValueError:
        print("It appears that the number of SNPs is not as expected. We suggest lowering the missing genotype filter threshold, and additional debugging.")   

def hashMapBuild(ss_snps, plink_snps, ret_tap):
    """
    Create a mapping of SNPs in the SS data to the corresponding SNP in the genotype data.
    """
    pref = dict()
    import re
    #Put all of the plink genotype SNPs into a dictionary for fast lookup
    for i in range(0, len(plink_snps)):
        if VAR_IN_ID:
            new_key = re.sub("_[ACGTN]+$", "",plink_snps[i])
            pref[new_key] = i 
            #pref[plink_snps[i][:-2]] = i
        else:
            print("This portion of code is in development. Please finish, ashton")
            print("Variants have not be included in ids. Please set the genotype variant ids to be of the form chr:pos:ref:alt, or let this program do it automatically by specifing the --reset_ids tag")          
    ss_snp_names = ss_snps.values
    for j in range(0, len(ss_snps)):
        curr_id = ss_snp_names[j]
        try:
            ret_tap[j] = int(pref[curr_id])
        except KeyError:
            print("Unable to find a match for ", curr_id)
            updateLog(str(pref)) 
            ret_tap[j] = -1
            updateLog("Found unmatched ID, ", curr_id, "Please investigate.")
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
        if plink_order[rand][:-2] != ss_order.iloc[rand]:
            print("SNP order not aligned, realigning now...")
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
    return ret
    
def linearGenoParse(snp_matrix, snp_indices,pvals, scores, ss_ids):
    """Parse the patient file
    snp_matrix -- the path to the massive genotype file
    snp_indieces -- the map with the indices for each pvalue we care about
    scores -- what is returned, stores the scores
    pvals -- which pvalues we measure at
    ss_ids -- the list of ids from the summary stats list, for ordering purposes
    """
    patient_ids = list()
    report_index = 100
    pvals.sort()
    with open(snp_matrix + ".raw", 'r') as istream:
        line_counter = 0
        dat = istream.readline().strip().split('\t')
        line_len = len(dat)
        while dat:
            if line_counter == 0:
                var_map = buildVarMap(dat[6:], ss_ids)
                if DEBUG:
                    DEBUG_DAT[REFID] = dat[6:]
                    #print(len(DEBUG_DAT[REFID])) 
                print("Proceeding with line-by-line calculation now...")
                first_pval = True
                prev = 0
                for p in pvals:
                    if first_pval:
                        updateLog("Number of final SNPs used at ", str(p), "pvalue level:", str(len(snp_indices[p][INDEX]))) 
                        first_pval = False
                    else:
                        updateLog("Number of final SNPs used at ", str(p), "pvalue level:", str(len(snp_indices[p][INDEX]) + prev))
                    prev = prev + len(snp_indices[p][INDEX])
                if DEBUG:
                    print("Writing out requested debug information")
                    DEBUG_DAT[BETA] = dict()
                    DEBUG_DAT[SNP] = dict()
                    prev_append_snp = []
                    prev_append_beta = []
                    for p in pvals:
                        snp_index = snp_indices[p][INDEX]
                        sel = var_map[snp_index] 
                        DEBUG_DAT[BETA][p] = snp_indices[p][B]
                        DEBUG_DAT[SNP][p] = ((itemgetter(*sel)(dat[6:])))
                        DEBUG_DAT[BETA][p] = np.concatenate((DEBUG_DAT[BETA][p],prev_append_beta))
                        DEBUG_DAT[SNP][p] = np.concatenate((DEBUG_DAT[SNP][p], prev_append_snp))
                        prev_append_snp = DEBUG_DAT[SNP][p]
                        prev_append_beta = DEBUG_DAT[BETA][p]
                    writeScoresDebug(DEBUG_DAT[SNP], "debug_vals_snps.tsv")
                    writeScoresDebug(DEBUG_DAT[BETA], "debug_vals_betas.tsv")
            else:
                patient_ids.append(dat[0])
                
                first_pval = True
                rel_data = dat[6:]
                prev_score = 0
                new_score = 0
                for p in pvals:
                    snp_index = snp_indices[p][INDEX]
                    #sel = var_map[snp_index]#should also use item getter here
                    sel = (itemgetter(*snp_index)(var_map))
                    new_genos = np.array(itemgetter(*sel)(rel_data), dtype = np.float64)
                    if first_pval:
                        new_score = scoreCalculation(new_genos, snp_indices[p][B])
                        scores[p].append(new_score)
                        first_pval = False
                    else:
                        new_score = scoreCalculation(new_genos, snp_indices[p][B]) + prev_score
                        scores[p].append(new_score)
                        if DEBUG:
                            updateLog(str(patient_ids[-1]), str(new_score))
                    prev_score = new_score     
            try:   
                dat = istream.readline().split('\t')
            except: #end of the file
                dat = None
            line_counter += 1

            if line_counter % report_index == 0:
                print("Currently at individual/sample", line_counter)
            if len(dat) <= 1:
                print("Finished all the samples!")
                break
            
    return scores, patient_ids, DEBUG_DAT #basically a pval:[list of snps]  
    if line_counter == 1:
        print("There appears to be some error with the genotype plink matrix. Please ensure the genotype data is correct or the plink call was not interrupted.")
        sys.exit()
    return scores, patient_ids, None


def calculatePRS(pvals, snp_matrix, snp_indices, snp_list):
    """
    Actually calculates the PRS Scores, and times the process
    Returns a list of scores and the corresponding patient_ids, as well as any relevant debug informtion if specified
    @param snp_matrix- the path to the enormous snp matrix
    @param snp_indices: the information of snp_indices at each p-value
    @param snp_list: a series object containing the name of each snp
    """
    #Seraj shoutout
    start = time.time()
    scores = { p : [] for p in pvals}
    scores, patient_ids, debug_tab = linearGenoParse(snp_matrix, snp_indices,pvals, scores, snp_list)
    end = time.time()
    print("Calculation and read in time:", str(end-start))
    return scores, patient_ids, debug_tab

def alignReferenceByPlink(old_plink,plink_version, ss, ss_type):
    """
    Make the reference file listing <IDs> <REF>
    Call plink to swap them
    Write this out to a new_plink file
    """
    print("Correcting for strand flips")
    if ss_type == "SAIGE":
        first = '''awk '(NR >1) {print $1":"$2"\t"$4} ' ''' + ss + " | uniq > aligner.t"
        check_call(first, shell = True)
    elif ss_type == "NEAL":
        #Separate the major allele from the other info
        first = "cut -f 1 " + ss + " | tail -n +2 | cut -f 1,2 -d ':' > ids.t"
        second = "cut -f 1 " + ss + " | tail -n +2 | cut -f 3 -d ':' > ref.t"
        final = "paste ids ref | uniq> aligner.t"
        check_call(first, shell = True)
        check_call(second, shell = True)
        check_call(final, shell = True)
    elif ss_type == "HELIX":
        #ID CHR BP A1(effect) A2 P BETA
        quick_call = "awk ' ( NR > 1) {print $1,$5}' " + ss + " | uniq > ./aligner.t"
        check_call(quick_call, shell = True) 
    elif ss_type == "DEFAULT": #ID REf ALT BETA Pvalue
        call = "cut -f 1 " + ss + ''' | awk -F ":" '(NR > 1) {print $1":"$2"\t"$3}' | uniq > ./aligner.t'''
        #call = "awk ' (NR > 1) {print $1,$2}' ./" + ss + " | uniq > ./aligner.t"
        check_call(call, shell = True)
    else:
        print("Have not extended the method for type", ss_type, ". Please align Reference information manually")
        sys.exit()
    ftype = " --pfile "
    ext = ".pvar"
    remove_multi = "cut -f 3 " + old_plink + ext + " | uniq -d > remove_multi.t"
    if plink_version == 1:
        ftype = " --bfile "
        ext = ".bim"
        remove_multi = "cut -f 2 " + old_plink + ext + " | uniq -d > remove_multi.t"
    check_call(remove_multi, shell = True)
    #plink_call = "plink2" + ftype + old_plink  + " --not-chr X --out new_plink --ref-allele aligner.t --make-pgen --max-alleles 2" 
    #The above fails to call it correctly, unsurprisingly
    try:
        #12/12 adding an additional call to filter stuff
        #Filter out multi-allelic snps and X
        
        plink_call_first = "plink2" + ftype + old_plink + " --make-pgen --out filtered_target --exclude remove_multi.t --not-chr X"
        check_call(plink_call_first, shell = True)
        plink_call = "plink2 --pfile filtered_target --ref-allele force aligner.t --make-pgen --out new_plink"
        #plink_call = "plink2" + ftype + old_plink + " --ref-allele aligner.t --make-pgen --out new_plink"
        check_call(plink_call, shell = True)
    except:
        print("It appears there was an error in aligning the reference, likely due to known alleles that we attempted to swap. We suggest correcting these alleles or aligining the reference strands indepndently using plink --ref-allele")
        sys.exit()
    if not DEBUG:
        clean_up = "rm aligner.t"
        check_call(clean_up, shell = True)
    return "new_plink"

def filterSumStatSNPs(ss, ss_type):
    """
    Filters summary statistics based on the specified file type, as follows:
    1) Remove anything more than bi-allelic
    2) Limit it to SNPs only
    3) Limit to nonambiguous SNPs, if possible keep ambiguous snps MAF > 0.4
    4) Remove null beta values
    5) Remove X chromosome snps
    6) Write it out in the default format so its easier to digest downstream. (ID REF ALT BETA PVAL), ID standardized to chr:pos:ref:alt
    """
    filter_list = {"AG", "AC", "AT", "TA", "TC", "TG", "GC", "GT", "GA", "CG", "CT", "CA"}
    noambig_list = {"AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG"}
    ambig_list = {"AT", "TA", "GC", "CG"}
    sizes = {"start": 0,"non_bi": 0, "ambig":0, "na":0, "x":0}
    comment = {"ambig": "Ambiguous SNPs along with INDELS removed:", "na": "Variants with missing effect sizes removed:", "x": "Genes on the X chromosome removed:", "non_bi": "Non-biallelic SNPs removed:"}
    ss_t = None
    out_name = "filtered_" + ss_type + ".ss"
    #Some room to beautify the code here...
    if ss_type == "HELIX":
        #Specify dtypes to speed up read in. Or better yet, select just the columsn we want from the get go
        #SNP    CHR BP  A1  A2  P   BETA
        #dtype= { "SNP":'string_',"CHR": 'category', "BETA":'float64',"P":'float64', "A1":'category', "A2":'category'}
        #usecols = ["SNP", "CHR", "A1", "A2", "BETA", "P"]
        ss_t = pd.read_csv(ss, sep = '\t', dtype = {"CHR":"string_"})
        sizes["start"] = ss_t.shape[0]
        #Remove non-biallelic snps
        ss_t.drop_duplicates(subset = "SNP", keep = False, inplace = True)  
        sizes["non_bi"] = ss_t.shape[0]
        #Remove ambiguous SNPs
        ss_t["comb"] = ss_t.A1 + ss_t.A2 
        ss_t = ss_t[ss_t.comb.isin(noambig_list)]
        sizes["ambig"] = ss_t.shape[0]
        #Also remove SNPs with NA
        ss_t = ss_t[ss_t.BETA.notnull()]
        sizes["na"] = ss_t.shape[0]
        #Also remove X chromosome
        ss_t = ss_t[ss_t.CHR != "X"]
        sizes["x"] = ss_t.shape[0]
        ss_t["ID"] = ss_t.SNP + ":" + ss_t.A2 + ":" + ss_t.A1
        ss_t.rename(columns = {"BETA":BETA, "P":PVAL, "A2":REF, "A1":ALT}, inplace = True)
        ss_t = ss_t[["ID", REF, ALT, BETA, PVAL]]    
    elif ss_type == "NEAL":
        #variant    minor_allele    minor_AF    low_confidence_variant  n_complete_samples  AC  ytx beta    se  tstat   pval
        ss_t = pd.read_csv(ss, sep = '\t',usecols = ["variant", "beta", "pval", "minor_AF"], dtype = {"variant": 'string_', "beta":'float64', "pval":'float64', "minor_AF":"float64"})
        sizes["start"] = ss_t.shape[0]
        id_data = ss_t.variant.str.split(":", expand = True)  #gives chr:pos:ref:alt
        ss_t["comb"] = id_data[2] + id_data[3]#last 2 characters at the back end of the variant
        ss_t["CHR"] = id_data[0]
        ss_t["loc"] = id_data[0] + ":" + id_data[1]
        #Removing non_bi-allelic SNPs
        ss_t.drop_duplicates(subset = "loc", keep = False, inplace = True)
        sizes["non_bi"] = ss_t.shape[0]
        #Removing ambiguous SNPs
        ss_t = ss_t[ss_t.comb.isin(filter_list)]
        ss_t = ss_t[(ss_t.comb.isin(ambig_list) & (ss_t.minor_AF > 0.4)) | (ss_t.comb.isin(noambig_list))]
        #remove ones that are in the ambiguous list unless they have a MAF > 0.4
        sizes["ambig"] = ss_t.shape[0]
        #Also remove SNPs with NA
        ss_t = ss_t[ss_t.beta.notnull()] #here its NaN
        sizes["na"] = ss_t.shape[0]
        #Also remove X chromosome
        ss_t = ss_t[ss_t.CHR != "X"]
        sizes["x"] = ss_t.shape[0]
        #Set it up for write_out ID, ref, alt, beta, pval
        ss_t.rename(columns = {"variant":"ID", "beta":BETA, "pval":PVAL}, inplace = True)
        ss_t[REF] = id_data[2]
        ss_t[ALT] = id_data[3]
        ss_t = ss_t[["ID", REF, ALT, BETA,PVAL]]
    elif ss_type == "DEFAULT":
        print("Are you sure you want to re-filter the summary stats? It appears they may have already been filtered...")
        input()
        #[ID REF ALT beta pval]
        ss_t = pd.read_csv(ss, sep = '\t')
        sizes["start"] = ss_t.shape[0]
        ss_t["comb"] = ss_t[REF] + ss_t[ALT]
        #Removing non_bi-allelic SNPs
        id_data = ss_t.ID.str.split(":", expand = True)
        ss_t["new_id"] = id_data[0] + ":" + id_data[1] 
        ss_t.drop_duplicates(subset = "new_id", keep = False, inplace = True)
        sizes["non_bi"] = ss_t.shape[0]
        #Removing ambiguous SNPs
        ss_t = ss_t[ss_t.comb.isin(noambig_list)] #Drawback of this format- doesn't allow for AF
        sizes["ambig"] = ss_t.shape[0]
        #Also remove SNPs with NA
        ss_t = ss_t[ss_t.beta.notnull()] #here its NaN
        sizes["na"] = ss_t.shape[0]
        #Also remove X chromosome
        ss_t['CHR'] = id_data[0]
        ss_t = ss_t[ss_t.CHR != "X"]
        sizes["x"] = ss_t.shape[0]
        #Set it up for write_out ID, ref, alt, beta, pval
        ss_t = ss_t[["ID", REF, ALT, BETA,PVAL]] 
    else:
        call_str = "echo Nothing given"
    end_size = ss_t.shape[0]
    #Write size changes out to log file
    prev = "start"
    for t in sizes:
        if t == "start":
            updateLog(str(sizes[t]) + " variants submitted in summary statistics file.")
        else:
            updateLog(comment[t], str(sizes[prev]-sizes[t]))
            prev = t
    #Write out to file #ID, ref, alt, beta, pvalue
    #Check to see if any changes even made
    if ss_t.shape[0] == sizes["start"]:
        updateLog("Summary stats alreadya appear to be filtered. No changes were made")
    ss_t.to_csv(out_name, sep = '\t', index = False) 
    return out_name

def plinkToMatrix(snp_keep, args, local_pvar, vplink):
    """
    Create a matrix with scores we can parse
    @param snp_keep is the path to the file with the SNP list. Currently just a single column, TODO check if this will work for plink filtering.
    """
    #Write out the list of SNPs to keep
    #if os.path.isfile("mat_form_tmp.bed"):
    if os.path.isfile("mat_form_tmp.raw"):
        print("Using currently existing genotype matrix...")
    else:
        syscall = " --pfile "
        annot =  " --pvar "
        if vplink == 1:
            syscall = " --bfile "
            annot = " --bim "
        call = "plink2 " + syscall  + args + annot + local_pvar + " --geno --not-chr X --extract " + snp_keep + " --out mat_form_tmp --export A --max-alleles 2"
        check_call(call, shell = True) 
        #--geno removes snps with more than 10% missing from the data.
        print("User-readable genotype matrix generated")
        if not os.path.isfile("mat_form_tmp.raw"): 
            print("PLINK file was not correctly created. Program will terminate.")
            sys.exit()
    return "mat_form_tmp"

#TODO- you can do this once with the 1000 genomes cohort and then just repeat with that list of "best" genes.
#Doesn't need to be done denovo every time!
#TODO: fix the id system for clumping, make sure consisent
def plinkClump(reference_ld, clump_ref,clump_r, geno_ids, ss):
    """
    Run plink and clump, get the new list of ids.
    Select just these from the summary statistics.
    """
    pre_clump_count = getEntryCount(geno_ids, False)
    if clump_ref == "NA":
        #Build the file for plink to run on 
        command = "echo SNP P > clump_ids.tmp"
        check_call(command, shell = True) 
        command = "awk ' {print $1,$NF} ' " + ss + ">> clump_ids.tmp" #Gives SNP id and p-value to assist in clumping
        check_call(command, shell = True)
        #Run plink clumping
        plink_command = "plink --bfile " + reference_ld + " --clump clump_ids.tmp --clump-best --clump-p1 1 --clump-2 1 --clump-r2 " + str(clump_r) + " --clump-kb 250"
        check_call(plink_command, shell = True)
        clump_file = "plink.clumped.best" #The default output name with flag clump-best
    else:
        print("Using provided reference clump...")
        clump_file = clump_ref
    #Get a way to map from one to the other
    temp_ids = '''cut -f 1,2 -d ":" ''' + ss + " > t && paste  t " + ss + " > id_mapper.tmp && rm t"
    check_call(temp_ids, shell = True)
    command = "awk '(FNR == NR) {a[$1];next} ($1 in a) {print $0}' " + clump_ref + " id_mapper.tmp | cut -f 2,3,4,5,6  > t && mv t "+ ss
    check_call(command, shell = True)
    command = "rm id_mapper.tmp"
    check_call(command, shell = True)
    #command = "awk '(FNR == NR) {a[$1];next} ($1 in a) {print $1}' " + clump_ref + " " + geno_ids + " > t && mv t " + geno_ids
    command = "cut -f 1 -d ' ' " + ss +" > " + geno_ids
    check_call(command, shell = True)
    #Log reporting
    post_clump_count = getEntryCount(geno_ids, False)
    updateLog("Number of SNPs removed by clumping (selecting only top clump variant):", str(pre_clump_count - post_clump_count)) 
    return ss, geno_ids
    

def getPatientIDs(snp_matrix):
    tnames = list(snp_matrix.columns)[6:]
    return [int(i.split("_")[0]) for i in tnames]


def writeScores(scores, ids, destination, debug = False):
   #May be able to just do it faster in PANDAS:
    scores["IID"] = ids
    tab_out = pd.DataFrame.from_dict(scores)
    tab_out = tab_out.set_index("IID")
    tab_out.to_csv(destination, sep = '\t') 
    return

def writeScoresDebug(debug_tab, destination):
    with open(destination, 'w') as ostream:
        p_len = dict()
        p_string = ""
        for p in debug_tab:
            p_string  = p_string + str(p) + '\t'
            p_len[p] = len(debug_tab[p])
        ostream.write(p_string[:-1] + '\n')
        max_val = max(p_len.values())
        for i in range(0,max_val):
            write_s = "" 
            for p in debug_tab:
                if i < p_len[p]:
                    write_s = write_s + str(debug_tab[p][i]) + '\t'
                else:
                    write_s  = write_s + '\t'
            ostream.write(write_s[:-1] + '\n')

def writeSNPIDs(ids):
    ids.to_csv("snp_ids.tmp", header = False, index = False, sep = '\t') #(should have 2 columns of the family/infamily ids, which appear to be the same.)
    return "snp_ids.tmp"

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Basic tools for calculating rudimentary Polygenic Risk Scores (PRSs). This uses PLINK bed/bim/fam files and GWAS summary stats as standard inputs. Please use PLINK 2, and ensure that there are no multi-allelic SNPs")
    parser.add_argument("-snps", "--plink_snps", help = "Plink file handle for actual SNP data (omit the .extension portion)")
    parser.add_argument("--preprocessing_done", help = "Specify this if you have already run the first steps and generated the plink2 matrix with matching IDs. Mostly for development purposes.")
    parser.add_argument("-ss", "--sum_stats", help = "Path to summary stats file. Be sure to specify format if its not DEFAULT. Assumes tab delimited", required = True)
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE", "NEAL", "HELIX"])
    parser.add_argument("-p", "--pval", default = 5e-8, type = float, help = "Specify what pvalue threshold to use.")
    parser.add_argument("--pvals", help= "Use this if you wish to specify multiple at once, in a list separated by commas")
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--prefiltered_ss", default = False, action = "store_true", help = "Specify this option if you have already filter summary stats for ambiguous SNPs, bi-alleleic variants, NAs. Recommended for speed.")
    parser.add_argument("--reset_ids", action = "store_true", help = "By default, program assumes genotype ids are in the form chr:loc:ref:alt. Specify this argument if this is NOT the case to.")
    parser.add_argument("--clump", help = "Specify if you wish to do clumping. Option to pass in a plink clumping.best file or the reference panel plink data.", action = "store_true", required = "--ld_ref" in sys.argv or "--clump_ref" in sys.argv)
    parser.add_argument("--ld_ref", default = None, help = "Specify an LD reference panel, required if you want to do clumping. Pass in plink1.9 format please.")
    parser.add_argument("--clump_r2", default = 0.25, help = "Specify the clumping threshold.")
    parser.add_argument("--debug", help = "Specify this if you wish to output debug results", default = False, action = "store_true")
    parser.add_argument("--clump_ref", help = "Specify this option if you wish to perform clumping and already have a plink clump file produced (the output of plink1.9 with the --clump-best argument). Provide the path to this file.", default = "NA")
    parser.add_argument("-pv", "--plink_version", default = 2, help = "Specify which version of plink files you are using for the reference data.", type = int, choices = [1,2])
    parser.add_argument("--align_ref", action = "store_true", help = "Specify this if you would like us to align the reference alleles between the target and training data. If this is already done, we recommend omitting this flag since it will create an additional copy of the genotype data.")
    parser.add_argument("--only_preprocess", action = "store_true", help = "Specify this option if you would like to terminate once preprocessing of SNPs (including aligning reference, removing ambiguous SNPs, preparing summary stat sublist) has completed.")
    parser.add_argument("--OVERRIDE", help = "Use for deubgging- pass the matrix in")
    args = parser.parse_args()
    pvals = [args.pval]
    DEBUG = args.debug
    start = time.time()    
    #Extract the pvalues to asses at
    if args.pvals:
        t = (args.pvals).strip().split(",")
        pvals = [float(x) for x in t]
        updateLog("Pvals to asses:", str(args.pvals))
    if args.reset_ids:
        VAR_IN_ID = False #Variant information is not included in genotype id
    if args.preprocessing_done:
        print("Preprocessing has already been completed. Proceeding with PRS calculations...")
        print("Not implemented error, will run as normal")
        #Unload everything in the pickle, which contains the path information for everything:
            #The snp_indices data, the snp_matrix path
    if not args.prefiltered_ss:
        print("Ambiguous SNPs and indels have not been filtered out from SS data. Filtering now...")
        args.sum_stats = filterSumStatSNPs(args.sum_stats, args.ss_format)
        args.ss_format = "DEFAULT" #ID REF ALT BETA PVAL
        print("Ambiguous SNPs and indels removed. New file is", args.sum_stats)
    if args.align_ref:
        args.plink_snps = alignReferenceByPlink(args.plink_snps, args.plink_version, args.sum_stats, args.ss_format) #Make a local copy of the data with the updated genotype
        print("A new plink2 fileset with matched reference sequences has been generated", args.plink_snps)
        print("Strand flips have been corrected for. Proceeding with analysis")
        args.plink_version =2 #We have updated the type to 2
    
    print("Selecting IDs for analysis...")
    local_pvar, geno_ids, ss_parse = prepSNPIDs(args.plink_snps, args.sum_stats,args.ss_format, vplink = args.plink_version, max_p = max(pvals))
    num_pat = getPatientCount(args.plink_snps, args.plink_version)
    updateLog("Number of patients detected in sample:", str(num_pat))
    print("Time for preprocessing (SNP filtering, reference alignment if specified, etc.", str(time.time() - start))
    print("Reading data into memory (only done on first pass)")
    if args.clump:
        ss_parse, geno_ids = plinkClump(args.ld_ref, args.clump_ref,args.clump_r2, geno_ids, ss_parse)
    print("Preprocessing complete")
    if args.only_preprocess:
        sys.exit()
    
    updateLog("Total number of SNPs included in analysis", str(getEntryCount(geno_ids,False)))  
    stats_complete = readInSummaryStats(ss_parse)
    snp_indices = indexSSByPval(pvals,stats_complete)
    if args.OVERRIDE:
        snp_matrix = args.OVERRIDE
    else:
        snp_matrix = plinkToMatrix(geno_ids, args.plink_snps, local_pvar, args.plink_version) 
    print("Parsing the genotype file...")
    scores, patient_ids, debug_dat = calculatePRS(pvals, snp_matrix, snp_indices,stats_complete[REFID])
    
    writeScores(scores, patient_ids, args.output + ".tsv")
    print("Scores written out to ", args.output + ".tsv")
    stop = time.time()
    print("Total runtime:", str(stop - start))
    #TODO: Add cleanup functionality, remove the massive files you've made.


