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
from subprocess import check_call 
from subprocess import call
import os
import time
import multiprocessing as mp
from operator import itemgetter
#from plinkio import plinkfile
#General reference things
CHR = "chr"
POS = "pos"
REFID = "ref_id"
SNP = "snp"
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
    ss = pd.read_csv(s_path, sep = ' ', header = None,  names = [REFID, REF,ALT, BETA, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64', SEBETA:'float64', TSTAT:'float64'})
    #ss = pd.read_csv(s_path, sep = ' ', header = None,  names = [REFID, REF,ALT, BETA,SEBETA,TSTAT, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64', SEBETA:'float64', TSTAT:'float64'})
    print("Data in memory...")
    #ss.info(memory_usage='deep') 
    
    return ss

def filterSumStats(pvals,ss):
    pvals_ref = dict() #dictionary containing a pointer for each p value
    pvals.sort() #smallest one first
    
    for p in pvals:
        samp = ss[(ss[PVAL] < float(p))]
        pvals_ref[p] = [samp.index.tolist(), samp[BETA].values] #Get just the pvalues and the indices.
    """
    #Alternative version:
    prev = -1
    for p in pvals:
        samp = ss[(ss[PVAL] < float(p)) & (ss[PVAL] >= prev)]
        pvals_ref[p] = [samp.index.tolist(), samp[BETA].values]
        prev = float(p)
    """
    return pvals_ref
        
def prepSummaryStats(geno_ids, sum_path, pval_thresh):
    """
    This function uses awk to filter the summary stats by pvalue, and then reads them into a dask data structure.
    """
    #If the id is in the geno_ids file, and it passes the p-value threshold, keep it
    command = "awk '($7 <= " + str(pval_thresh) + ") {print $0}' " + sum_path + " > ss.tmp"
    check_call(command, shell = True) 
    #Pro tip- in benchmarks, cut performs better here than awk
    command = "cut -f 1 -d ' ' ss.tmp > ss_ids.tmp"
    check_call(command, shell = True)
    #TODO: determine if the above awk/cut commands are faster than dask filtering + computing the IDs. Or maybe just doing that with pandas.
    return "ss.tmp", "ss_ids.tmp"

#This will extract the ids we want to be working with.
def prepSNPIDs(snp_file, ss_file, ss_type, ambig, id_type, vplink = 2, max_p = 1):
    """
    This extracts all the IDs from the genotype data (pvar file) that we wish to be using and selects just the data from the summary stats data we want
    @return path to the genotype ids
    @return path to the summary stats.
    """
    local_geno = snp_file + ".pvar"
    max_p = str(max_p)
    if id_type:
        local_geno = "local_geno.pvar" #We will be making our own pvar
        if vplink == 1:
            local_geno = "local_geno.bim"

    #Remove all the header lines, just get what we need.
    if not os.path.isfile('geno_ids.f'):
        #Step 1: format the PVAR file to have the IDs as we wouldl like them
        #TODO: update id_type flag to bypass this step if not necessary.
        #By default, assume the IDs are in the right way.
        command = '''awk '(!/##/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' ''' + snp_file + ".pvar  > local_geno.pvar"
        #Simplifying things, we aren't going to re-write IDs....
        if id_type: #This argument means the IDs in the genotype data are NOT as we wish them to be. Need to modify 
            command = ''' awk '(!/##/ && /#/) {print $1"\t"$2"\tID\t"$4"\t"$5} (!/#/) {print $1"\t"$2"\t"$1":"$2":"$4":"$5"\t"$4"\t"$5}' ''' + snp_file + ".pvar > local_geno.pvar"
            if vplink == 1:
                #Note that this puts the allele in chr:num:minor:major
                command = '''awk '{print $1"\t"$1":"$4":"$5":"$6"\t"$3"\t"$4"\t"$5"\t"$6}' ''' + snp_file + ".bim > local_geno.bim"
                updateLog("Number of variants detected in original genotype file", str(getEntryCount(snp_file + ".bim", False)))
            else:
                updateLog("Number of variants detected in original genotype file", str(getEntryCount(snp_file + ".pvar", True)))
        try:
            check_call(command, shell = True)
        except subprocess.CalledProcessError:
            print("Unable to generate modified pvar/bim file. Are you sure you specified the correct plink version?")
            print("Failed on command", command)
            sys.exit()
        #Step 2: Generate file of the IDs from genotype data for filtering of the summary statistics
        if vplink == 2:
            command_n = "awk ' (NR> 1 && $1 !~ /X/) {print $3}'  local_geno.pvar > geno_ids.f"
        else:
            command_n = "awk ' ($1 !~ /X/) {print $2}'  local_geno.bim > geno_ids.f"
        check_call(command_n, shell = True)
    if not os.path.isfile("ss_filt.f"):
        #Note that for all types, we expect IDs to be of the form chr:#:ref:alt
        #CHANGE on 10-22: only outputing the pare minimum, ID, ref, alt, beta, pvalue. The rest isn't needed
        #In cases where the REF and ALT alleles are not specified, it is ID, major, minor/effect, beta, pvalue
        #11/14- adding functionality to filter out SNPs from our max p-value and below at this point.
        if ss_type == "SAIGE":
            #ID, REF, ALT, BETA, SEBETA, Tstat, pval
             command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($3 in a && $13 <= ''' + max_p + ''' ) {print $3, $4,$5,$10, $13}' geno_ids.f ''' + ss_file + " > ss_filt.f"  
            #command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($3 in a) {print $3, $4,$5,$10,$11,$12,$13}' geno_ids.f ''' + ss_file + " > ss_filt.f"  
        elif ss_type == "NEAL":
            command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a && $11 <= ''' + max_p + ''' ) {print $1, "N",$2,$8,$11}' geno_ids.f ''' + ss_file + " > ss_filt.f"    
            #command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a) {print $1, "N",$2,$8,$9,$10,$11}' geno_ids.f ''' + ss_file + " > ss_filt.f"   
        elif ss_type == "TAIBO": #SNP chr BP, A1 A2 P Beta --> ID ref alt beta pvalue
            command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1":"$5":"$4 in a && $6 <= ''' + max_p + ''' ) {print $1":"$5":"$4, $5,$4,$7,$6 }' geno_ids.f ''' + ss_file + " > ss_filt.f"
        elif ss_type == "DEFUALT": #The kind of output we like, with a header ID ALT REF BETA SNP
            command = '''awk '(FNR == NR) {a[$1];next} (FNR == 1 && NR == 2) {next;} ($1 in a && $4 <= ''' + max_p + ") {print $0}' geno_ids.f " + ss_file + " > ss_filt.f"
        else:
            print("The type of ss file hasn't been specified. Please specify this.")
    
        #command = '''awk '(FNR == 1) {next;} {print $1":"$2, $4,$5,$10,$11,$12,$13}' ''' + ss_file + " > ss_filt.f"
        check_call(command, shell = True)
        #reset the ids we use downstream
        command = "cut -f 1 -d ' ' ss_filt.f > geno_ids.f"
        check_call(command, shell = True)
    return local_geno,"geno_ids.f", "ss_filt.f"

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

#There must be a better way than this...
def hashMapBuild(ss_snps, plink_snps, ret_tap): #we need to deal with this issue of multiallelic snps
    pref = dict()
    
    #Put all of the plink genotype SNPs into a dictionary for fast lookup
    for i in range(0, len(plink_snps)):
        if VAR_IN_ID:
            #id_name= str(plink_snps[i].name)
            pref[plink_snps[i][:-2]] = i
            #print(id_name)
        else:
            print("This portion of code is in development. Please finish, ashton")
            id_name = str(plink_snps[i]) + ":" + str(plink_snps[i]) + ":" + str(plink_snps[i]) #why is this order so? nervous making
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
        if plink_order[rand][:-2] != ss_order.iloc[rand]:
        #if plink_order[rand].name != ss_order.iloc[rand]:
            #print("It didn't work as we wanted.....", plink_order[rand].name, ss_order.iloc[rand])
            print("SNP order not aligned, realigning now...")
            print(plink_order[rand])
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
    return ret
    
def linearGenoParse(snp_matrix, snp_indices,pvals, scores, ss_ids):
    """Parse the patient file

    Currently in development, Oct 31
    snp_matrix -- the path to the massive genotype file
    snp_indieces -- the map with the indices for each pvalue we care about
    scores -- what is returned, stores the scores
    pvals -- which pvalues we measure at
    ss_ids -- the list of ids from the summary stats list, for ordering purposes
    """
    from operator import itemgetter
    patient_ids = list()
    report_index = 100
    pvals.sort()
    with open(snp_matrix + ".raw", 'r') as istream:
        line_counter = 0
        dat = istream.readline().split('\t')
        line_len = len(dat)
        while dat:
            if line_counter == 0:
                var_map = buildVarMap(dat[6:], ss_ids)
                if DEBUG:
                   DEBUG_DAT[REFID] = dat[6:] 
                print("Proceeding with line-by-line calculation now...")
                for p in pvals:
                    updateLog("Number of final SNPs used at ", str(p), "pvalue level:", str(len(snp_indices[p][INDEX]))) 
                #Need to check that this goes the way we expect.

            else:
                patient_ids.append(dat[0])
                
                for p in pvals:
                    snp_index = snp_indices[p][INDEX]
                    sel = var_map[snp_index]
                    try:
                        res_list = np.array(itemgetter(*sel)(dat[6:]), dtype = np.float64)
                        scores[p].append(scoreCalculation(res_list, snp_indices[p][B])) 
                    
                    except:
                        print("Unable to read line", line_counter)
                        print(dat)
                        scores[p].append("NA")
                # faster way to implement?:
                """
                first_pval = True
                rel_data = dat[6:]
                prev_score = 0
                new_score = 0
                for p in pvals:
                    if first_pval:
                        snp_index = snp_indices[p][INDEX]
                        sel = var_map[snp_index]
                        new_genos = np.array(itemgetter(*sel)(rel_data), dtype = np.float64)
                        new_score = scoreCalculation(new_genos, snp_indices[p][B])
                        scores[p].append(new_score)
                        first_pval = False
                    else:
                        snp_index = snp_indices[p][INDEX]
                        sel = var_map[snp_index]
                        new_genos = np.array(itemgetter(*sel)(rel_data), dtype = np.float64)
                        new_score = scoreCalculation(new_genos, snp_indices[p][B]) + prev_score
                        scores[p].append(new_score)    
                    prev_score = new_score     
                    #faster because only extracting each snp once instead of p times for some of them
                    #only doing m dots instead of m + m-p + m-p2 + ....
            """ 
            try:   
                dat = istream.readline().split('\t')
            except: #end of the file
                dat = None
            line_counter += 1

            if line_counter % report_index == 0:
                print("Currently at individual/sample", line_counter)
                updateLog("Currently at individual/sample", str(line_counter))
            if len(dat) <= 1:
                print("Finished all the samples!")
                break
    if DEBUG:
        print("Writing out requested debug information")
        DEBUG_DAT[BETA] = dict()
        DEBUG_DAT[SNP] = dict()
        for p in pvals:
            snp_index = snp_indices[p][INDEX]
            sel = var_map[snp_index] 
            #debug_betas[p] = snp_indices[p][B]
            #debug_snps[p] = sel #just the ids, not the names....
            DEBUG_DAT[BETA][p] = snp_indices[p][B]
            DEBUG_DAT[SNP][p] = ((itemgetter(*sel)(DEBUG_DAT[REFID])))
    return scores, patient_ids, DEBUG_DAT #basically a pval:[list of snps]  
    if line_counter == 1:
        print("There appears to be some error with the genotype plink matrix. Please ensure the genotype data is correct or the plink call was not interrupted.")
        sys.exit()
    return scores, patient_ids, None


def calculatePRS(pvals, snp_matrix, snp_indices, snp_list):
    #Seraj shoutout
    start = time.time()
    scores = { p : [] for p in pvals}
    scores, patient_ids, debug_tab = linearGenoParse(snp_matrix, snp_indices,pvals, scores, snp_list)
    end = time.time()
    print("Calculation and read in time:", str(end-start))
    return scores, patient_ids, debug_tab

def sameAlleleAssesment(s1, s2):
    t = (s1 != s2).sum()
    if t == 0:
        return True
    else:
        print(s1, s2)
        print("^Unusual case")
        return False

def alignReferenceByPlink(old_plink,plink_version, ss, ss_type):
    #Make the reference file with major 
    #get the reference information you want
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
    elif ss_type == "TAIBO":
        #ID CHR BP A1(effect) A2 P BETA
        quick_call = "awk ' ( NR > 1) {print $1,$5}' " + ss + " | uniq > ./aligner.t"
        check_call(quick_call, shell = True) 
    elif ss_type == "DEFAULT": #ID REf ALT BETA Pvalue
        call = "cut -f 1,2 | uniq > aligner.t"
        check_call(call)
    else:
        print("Have not extended the method for type", ss_type, ". Please align Reference information manually")
        sys.exit()
    ftype = " --pfile "
    if plink_version == 1:
        ftype = " --bfile "
    #plink_call = "plink2" + ftype + old_plink  + " --not-chr X --out new_plink --ref-allele force aligner.t --make-pgen" 
    plink_call = "plink2" + ftype + old_plink  + " --not-chr X --out new_plink --ref-allele aligner.t --make-pgen --max-alleles 2" 
    check_call(plink_call, shell = True)
    clean_up = "rm *.t"
    check_call(clean_up, shell = True)
    return "new_plink"

def filterAmbiguousSNPs(ss, ss_type, thresh = str(0.4)):
    #okay, I'm biting the bullet and reading it in...
    #there are 4 stages:
    #0) Remove anything more than bi-allelic
    #1) Limit it to SNPs only
    #2) Limit to nonambiguous SNPs, if possible
    #3) Remove null beta values
    #4) Remove X chromosome
    #5) Write it out in the default format so its easier to digest downstream. (
    filter_list = {"AG", "AC", "AT", "TA", "TC", "TG", "GC", "GT", "GA", "CG", "CT", "CA"}
    noambig_list = {"AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG"}
    ambig_list = {"AT", "TA", "GC", "CG"}
    sizes = {"start": 0,"non_bi": 0, "ambig":0, "na":0, "x":0}
    comment = {"ambig": "Ambiguous SNPs along with INDELS removed:", "na": "Variants with missing effect sizes removed:", "x": "Genes on the X chromosome removed:", "non_bi": "Non-biallelic SNPs removed:"}
    ss_t = None
    out_name = "noAmbig_lowfilt.ss"
    if ss_type == "TAIBO":
        #Specify dtypes to speed up read in. Or better yet, select just the columsn we want from the get go
        #SNP    CHR BP  A1  A2  P   BETA
        #dtype= { "SNP":'string_',"CHR": 'category', "BETA":'float64',"P":'float64', "A1":'category', "A2":'category'}
        #usecols = ["SNP", "CHR", "A1", "A2", "BETA", "P"]
        ss_t = pd.read_csv(ss, sep = '\t')
        sizes["start"] = ss_t.shape[0]
        #Remove non-biallelic snps
        s_t.drop_duplicates(subset = "SNP", keep = False, inplace = True)  
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
        ss_t = ss_t["ID", REF, ALT, BETA, PVAL]    
    elif ss_type == "NEAL":
        #usecols= ["variant", "beta", "pval"]
        #dtype = {"variant": 'string_', "beta":'float64', "pval":'float64'} 
        #variant    minor_allele    minor_AF    low_confidence_variant  n_complete_samples  AC  ytx beta    se  tstat   pval
        ss_t = pd.read_csv(ss, sep = '\t')#there are ways to speed this up, consider them
        sizes["start"] = ss_t.shape[0]
        id_data = ss_t.variant.str.split(":", expand = True) 
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
        ss_t = ss_t["ID", REF, ALT, BETA,PVAL]
    else:
        call_str = "echo Nothing given"
    #call_str = call_str.replace("0.4", thresh)
    #check_call(call_str, shell = True)
    end_size = ss_t.shape[0]
    prev = "start"
    for t in sizes:

        if t == "start":
            updateLog(str(sizes[t]) + " variants submitted in summary statistics file.")
        else:
            updateLog(comment[t], str(sizes[prev]-sizes[t]))
            prev = t
    #Write out to file #ID, ref, alt, beta, pvalue
    ss_t = ss_t.drop("comb", axis = 1)
    ss_t.to_csv(out_name, sep = '\t', index = False) 
    return out_name

def plinkToMatrix(snp_keep, args, local_pvar, vplink):
    """
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
            #check_call(["plink2", "--bfile", args, "--bim", local_pvar, "--not-chr", "X", "--extract", snp_keep, "--out", "mat_form_tmp", "--export", "A", "--max-alleles", "2"], shell = True)
        call = "plink2 " + syscall  + args + annot + local_pvar + " --geno --not-chr X --extract " + snp_keep + " --out mat_form_tmp --export A --max-alleles 2"
        check_call(call, shell = True) 
        #--geno removes snps with more than 10% missing from the data.
        
        print("New plink matrix made!")
        if not os.path.isfile("mat_form_tmp.raw"): #C
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
        check_call(command, shell = True) 
        command = "awk ' {print $1,$NF} ' " + ss + ">> clump_ids.tmp"
        check_call(command, shell = True)
        #Run plink clumping
        plink_command = "plink --bfile " + reference_ld + " --clump clump_ids.tmp --clump-best"
        check_call(plink_command, shell = True)
        clump_file = "plink.clumped.best"
    else:
        print("Using provided reference clump...")
        clump_file = clump_ref
    #Refilter summary stats....
    #Actually, I should be able to just filter the list based on the plink input... buttt then we still have an order issue. Which I guess doesn't matter but.
    command = "awk '(FNR == NR) {a[$1];next} ($1 in a) {print $0}' " + clump_ref + " " + ss + " > t && mv t " + ss
    check_call(command, shell = True)
    #command = "awk '(FNR == NR) {a[$1];next} ($1 in a) {print $1}' " + clump_ref + " " + geno_ids + " > t && mv t " + geno_ids
    command = "cut -f 1 -d ' ' " + ss +" > " + geno_ids
    check_call(command, shell = True)
    return ss, geno_ids
    

def getPatientIDs(snp_matrix):
    tnames = list(snp_matrix.columns)[6:]
    return [int(i.split("_")[0]) for i in tnames]


def writeScores(scores, ids, destination, debug = False):
    with open(destination, 'w') as ostream:
        p_string = "IID\t"
        num = 0
        for p in scores:
            p_string = p_string + str(p) + '\t'
            num = len(scores[p]) #number of patients
        ostream.write(p_string[:-1] + '\n')
        
        for i in range(0, num):
            #write_s = str(ids[i].iid) + '\t'
            write_s = str(ids[i]) + '\t'
            for p in scores:
                write_s = write_s + str(scores[p][i]) + '\t'
            write_s = write_s[:-1] + '\n'
            ostream.write(write_s)

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
            for p in scores:
                if i < p_len[p]:
                    write_s = write_s + str(debug_tab[p][i]) + '\t'
                else:
                    write_s  = write_s + '\t'
            ostream.write(write_s[:-1] + '\n')

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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Basic tools for calculating rudimentary Polygenic Risk Scores (PRSs). This uses PLINK bed/bim/fam files and GWAS summary stats as standard inputs. Please use PLINK 2, and ensure that there are no multi-allelic SNPs")
    parser.add_argument("-snps", "--plink_snps", help = "Plink file handle for actual SNP data (omit the .extension portion)")
    parser.add_argument("--preprocessing_done", help = "Specify this if you have already run the first steps and generated the plink2 matrix with matching IDs. Mostly for development purposes.")
    parser.add_argument("-ss", "--sum_stats", help = "Path to summary stats file. Be sure to specify format if its not DEFAULT. Assumes tab delimited", required = True)
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE", "NEAL", "TAIBO"])
    parser.add_argument("-p", "--pval", default = 5e-8, type = float, help = "Specify what pvalue threshold to use.")
    parser.add_argument("--pvals", help= "Use this if you wish to specify multiple at once, in a list separated by commas")
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--no_ambig", default = False, action = "store_true", help = "Specify this option if you have already filter for ambiguous SNPs and bi-alleleic variants. Recommended for speed.")
    parser.add_argument("--var_in_id", action = "store_true", help = "Specify this if you wish to use more specific id of the format chr:loc:ref:alt. For default SNP in file, just leave this.")
    parser.add_argument("--clump", help = "Specify if you wish to do clumping. DEfault threshold is 0.5", action = "store_true", required = "--ld_ref" in sys.argv or "--clump_ref" in sys.argv)
    parser.add_argument("--ld_ref", help = "Specify an LD reference panel, required if you want to do clumping. Pass in plink1.9 format please.")
    parser.add_argument("--debug", help = "Specify this if you wish to output debug results", default = False, action = "store_true")
    parser.add_argument("--clump_ref", help = "Specify this option if you wish to perform clumping and already have a plink clump file produced (the output of plink1.9 with the --clump-best argument). Provide the path to this file.", default = "NA")
    parser.add_argument("-pv", "--plink_version", default = 2, help = "Specify which version of plink files you are using for the reference data.", type = int, choices = [1,2])
    parser.add_argument("--align_ref", action = "store_true", help = "Specify this if you would like us to align the reference alleles between the target and training data. If this is already done, we recommend omitting this flag since it will create an additional copy of the genotype data.")
    parser.add_argument("--only_preprocess", action = "store_true", help = "Specify this option if you would like to terminate once preprocessing of SNPs (including aligning reference, removing ambiguous SNPs, preparing summary stat sublist) has completed.")
    args = parser.parse_args()
    pvals = [args.pval]
    DEBUG = args.debug
    start = time.time()    
    if args.pvals:
        pvals = [float(x) for x in (args.pvals).split(",")]
        updateLog("Pvals to asses:", str(args.pvals))
    if args.var_in_id:
        VAR_IN_ID = True
    #there are ways to improve this....
    if args.preprocessing_done:
        print("Preprocessing has already been completed. Proceeding with PRS calculations...")
        #Unload everything in the pickle, which contains the path information for everything:
            #The snp_indices data, the snp_matrix path
    if not args.no_ambig:
        print("Ambiguous SNPs and indels have not been filtered out from SS data. Filtering now...")
        args.sum_stats = filterAmbiguousSNPs(args.sum_stats, args.ss_format)
        args.ss_format = "DEFAULT" #ID REF ALT BETA PVAL
        print("Ambiguous SNPs and indels removed. New file is", args.sum_stats)
    if args.align_ref:
        args.plink_snps = alignReferenceByPlink(args.plink_snps, args.plink_version, args.sum_stats, args.ss_format) #Make a local copy of the data with the updated genotype
        print("A new plink2 fileset with matched reference sequences has been generated", args.plink_snps)
        print("Strand flips have been corrected for. Proceeding with analysis")
        args.plink_version =2 #We have updated the type to 2
    
    print("Selecting IDs for analysis...")
    max_p = max(pvals)
    local_pvar, geno_ids, ss_parse = prepSNPIDs(args.plink_snps, args.sum_stats,args.ss_format, args.no_ambig, args.var_in_id, args.plink_version, max_p)
    num_pat = getPatientCount(args.plink_snps, args.plink_version)
    updateLog("Number of patients detected in sample:", str(num_pat))
    updateLog("Total number of SNPs included in analysis", str(getEntryCount(geno_ids,False)))  
    print("Time for preprocessing (SNP filtering, reference alignment if specified, etc.", str(time.time() - start))
    print("Reading data into memory (only done on first pass)")
    #Selecting out SNPs that are likely in LD
    #at this point, the geno ids and ss parse list  has what we want. We will reduce it more though through clumping via plink.
    if args.clump:
        ss_parse, geno_ids = plinkClump(args.ld_ref, args.clump_ref, geno_ids, ss_parse)
    if args.only_preprocess:
        print("Preprocessing complete")
        sys.exit()
   
    stats_complete = readInSummaryStats(ss_parse)
    snp_indices = filterSumStats(pvals,stats_complete)
    #time for benchmarking
    snp_matrix = plinkToMatrix(geno_ids, args.plink_snps, local_pvar, args.plink_version) 
    print("Parsing the genotype file...")
    scores, patient_ids, debug_dat = calculatePRS(pvals, snp_matrix, snp_indices,stats_complete[REFID])
    writeScores(scores, patient_ids, args.output + ".tsv")
    if DEBUG:
        writeScoresDebug(debug_dat[SNP], "debug_vals_snps.tsv")
        writeScoresDebug(debug_dat[BETA], "debug_vals_betas.tsv")
    #Match SNPs by address and/or reference ID.
    print("Scores written out to ", args.output + ".tsv")
    stop = time.time()
    print("Total runtime:", str(stop - start))
    #TODO: Add cleanup functionality, remove the massive files you've made.


