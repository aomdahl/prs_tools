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
import datetime
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
START_SEQ="Start of run:"
INDEX = 0
B = 1
global VAR_IN_ID
MEMORY = ""
VAR_IN_ID = True #Variants are included in ID, assumed by default
DEBUG = True        
DEBUG_DAT = dict()
LOG="log_file.txt"
################
class MissingGenoHandler:
    def __init__(self, ss_threshold_lengths, nsamps, pvals): #Default constructor
        self.running_totals = dict()
        self.num_dropped = dict()
        self.num_samps = nsamps
        self.missingness_index = dict() #Object will have structure of pvalue:[[[5,6,7],[beta1,beta2,beta3]],[[1,2,3],[beta1,beta2,beta3]]]
        self.pvals = pvals.sort() #smallest to largest.
        self.sample_means = dict()
        self.score_modifier = dict()
        self.pvals = pvals
        for p in pvals:
            self.missingness_index[p] = list()
            self.running_totals[p] = np.zeros(ss_threshold_lengths[p])
            self.sample_means[p] = None
            self.score_modifier[p] = None
            self.num_dropped[p] = np.zeros(ss_threshold_lengths[p])

    #getters and setters
    def setSNPCount(self, nsnps):
        self.running_totals = np.zeros(nsnps)

    def getNACount(self):
        drops = self.num_dropped.values()
        return np.sum([np.count_nonzero(x) for x in drops])

    def __calculateMeans(self):
        if self.num_samps == 0:
            print("Number of snps hasn't been set, please set")
            sys.exit()
        for p in self.pvals:
            self.sample_means[p] = self.running_totals[p]/(self.num_samps - self.num_dropped[p])

    def updateTotal(self, genos, p):
        self.running_totals[p]  += genos

    def updateMissingList(self, na_indices, weights, p):
        self.missingness_index[p].append([na_indices, weights])
        
    def handleNAs(self, sample_genos, weights, pval):
        NA_I=0 #Index for NA list
        na_locations = np.where(np.isnan(sample_genos)) #based on the shared weight/snp index.
        if len(na_locations[NA_I]) < 1: #NO Nas
            self.updateTotal(sample_genos, pval)
            self.updateMissingList([],[],pval)
        else: #Omit that patient.
            self.updateTotal(np.nan_to_num(sample_genos), pval)
            self.updateMissingList(na_locations[NA_I], itemgetter(*na_locations)(weights),pval) #should choose out correct weights.
            self.num_dropped[pval][na_locations[NA_I]] += 1 #increment the dropped at just those locations
        return np.nan_to_num(sample_genos)

    def __calculateImputedGains(self, p): #Private
        """
        Calculate the score change at each p-value
        """
        patient_gains = np.zeros(self.num_samps)
        for i in range(0, len(self.missingness_index[p])):
            ind = self.missingness_index[p][i]
            #NA locations refer to mean snps
            if len(ind[INDEX]) == 0 and len(ind[B]) == 0: #this is our error.
                continue
            else:
                imputed_values = itemgetter(*ind[INDEX])(self.sample_means[p])
                patient_gains[i] = np.dot(imputed_values, ind[B])
                if DEBUG:
                    print(p,i)
                    print("Index of missing snps:",ind[INDEX])
                    print("Corresponding imputed values", imputed_values)
                    print("Correspondinb betas:", ind[B])
        return patient_gains

    def __buildImputedAdjustor(self):
        """
        Function builds a dictionary with vectors to add to scores at each p-value.
        """
        first = True
        prev = None
        for p in self.pvals:
            gains = self.__calculateImputedGains(p)
            if first:
                self.score_modifier[p] = gains
                first = False
            else:
                self.score_modifier[p]  = gains + prev
            prev = self.score_modifier[p]

            
    def updateScores(self, scores_):
        self.__calculateMeans()
        self.__buildImputedAdjustor()
        if(len(scores_) > 1):
            updateLog("Sublisting functionality for missing data has not been implemented yet, will not be accounted for here.")
        for p in self.pvals:
            #TODO add rounding here. Wait until finished testing for other purposes, but do it.
            #scores[p] = np.around(self.score_modifier[p] + scores[p],decimals = 5)
            scores_[p] = self.score_modifier[p] + scores_[p] 
        return scores_          
######################################
class SNPIndexManager:
    def __init__(self, pvals, ss = None): #Default constructor
        self.snp_indices = dict()
        self.var_map = None #maps summary stat variants (index) to genotype/plink file variants (value)
        self.var_map_ids = None #maps variant ids to genotype/plink index
        self.ss_map_ids = None #maps variatn ids to summary stat index
        self.variant_to_pval = dict() #Lists the p-value category of each variant
        self.geno_extraction_map = dict() #for each p values, the indices to get the genotype information out 
        self.geno_subset_extraction_maps = list()
        self.pvals = pvals
        self.ss_counts = dict()
        self.subsets = list() #Store subsets we wish to score separately on, if they exist.
        self.subsets_var_map = list()
        self.var_map_index = { p : dict() for p in self.pvals} #this has the 
        self.extractionListsBuilt = False #Help switch to track if we have built extraction lists 
        if ss is not None:
            self.indexSSByPval(ss) #this updates snp_indices
            for p in self.pvals:
                self.ss_counts[p] = len(self.snp_indices[p][INDEX])
    
    def addSubsets(self, subs):
        if len(subs) == 0 or subs[0] == "":
            ##self.subsets_var_map.append({ p : [] for p in self.pvals})
            print("No subsets of SNPs included")
            return
        self.subsets = subs
        for i in subs:
            #self.subsets.append(i)
            self.subsets_var_map.append({ p : [] for p in self.pvals}) #Add a dictionary for each subset, indexed in order.
        #print(self.subsets)

    def numVariantSubsets(self):
        #print("we've got subsets,", len(self.subsets), "in addition to the full list of variants")
        return len(self.subsets) + 1 #The main subset must be counted, so always at least 1
    def __syncSubsets(self):
        """
        Map each subset to the main dataset- store the indices for all the entries
        """
        if self.var_map is None: #Only if we have already done 
            raise Exception("Variant map has not been set. Program will terminate")
            sys.exit()
        if len(self.subsets) == 0:
            return
        if self.subsets[0] != "":
            self.geno_subset_extraction_maps = [None] * len(self.subsets)
        for i, subset_path in enumerate(self.subsets):
            if subset_path == "":
                continue
            with open(subset_path, 'r') as istream:
                for line in istream:
                    snp = line.strip() 
                    if snp not in self.ss_map_ids:
                        updateLog("Subset SNP " + str(snp) +   " not in full list and will be omitted", False)
                    else:
                        try:
                            pval_cat= self.variant_to_pval[snp] #get the pval group the snp is in
                            index = self.var_map_index[pval_cat][snp] #get the index in the var_map its in
                            #index = self.ss_map_ids[snp]
                            self.subsets_var_map[i][pval_cat].append(index) #this gives the index in the variant map we want for each pvalue
                            #pval:[1,5,1123,55], pval2:[55,12,125]...
                        except:
                            print("Missing the snp, will skip")

    def getSNPCount(self, p):
        try:
            return self.ss_counts[p]
        except IndexError:
            print("Cannot return counts; summary stats have not been initialized")

       #return len(self.snp_indices[p][INDEX])

    def getSSTotals(self):
        return self.ss_counts
        
    def getIndexedBetas(self, p):
        return self.snp_indices[p][B]

    def getAllIndexedBetas(self, p):
        ret_list = [self.snp_indices[p][B]]
        for subset in  self.geno_subset_extraction_maps:
            ret_list.append(subset[p][B])
        return ret_list

    def resetSNPIndex(self):
        self.var_map = None
        self.var_map_ids = None

    #Method for determine if IDs are in the right order
    #Doesn't check every id, just samples a few to improve our confidence that we have the same order
    #By checking >= 5 locations at random and ensuring they are the same. 
    #If each position is equally likely for each entry, then the probability of having different entries by chance
    #For this many checks is really small (i.e. (1/number of snps)^5).
    def __quickIDCheck(self, ss_order, plink_order):
        #Try a number of samples to you get a high certainty
        PROB_ACC = 0.0000001
        import random
        try:
            min_run = (np.log(PROB_ACC) / np.log(1/float(len(ss_order))))
        except ZeroDivisionError:
            print("Error in reading in ids from summary stats >_<")
            sys.exit()
        for i in range(0, int(min_run)+5): #I want to do at least 5. Actually a better way may be to sum up the ids across an interval and see if those are the same
            rand = random.randint(0, len(ss_order))
            try:
                if plink_order[rand][:-2] != ss_order.iloc[rand]:
                    print("SNP order not aligned, realigning now...")
                    return False
            except IndexError:
                #print(plink_order, rand, ss_order)
                print("Note that plink order and ss_order have different lengths")
                return False
        return True

    def mapSStoGeno(self, geno_order, ss_order): #Previously buildVarMap
        """
            Return a numpy array that maps an index in the ss file to a number in plink
            TODO: Speed up this step, its too slow
        """
        #if self.__quickIDCheck(ss_order, geno_order):
        #    self.var_map =  np.array(range(0, len(ss_order) + 1)) 
        #else:
        print("mapping in all cases now. Update when time is right.")
        self.var_map, self.var_map_ids, self.ss_map_ids = self.__hashMapBuild(ss_order, geno_order) 

    def __hashMapBuild(self, ss_snps, plink_snps): #add functionality of a map too
        """
        Create a mapping of SNPs in the SS data to the corresponding SNP in the genotype data.
        @return ret_tap:  numpy array where index corresponds to summary stat variant, and the entry is the index of the plink data genotype
        @return pref: a dictionary of snp ids: index (from the plink data). 
        @r
        """
        ret_tap = np.zeros(len(ss_snps), dtype = np.int32)
        pref = dict()
        ss_pref = dict()
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
            #Modified on march 4- if we are using the subsetting argument, there will be times we fail to find a match
            ss_pref[curr_id] = j
            if curr_id in pref:
                ret_tap[j] = int(pref[curr_id])
            else:
                ret_tap[j] = -1
            """
            TODO: add ability to write those out to another file
            try:
                ret_tap[j] = int(pref[curr_id])
            except KeyError:
                print("Unable to find a match for ", curr_id)
                ret_tap[j] = -1
                updateLog("Found unmatched ID, ", curr_id, ". May be due to SNPs removed that don't pass the call rate (--geno); see mat_form_tmp.log file. If this is not the case, please investigate further.")
                """
        print("Re-indexing complete")
        return ret_tap, pref, ss_pref

    def indexSSByPval(self, ss, method = "split"):
        """
        Generates an index that lists the location of pvalues between a given range
        i.e. betas < p= 0.001 [3,2,5,11]
        i.e. betas from 0.001 < p < 0.01 are at [12,14,15,161]
        i.e. betas from 0.01 < p < 1 are at [1,4,17,18...]
        listed relative to the next smallest SNP, for the way we calculate downstream 
        """
        self.pvals.sort() #smallest one first
        if method != "split": #don't bin the p-values, default is to bin. We could end it here though.
            for p in self.pvals:
                samp = ss[(ss[PVAL] <= float(p))]
                self.snp_indices[p] = [samp.index.values.tolist(), samp[BETA].values] 
                self.variant_to_pval.update(zip(samp[REFID].tolist(),[p]*samp.shape[0]))
     
        prev = -1
        for p in self.pvals:
            samp = ss[(ss[PVAL] <= float(p)) & (ss[PVAL] > prev)]
            self.snp_indices[p] = [samp.index.tolist(), samp[BETA].values]
            #Keep track of where in the list each variant is 
            self.var_map_index[p] = dict(zip(samp[REFID], list(range(0, len(samp))))) #Each Id has its location
            prev = float(p) 
            self.variant_to_pval.update(zip(samp[REFID].tolist(),[p]*samp.shape[0]))

    def getSumStatSNPIndex(self, p):
        return self.snp_indices[p][INDEX]

    def __buildExtractionList(self):
        ss_snps = dict()
        for p in self.pvals:
            ss_snps[p] = self.snp_indices[p][INDEX]
            if len(ss_snps[p]) > 0:
                self.geno_extraction_map[p] = (itemgetter(*ss_snps[p])(self.var_map))
            else:
                self.geno_extraction_map[p] = []

    def getGenoIndices(self, p): #formerly buildExtractionlist
        if not self.geno_extraction_map: #its empty
            self.__buildExtractionList()
            self.__syncSubsets()
        return self.geno_extraction_map[p]

    def __buildSubsetExtractionLists(self):
        """
        Update the subset lists to contain the relevant extracted values for scoring along with the Betas.
        """
        if len(self.subsets_var_map) > 0:
            print("Building the subset extraction lists....")
        for i,subset in enumerate(self.subsets_var_map): #ordered subsets.
            self.geno_subset_extraction_maps[i] = { p : [[],[]] for p in self.pvals} #make each a tuple, with betas in one and 
            for p in self.pvals:
                curr_p_indices = subset[p]
                #Make sure its not empty
                #POSSIBLE ERROR SORUCE HERE.
                if len(curr_p_indices) > 0 and len(self.snp_indices[p][INDEX]) > 0:
                    converted_indices = ((itemgetter(*curr_p_indices)(self.snp_indices[p][INDEX])))
                    converted_betas = ((itemgetter(*curr_p_indices)(self.snp_indices[p][B])))
                    if len(converted_indices) > 0:
                        ret_app = (itemgetter(*converted_indices)(self.var_map))
                    else:
                        ret_app =  []
                else:
                    ret_app = []
                    converted_betas = []

                self.geno_subset_extraction_maps[i][p][INDEX] = ret_app
                self.geno_subset_extraction_maps[i][p][B] = converted_betas
        return

    def getAllGenoIndices(self, p):
        #if not self.geno_extraction_map or not self.geno_subset_extraction_maps: #its empty
        if not self.extractionListsBuilt:
            self.__buildExtractionList()
            self.__syncSubsets()
            self.__buildSubsetExtractionLists()
            self.extractionListsBuilt = True
            print("Variant extraction lists made")
        ret_list = [self.geno_extraction_map[p]]
        for subset in self.geno_subset_extraction_maps: #ordered subsets.
            ret_list.append(subset[p][INDEX])
        return ret_list
######################################

def memoryAllocation(set_as = ""):
    global MEMORY   
    if set_as != "":
        MEMORY = " --memory " + str(set_as)
    return MEMORY

def updateLog(*prints):
    """
    Update the log file, can write out to console as well if last argument passed is true
    @param prints is list of things to write to the log on separate lines.
    """
    with open(LOG, 'a') as ostream:
        if START_SEQ == prints[0]:
            ostream.write("\n---------------- New BPRS Run ----------------\n")
        if prints[-1] == True:
            for i in range(0, len(prints)-1):
                print(prints[i])
                ostream.write(str(prints[i]) + ' ')
        else:
            for a in prints:
                if a != False:
                    ostream.write(str(a) + " ")
        ostream.write("\n")
    
def getEntryCount(fin, header):
    """
    Quickly count the number of entries in a file
    @param fin: file to count
    @param header: if the file has a header not.
    """
    if header:
        return sum(1 for line in open(fin)) -1
    else:
        return sum(1 for line in open(fin))

def firstColCopy(f1, f2):
    """
    Write the first column of a file to another file, using bash commands
    Typically used for ss_filt > geno_ids
    """
    command = " awk '{print $1}' " + f1 + " > " + f2
    check_call(command, shell = True)

def getPatientCount(snp_file, vers):
    try:
        if vers == 1:
            return getEntryCount(snp_file + ".fam", False)
        return getEntryCount(snp_file + ".psam", True)
    except:
        print("It appears the plink family file is missing or has the incorrect file extension. Please check it out!")
        sys.exit()

def setFtype(plink_version):
    if plink_version == 2:
        return " --pfile ", ".pvar"
    else:
        return " --bfile ", '.bim'

def readInSummaryStats(s_path):
    #ID, ref, alt, beta, pvalue. Th
    #BEWARE: may have issue with sep, depending..
    try:
        ss = pd.read_csv(s_path, sep = '\t', header = None,  names = [REFID, REF,ALT, BETA, PVAL], dtype= { REF:'category',ALT: 'category', BETA:'float64', PVAL:'float64'})
    except TypeError:
        print("It seems that there was an error with the summary statistics data. Did you omit --reset_ids?")
        print("Also try clearing out all the *.f files in the local directory.")
        sys.exit()
    print("Data in memory...")
    return ss


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
    try:
        if vplink == 2:
            local_geno = snp_file + ".pvar"
            updateLog("Number of variants detected in original genotype file", str(getEntryCount(snp_file + ".pvar", True)))
        else: #vplink == 1:
            local_geno = snp_file + ".bim"
            updateLog("Number of variants detected in original genotype file", str(getEntryCount(local_geno, False)))
    except FileNotFoundError:
        print("Are you sure you specified the correct version of plink? We are unable to find the plink file you specified...")
        print("Program will quit.")
        sys.exit() 
    geno_id_list = "geno_ids.f"
    inter_sum_stats = "ss_filt.f" #Intersected summary stats with genotype ids
    #If we have already generated our geno_id_list,skip this step:
        #Logic as follows: if the geno_ids.f file doesn't exist, or the variant IDs haven't been set to 1:123:R:A or the geno_ids file has size 0, start from scratch 
    if VAR_IN_ID and os.path.isfile("local_geno.pvar") and os.stat("local_geno.pvar") != 0 and args.preprocessing_done:
            print("Using the existing pvar file...")
            local_geno = "local_geno.pvar" 
    if not os.path.isfile(geno_id_list) or not VAR_IN_ID or os.stat(geno_id_list).st_size == 0 : #If we've already done this step, no point in doing more work....   
        if not VAR_IN_ID: #need to rework the IDs so they are chr:pos:ref:alt
            local_geno = "local_geno.pvar" #We will be making our own pvar
            command = ''' awk '(!/##/ && /#/) {print $1"\t"$2"\tID\t"$4"\t"$5} (!/#/) {print $1"\t"$2"\t"$1":"$2":"$4":"$5"\t"$4"\t"$5}' ''' + snp_file + ".pvar > " + local_geno
            if vplink == 1:
                #local_geno = "local_geno.bim"
                #command = '''awk '{print $1"\t"$1":"$4":"$5":"$6"\t"$3"\t"$4"\t"$5"\t"$6}' ''' + snp_file + ".bim > " + local_geno
                print("Following filtering, this utility requires files in plink2 format. Please convert your plink1 data to plink2 to proceed")
                sys.exit()
            try:
                check_call(command, shell = True)
            except subprocess.CalledProcessError:
                print("Unable to generate modified pvar/bim file. Are you sure you specified the correct plink version?")
                updateLog("Failed on command", command, True)
                sys.exit()
        if VAR_IN_ID and os.path.isfile("local_geno.pvar") and os.stat("local_geno.pvar") != 0:
            print("Using the existing pvar file...")
            local_geno = "local_geno.pvar" 
        if vplink == 2:
            command_n = "awk ' (NR> 1 && $1 !~ /X/) {print $3}' " + local_geno + " > " + geno_id_list
        else:
            command_n = "awk ' ($1 !~ /X/) {print $2}' " + local_geno + " > " + geno_id_list
        check_call(command_n, shell = True)
        VAR_IN_ID = True#we have imposed the variant information into the IDs
    #If the ss_filt.f file doesn't exist or it has a different length than the geno_ids file,
    if not os.path.isfile(inter_sum_stats) or (os.path.isfile(inter_sum_stats) and (getEntryCount(geno_id_list, False) != getEntryCount(inter_sum_stats, False))):
        #Choose only ids below the highes p-value threshold that appear both in the genotype and summary stat data. Trying to minimizize write out time later on.
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
        firstColCopy(inter_sum_stats, geno_id_list)
    return local_geno, geno_id_list, inter_sum_stats

def scoreCalculation(geno, betas):
    if len(geno) < 2:
        print("an empty list")
        return 0
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
        print(len(geno))
        print(len(betas))
        print("It appears that the number of SNPs is not as expected. We suggest lowering the missing genotype filter threshold, and additional debugging.")   

def recordSNPsWeights(curr_p,im, dat, SNPS_):
    """
    Helper function to just update the DEBUG data with weights and such
    @param curr_p: pvalue on current iteration
    @param im: index manager object
    @param dat: the line of data we are considering
    @param SNPS_: index of SNP data onward.
    """
    DEBUG_DAT[BETA] = dict()
    DEBUG_DAT[SNP] = dict()
    prev_append_snp = []
    prev_append_beta = []
    for p in pvals:
        DEBUG_DAT[BETA][curr_p] = im.getIndexedBetas(curr_p)
        if len(im.getAllGenoIndices(p)[0]) == 0:
            DEBUG_DAT[SNP][curr_p] = []
        else:
            #sel = (itemgetter(*snp_index)(var_map))  #Index manager
            #sel = index_manager.getGenoIndices(p)
            sel = im.getAllGenoIndices(curr_p)[0] #I would prefer to keep it the old way, but for working sake.....
            try:
                DEBUG_DAT[SNP][curr_p] = ((itemgetter(*sel)(dat[SNPS_:])))
            except TypeError as e:
                print(sel)
                print(e)
                print("SNP ID", SNP)
                print("Pvalue", curr_p)
                print("Errors detected in matrix header. Will not print snps used.")
                continue
        DEBUG_DAT[BETA][curr_p] = np.concatenate((DEBUG_DAT[BETA][curr_p],prev_append_beta))
        DEBUG_DAT[SNP][curr_p] = np.concatenate((DEBUG_DAT[SNP][curr_p], prev_append_snp))
        prev_append_snp = DEBUG_DAT[SNP][curr_p]
        prev_append_beta = DEBUG_DAT[BETA][curr_p]

def linearGenoParse(snp_matrix, index_manager,pvals, scores, ss_ids, nsamps,args, writeSNPListsOnly = True):
    """
    Parse the patient file
    @param snp_matrix -- the path to the massive genotype file
    @param snp_indieces -- the map with the indices for each pvalue we care about
    @param scores -- what is returned, stores the scores
    @param pvals -- which pvalues we measure at
    @param ss_ids -- the list of ids from the summary stats list, for ordering purposes
    @param writeSNPLists - basically a dummy variable to .
    """
    #############Set up data
    no_nas = args.no_na
    patient_ids = list()
    report_index = 100
    #Plink (matrix) data index for readability
    PDI = {"FID":0, "IID":1, "PAT":2, "MAT":3, "SEX":4, "PHENOTYPE" : 5, "SNPS":6}
    SNPS_ = PDI["SNPS"]
    imputers = list()
    for i in range(0, index_manager.numVariantSubsets()):
        imputers.append(MissingGenoHandler(index_manager.getSSTotals(),nsamps,pvals)) #create the object to handle this....
    
    ############Parse the matrix file
    print("Opening", snp_matrix)
    import mmap; import contextlib
    na_set = set(["NA", "-", "na", "NaN"])
    with open(snp_matrix + ".raw", 'r') as f:
        with contextlib.closing(mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)) as istream:
            line_counter = 0
            dat = istream.readline().strip()
            pvals.sort()
            while dat:
                sample_counter = line_counter - 1 #Track which sample we are on, 1 - the line number
                dat = (dat.decode('ascii')).split('\t')
                if line_counter == 0:
                    index_manager.mapSStoGeno(dat[SNPS_:], ss_ids) #builds a map of variants
                    DEBUG_DAT[REFID] = dat[SNPS_:] #track the variants for debug output
                    print("Proceeding with line-by-line calculation now...")
                    ########### Update the log, write out SNPs and weights
                    prev = 0
                    for p in pvals:
                        num_snps = index_manager.getSNPCount(p)
                        updateLog("Number of final SNPs used at ", str(p), "pvalue level:", str(num_snps + prev))
                        prev = prev + num_snps
                        recordSNPsWeights(p,index_manager, dat, SNPS_)
                    writeScoresDebug(DEBUG_DAT[SNP], os.path.basename(snp_matrix) + ".debug_vals_snps.tsv")
                    writeScoresDebug(DEBUG_DAT[BETA], os.path.basename(snp_matrix) + ".debug_vals_betas.tsv")
                    print('SNP order established. Scoring will now begin.')        
                else:
                    rel_data = dat[SNPS_:]
                    patient_id = str(dat[PDI["IID"]])
                    patient_ids.append(patient_id)
                    prev_score = 0; new_score = 0
                    for j, p in enumerate(pvals):
                        i = 0
                        for sel, beta in zip(index_manager.getAllGenoIndices(p), index_manager.getAllIndexedBetas(p)):
                            prev_score = 0
                            if i == 0 and args.subsets_only: #we don't want the full list, just the sublist 
                                i += 1; continue
                            if j > 0:
                                prev_p = pvals[j-1]
                                prev_score = scores[i][prev_p][sample_counter]
                            if len(sel) == 0: #There are 0 snps at this significance level (!) #use the manager
                                new_score = prev_score + 0
                                scores[i][p][sample_counter] = new_score
                                if DEBUG: updateLog(str(patient_ids[-1]), str(new_score))
                                continue
                            try:
                                rel_genos  = itemgetter(*sel)(rel_data) #make sure this is the best way to do it.
                                if no_nas:
                                    new_genos = np.array(rel_genos, dtype = np.float64)
                                else:
                                    if any([x in rel_data for x in na_set]):
                                        new_genos = np.genfromtxt(rel_genos, dtype = np.float64)
                                    else:
                                        new_genos = np.array(rel_genos, dtype = np.float64)
                                    if len(index_manager.getAllGenoIndices(p)) > 1: #If we are snp subsetting.
                                        print("Please ensure genotype SNP file has no missing data before proceeding. BPRS is unable to handle snp sublists and NAs in matrix data.")
                                        sys.exit()
                                    new_genos = imputers[i].handleNAs(new_genos, beta, p)
                            except IndexError:
                                updateLog("Patient ", patient_id, "appears to have missing entries and will not be scored at",p,"threshold", True)
                                #Find the NA
                                scores[i][p][sample_counter] = float('nan')
                                continue #don't calculate the score, simply proceed to the next iteration of p. 
                            new_score = scoreCalculation(new_genos, beta) + prev_score
                            try:
                                scores[i][p][sample_counter] = new_score
                            except IndexError:
                                print("Subset number", i)
                                print("pvalue", p)
                                print("Sample count", sample_counter)
                                print(len(scores[i][p][sample_counter]))
                            if DEBUG: updateLog(str(patient_ids[-1]), str(new_score))
                            i += 1
                dat = istream.readline()
                line_counter += 1
                if sample_counter % report_index == 0:
                    print("Currently at individual/sample", sample_counter)
    print("Finished parsing sample genotype data!")
    if not no_nas:
        print("Mean imputing missing genotype data, if any...")
        for i, imp in enumerate(imputers):
            scores[i] = imp.updateScores(scores[i])
    updateLog("Number of genes with missing entries detected in genotype data: " +str(imputers[0].getNACount()), True)
    if sample_counter <= 1:
        updateLog("There appears to be some error with the genotype plink matrix- no entries were found. Please ensure the genotype data is correct or the plink call was not interrupted.  Consider deleting *.raw and re-running.", True)
        sys.exit()     
    return scores, patient_ids #basically a pval:[list of snps]  

def calculatePRS(pvals, snp_matrix, index_manager, snp_list, sample_count, args):
    """
    Actually calculates the PRS Scores, and times the process
    Returns a list of scores and the corresponding patient_ids, as well as any relevant debug informtion if specified
    @param snp_matrix- the path to the enormous snp matrix
    @param snp_indices: the information of snp_indices at each p-value
    @param snp_list: a series object containing the name of each snp
    """
    #Seraj shoutout
    start = time.time()
    #scores = [{ p : ([0]*sample_count) for p in pvals}]*index_manager.numVariantSubsets()
    #Doesn't quite work, point to same things.
    scores = list()
    for i in range(0, index_manager.numVariantSubsets()):
        scores.append({ p : ([0]*sample_count) for p in pvals}) #Error here, haven't specifid the number of samples yet.
    scores, patient_ids = linearGenoParse(snp_matrix, index_manager,pvals, scores, snp_list, sample_count, args)
    end = time.time()
    print("Calculation and read in time:", str(end-start))
    return scores, patient_ids

def buildAlignHelper(ss, ss_type):
    """
    Create the "aligner.t" file that will be used to make sure the reference and target data line up.
    """
    if ss_type == "SAIGE":
            first = '''awk '(NR >1) {print $1":"$2"\t"$4} ' ''' + ss + " | uniq > aligner.t"
            check_call(first, shell = True)
    elif ss_type == "NEAL":
        #Separate the major allele from the other info
        first = "cut -f 1 " + ss + " | tail -n +2 | cut -f 1,2 -d ':' > ids.t"
        second = "cut -f 1 " + ss + " | tail -n +2 | cut -f 3 -d ':' > ref.t"
        final = "paste ids.t ref.t | uniq > aligner.t"
        check_call(first, shell = True)
        check_call(second, shell = True)
        check_call(final, shell = True)
    elif ss_type == "DEFAULT": #ID REf ALT BETA Pvalue
        call = "cut -f 1 " + ss + ''' | awk -F ":" '(NR > 1) {print $1":"$2"\t"$3}' | uniq > ./aligner.t'''
        check_call(call, shell = True)
    else:
        print("Have not extended the method for type", ss_type, ". Please align Reference information manually")
        sys.exit()
    return "aligner.t"

def plinkRemoveDupsHelper(old_plink, plink_version):
    """
    Create the "remove_multi.t" file that will be used to make sure no multi-allelic snps from target data are considered. 
    """
    ftype, ext = setFtype(plink_version)
    #remove_multi = "cut -f 3 " + old_plink + ext + " | uniq -d > remove_multi.t"
    remove_multi = '''awk  '/^[^#]/ { print $1":"$2}' ''' + old_plink + ext + " | uniq -d > remove_multi.t"
    if plink_version == 1:
        remove_multi = '''awk  '/^[^#]/ { print $1":"$4}' ''' + old_plink + ext + " | uniq -d > remove_multi.t"
    check_call(remove_multi, shell = True)
    return "remove_multi.t"

def plinkUpdateIDsHelper(old_plink, plink_version):
    """
    Create the std_ids.t file for resetting the input ids for target data.
    
    ext = ".pvar"
    id_comm = '''awk  '/^[^#]/ { print $3"\t"$1":"$2}' ''' + old_plink + ext + " > std_ids.t"
    if plink_version == 1:
        ext = ".bim"
        id_comm = '''awk  '/^[^#]/ { print $2"\t"$1":"$4}' ''' + old_plink + ext + " > std_ids.t"
    check_call(id_comm, shell = True)
    return "std_ids.t"
    """
    #New version: replicate the pvar file fully, simply with new IDs.
    ext = ".pvar"
    id_comm = '''awk  '(!/##/) { print $1"\t"$2"\t"$1":"$2"\t"$4"\t"$5"\t"$6"\t"$7}' ''' + old_plink + ext + ''' | sed "1s/.*/#CHROM\tPOS\tID\tREF\tALT\tFILTER\tINFO/" > std_ids.t'''
    if plink_version == 1:
        ext = ".bim"
        #id_comm = '''awk  '/^[^#]/ { print $2"\t"$1":"$4}' ''' + old_plink + ext + " > std_ids.t"
        id_comm = '''awk  '(!/##/) { print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' ''' + old_plink + ext + '''  > std_ids.t'''
    check_call(id_comm, shell = True)
    return "std_ids.t"

def alignReferenceByPlink(old_plink,plink_version, ss, ss_type, plink_path):
    """
    Make the reference file listing <IDs> <REF>
    Call plink to swap them
    Write this out to a new_plink file
    """
    print("Correcting for strand flips")
    #Remove duplicates
    remove_file = plinkRemoveDupsHelper(old_plink, plink_version)
    new_ids_file = plinkUpdateIDsHelper(old_plink, plink_version)
    ftype, ext = setFtype(plink_version)
    
        #In order to align, we need to make sure we have unique ids- so we account for multi-allelic SNPs and indels here.
        #Filter out multi-allelic snps, X-chr snps,ilter variants with missing call rates > 0.05%, and recode the names to a standard.
    try:
        plink_call_first = (plink_path + "plink2 " + ftype + old_plink +  " --pvar " + new_ids_file +  " --geno 0.05 --make-pgen --out filtered_target --exclude " + remove_file + " --not-chr X" + memoryAllocation())
        check_call(plink_call_first, shell = True)
        #" --update-ids " + new_ids_file + <modified the pvar directly
    except subprocess.CalledProcessError:
        updateLog("We have encountered an error calling plink2. Are all of your file paths corrected? Do you have plink2 installed. It can be found at https://www.cog-genomics.org/plink/2.0/", True)
        updateLog(plink_call_first)
        sys.exit() 

    #Build the aligner file
    aligner_ref = buildAlignHelper(ss, ss_type)
    plink_call = plink_path + "plink2 --pfile filtered_target --ref-allele force " + aligner_ref + " --make-pgen --out new_plink" + memoryAllocation()
    try:
        check_call(plink_call, shell = True)
    except:
        print("It appears there was an error in aligning the reference, likely due to known alleles that we attempted to swap. We suggest correcting these alleles or aligining the reference strands indepndently using plink --ref-allele")
        updateLog("Command yielding error:" + plink_call)
        sys.exit()
    if not DEBUG:
        if os.path.isfile("filtered_target.pvar"):
            check_call("rm filtered_target*", shell = True) 
        clean_up = "rm *.t"
        check_call(clean_up, shell = True)
    return "new_plink"

def filterSumStatSNPs(ss, ss_type, keep_ambig_filter):
    """
    Filters summary statistics based on the specified file type, as follows:
    1) Remove anything more than bi-allelic
    2) Limit it to SNPs only
    3) Limit to nonambiguous SNPs, if possible keep ambiguous snps MAF > 0.3
    4) Remove null beta values
    5) Remove X chromosome snps
    6) Write it out in the default format so its easier to digest downstream. (ID REF ALT BETA PVAL), ID standardized to chr:pos:ref:alt
    """
    filter_list = {"AG", "AC", "AT", "TA", "TC", "TG", "GC", "GT", "GA", "CG", "CT", "CA"}
    noambig_list =   {"AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG"}
    ambig_list = {"AT", "TA", "GC", "CG"}
    sizes = {"start": 0,"non_bi": 0, "ambig":0, "na":0, "x":0}
    comment = {"ambig": "Ambiguous SNPs along with INDELS removed:", "na": "Variants with missing effect sizes removed:", "x": "Genes on the X chromosome removed:", "non_bi": "Non-biallelic SNPs removed:"}
    #Set this to 0.3 from 0.4 on MArch 16
    ambig_filter = 0.3
    if keep_ambig_filter:
        ambig_filter = -1
    ss_t = None
    out_name = "filtered_" + ss_type + ".ss"
    #Some room to beautify the code here...
    if ss_type == "HELIX":
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
        #effect is always givein in terms of the alternate allele.  
        #variant    minor_allele    minor_AF    low_confidence_variant  n_complete_samples  AC  ytx beta    se  tstat   pval
        ss_t = pd.read_csv(ss, sep = '\t',usecols = ["variant","minor_allele", "beta", "pval", "minor_AF"], dtype = {"variant": 'string_', "beta":'float64', "pval":'float64', "minor_AF":"float64"})
        sizes["start"] = ss_t.shape[0]
        id_data = ss_t.variant.str.split(":", expand = True)  #gives chr:pos:ref:alt
        ss_t["comb"] = id_data[2] + id_data[3]#last 2 characters at the back end of the variant
        ss_t["CHR"] = id_data[0]
        ss_t["location"] = id_data[0] + ":" + id_data[1]
        ss_t[REF] = id_data[2]
        ss_t[ALT] = id_data[3]
        ss_t.drop_duplicates(subset = "location", keep = False, inplace = True)
        sizes["non_bi"] = ss_t.shape[0]
        #Removing ambiguous SNPs
        ss_t = ss_t[ss_t.comb.isin(filter_list)]
        #March 16, switch from > to <
        ss_t = ss_t[(ss_t.comb.isin(ambig_list) & (ss_t.minor_AF < ambig_filter)) | (ss_t.comb.isin(noambig_list))]
        #remove ones that are in the ambiguous list unless they have a MAF < 0.3, our default
        sizes["ambig"] = ss_t.shape[0]
        #Also remove SNPs with NA
        ss_t = ss_t[ss_t.beta.notnull()] #here its NaN
        sizes["na"] = ss_t.shape[0]
        #Also remove X chromosome
        ss_t = ss_t[ss_t.CHR != "X"]
        sizes["x"] = ss_t.shape[0]
        #Set it up for write_out ID, ref, alt, beta, pval
        ss_t.rename(columns = {"variant":"ID", "beta":BETA, "pval":PVAL}, inplace = True)
        ss_t["ID"] = ss_t.location + ":" +  ss_t.REF + ":" + ss_t.ALT
        ss_t = ss_t[["ID", REF, ALT, BETA,PVAL]]
    elif ss_type == "DEFAULT":
        #[ID REF ALT beta pval]
        ss_t = pd.read_csv(ss, sep = '\t')
        sizes["start"] = ss_t.shape[0]
        try:
            ss_t["comb"] = ss_t[REF] + ss_t[ALT]
        except KeyError:
            print("It seems that the input summary stats file of type DEFAULT doesn't have the right header or delimiter.")
            sys.exit()
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

def plinkToMatrix(snp_keep, args, local_pvar, vplink, plink_path):
    """
    Create a matrix with scores we can parse
    @param snp_keep is the path to the file with the SNP list. Currently just a single column, TODO check if this will work for plink filtering.
    """
    #Write out the list of SNPs to keep
    #Careful- just to check our SNP numbers
    if os.path.isfile("mat_form_tmp.raw"):
        #Check that its SNPs match the SNPs we have
        with open("mat_form_tmp.raw", 'r') as istream:
            header_line = istream.readline().strip()
        if len(header_line) - 6 == getEntryCount(snp_keep, False):
            print("Using currently existing genotype matrix...")
            return "mat_form_tmp"

    syscall = " --pfile "
    annot =  " --pvar "
    if vplink == 1:
        syscall = " --bfile "
        annot = " --bim "
    call = plink_path + "plink2 " + syscall  + args + annot + local_pvar + " --geno 0.05 --not-chr X --extract " + snp_keep + " --out mat_form_tmp --export A" + memoryAllocation()
    check_call(call, shell = True)
    #--geno removes snps with more than % missing from the data.
    print("User-readable genotype matrix generated")
    if not os.path.isfile("mat_form_tmp.raw"):
        print("PLINK file was not correctly created. Program will terminate.")
        sys.exit()
    return "mat_form_tmp"

#TODO: fix the id syst" + clumping, make sure consisent
def plinkClump(reference_ld, clump_ref,clump_r, maf_thresh, geno_ids, ss):
    """
    Run plink and clump, get the new list of ids.
    Select just these from the summary statistics.
    """
    pre_clump_count = getEntryCount(geno_ids, False)
    updateLog("Clumping...", True)
    updateLog("SNPs prior to clumping", str(pre_clump_count))
    if clump_ref == "NA":
        #Build the file for plink to run on 
        command = "echo SNP P > clump_ids.tmp"
        check_call(command, shell = True)
        command = 'cut -f 1 ' + ss + " | cut -f 1,2 -d ':' > only_ids.tmp"
        check_call(command, shell = True)
        command = 'cut -f 5 ' + ss + " > only_ps.tmp && paste only_ids.tmp only_ps.tmp >> clump_ids.tmp" 
         #Gives SNP id and p-value to assit in clumping
        check_call(command, shell = True)
        #Run plink clumping,do not use --clump-best
        plink_command = "plink --bfile " + reference_ld + " --clump clump_ids.tmp --clump-p1 1 --clump-p2 1 --clump-r2 " + str(clump_r) + " --clump-kb 250 --clump-field P --clump-snp-field SNP --maf " + str(maf_thresh) + memoryAllocation()
        updateLog(plink_command)
        check_call(plink_command, shell = True)
        clump_file = "plink.clumped" #The default output name with flag clump
    else:
        print("Using provided reference clump...")
        clump_file = clump_ref
    #Get a way to map from one to the other
    temp_ids = '''cut -f 1,2 -d ":" ''' + ss + " > t && paste  t " + ss + " > id_mapper.tmp && rm t"
    check_call(temp_ids, shell = True)
    command = "awk '(FNR == NR) {a[$3];next} ($1 in a) {print $0}' " + clump_file + " id_mapper.tmp | cut -f 2,3,4,5,6  > t && mv t "+ ss
    check_call(command, shell = True)
    if not DEBUG:
        command = "rm *.tmp"
        check_call(command, shell = True)
    #Update geno_ids for extraction.
    firstColCopy(ss, geno_ids)
    #Log reporting
    post_clump_count = getEntryCount(geno_ids, False)
    updateLog("Number of SNPs removed by clumping (selecting only top clump variant):", str(pre_clump_count - post_clump_count)) 
    updateLog("SNPs remaining:", str(post_clump_count))
    return ss, geno_ids

def writeScores(scores, ids,args): #subset_files,dest_path):
#outfile_name = args.output + os.path.splitext(matrices_for_parsing[i])[0] + ".scores.tsv"
    #Print the full one:
    subset_files = args.split_scores
    dest_path = args.output
    if not args.subsets_only:
        full = scores[0]
        full["IID"] = ids
        tab_out = pd.DataFrame.from_dict(full)
        tab_out = tab_out.set_index("IID")
        tab_out = tab_out.round(5)
        tab_out.to_csv(dest_path + '.full_prs.scores.tsv', sep = '\t') 

    subset_list = subset_files.split(',')
    if (len(subset_list) != len(scores)) and (len(subset_list) + 1 != len(scores)): #Subset_lists will always have 1, because its an empty string
        print("There was an error in scoring on subsets, scores were not recorded")
        print(len(subset_list))
        print(len(scores))
        sys.exit()
    c = 1
    for subset_name in subset_list:
        if subset_name == "":
            continue
        subset = scores[c] #we already did the base one, with a different name.
        subset["IID"] = ids
        tab_out = pd.DataFrame.from_dict(subset)
        tab_out = tab_out.set_index("IID")
        tab_out = tab_out.round(5)
        a = os.path.basename(subset_name)
        a = os.path.splitext(a)[0]
        destination = dest_path + "." + a + ".subset_prs." + "scores.tsv"
        tab_out.to_csv(destination, sep = '\t') 
        c += 1
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Basic tools for calculating rudimentary Polygenic Risk Scores (PRSs). This uses PLINK bed/bim/fam files and GWAS summary stats as standard inputs. Current version requires both plink2 and plink1")
    parser.add_argument("-snps", "--plink_snps", help = "Plink file handle for actual SNP data (omit the .extension portion)")
    parser.add_argument("--preprocessing_done", action = "store_true", help = "Specify this if you have already run the first steps and generated the plink2 matrix with matching IDs. Mostly for development purposes.", required = "--split_scores" in sys.argv)
    parser.add_argument("-ss", "--sum_stats", help = "Path to summary stats file. Be sure to specify format if its not DEFAULT. Assumes tab delimited", required = True)
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE", "NEAL"])
    parser.add_argument("--pvals", default = "ALL", help= "Specify which p-values to threshold at. To do multiple at once, list them separated by commas")
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--prefiltered_ss", default = False, action = "store_true", help = "Specify this option if you have already filter summary stats for ambiguous SNPs, bi-alleleic variants, NAs. Recommended for speed.")
    parser.add_argument("--reset_ids", action = "store_true", help = "By default, program assumes genotype ids are in the form chr:loc:ref:alt. Specify this argument if this is NOT the case to.")
    parser.add_argument("--clump", help = "Specify if you wish to do clumping. Option to pass in a plink clumping file or the reference panel plink data.", action = "store_true", required = "--ld_ref" in sys.argv or "--clump_ref" in sys.argv)
    parser.add_argument("--ld_ref", default = None, help = "Specify an LD reference panel, required if you want to do clumping. Pass in plink1.9 format please, and ensure that all SNP ids have the form chr:pos (i.e. 1:2345).")
    parser.add_argument("--clump_r2", default = 0.25, help = "Specify the clumping threshold.")
    parser.add_argument("--debug", help = "Specify this if you wish to output debug results", default = False, action = "store_true")
    parser.add_argument("--clump_ref", help = "Specify this option if you wish to perform clumping and already have a plink clump file produced (the output of plink1.9 clumping). Provide the path to this file.", default = "NA")
    parser.add_argument("-pv", "--plink_version", default = 2, help = "Specify which version of plink files you are using for the reference data.", type = int, choices = [1,2])
    parser.add_argument("--align_ref", action = "store_true", help = "Specify this if you would like us to align the reference alleles between the target and training data. If this is already done, we recommend omitting this flag since it will create an additional copy of the genotype data.")
    parser.add_argument("--preprocess_only", action = "store_true", help = "Specify this option if you would like to terminate once preprocessing of SNPs (including aligning reference, removing ambiguous SNPs, preparing summary stat sublist) has completed.")
    parser.add_argument("--plink_snp_matrix", help = "Pass in a previously calculated SNP matrixs")
    parser.add_argument("--no_ambiguous_snps", help = "Specify this if you wish to remove all ambiguous SNPs whatsover, regardless of MAF. Default is to keep those with MAF < 0.3", action = "store_true", default = False)
    parser.add_argument("--memory", help = "Specify memory allocation in MB for tasks called in plink. Default omits this argument", default = "")
    parser.add_argument("--plink2_path", help = "Specify the path to plink2 if its not in your PATH variable", default = "") 
    parser.add_argument("--split_scores", help = "Specify (a) list(s) of variants to split scores by after filtering has been done already.P;ease list separated by commas (no spaces.) Importantly, this happens after clumping: for other list inputs, do and then clump. Typically used when matrix has already been calculated", default = "")
    parser.add_argument("--serialize", help = "If you want to write out scores debugging purposes", action = "store_true")
    parser.add_argument("--no_na", action = "store_true", help = "Specify this if you know your data contains no NAs.", default = False, required = "--split_scores" in sys.argv)
    parser.add_argument("--debug_pickle", help = "Specify a pickled score data if its available. For debugging.")
    #TODO implement the path for subsetting by an input list with the following argument
    parser.add_argument("--select_vars", help = 'Specify a list of variants that you want in the analysis NOT YET IMPLENENTED')
    parser.add_argument("--subsets_only", help = "Score only on subsets indicated, not on the full list", action = "store_true", default = False)
    args = parser.parse_args()
    DEBUG = args.debug
    start = time.time()    
    #Extract the pvalues to asses at
    updateLog(START_SEQ, str(datetime.datetime.now()), '\n', "------------------------------")
    args_print = str(args).split(",")[1:]
    updateLog("Arguments", '\n'.join(args_print))
    if args.pvals != "ALL":
        t = (args.pvals).strip().split(",")
    else:
        t=["5e-8", "1e-6", "1e-5", "1e-4", "1e-3", "0.01", "0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "1"]
    pvals = [float(x) for x in t]
    updateLog("Pvals to asses:", str(",".join(t)))
    memoryAllocation(set_as = args.memory)
    if args.reset_ids:
        VAR_IN_ID = False #Variant information is not included in genotype id
    if args.preprocessing_done:
        if not args.clump:
            print("Please note that you have not provided clumping data, so no clumping will be performed.")
        updateLog("Assuming data preprocessing- including filtering by summary stats, aligning references, and clumping- has already been completed. Proceeding with PRS calculations...", True)
    if not args.prefiltered_ss and not args.preprocessing_done:
        print("Filtering ambiguous SNPs and indels from SS data...")
        args.sum_stats = filterSumStatSNPs(args.sum_stats, args.ss_format, args.no_ambiguous_snps)
        args.ss_format = "DEFAULT" #ID REF ALT BETA PVAL
        updateLog("Ambiguous SNPs and indels removed. Filtered summary statistics written out to", args.sum_stats, True)
    if args.align_ref and not args.preprocessing_done:
        args.plink_snps = alignReferenceByPlink(args.plink_snps, args.plink_version, args.sum_stats, args.ss_format, args.plink2_path) #Make a local copy of the data with the updated genotype
        print("A new plink2 fileset with matched reference sequences has been generated", args.plink_snps)
        updateLog("Strand flips have been corrected for. Proceeding with analysis", True)
        args.plink_version = 2 #We have updated the type to 2
    if not args.align_ref:
        print("Genotype ref/alt SNPs will not be aligned to GWAS ref/alt SNPs. Please use the --align_ref to align them")
        
    print("Selecting IDs for analysis...")
    local_pvar, geno_ids, ss_parse = prepSNPIDs(args.plink_snps, args.sum_stats,args.ss_format, vplink = args.plink_version, max_p = max(pvals))
    num_pat = getPatientCount(args.plink_snps, args.plink_version)
    updateLog("Number of patients detected in sample:", str(num_pat), True)
    updateLog("Time for preprocessing (SNP filtering, reference alignment if specified, etc.):", str(time.time() - start), True)


    if args.clump:
        ss_parse, geno_ids = plinkClump(args.ld_ref, args.clump_ref,args.clump_r2, args.maf, geno_ids, ss_parse)
    print("Preprocessing complete")
    if args.preprocess_only:
        updateLog("Plink file data has been updated into new_plink.*. Use this in downstream runs.", True)
        sys.exit()

    stats_complete = readInSummaryStats(ss_parse)
    snp_index_manager = SNPIndexManager(pvals,ss = stats_complete) #constructor allows for adding sum stats.
    snp_index_manager.addSubsets(args.split_scores.split(','))
    if args.plink_snp_matrix:
        snp_matrix = args.plink_snp_matrix
    else:
        snp_matrix = plinkToMatrix(geno_ids, args.plink_snps, local_pvar, args.plink_version, args.plink2_path)
    print("Parsing the genotype file...")
    if not args.debug_pickle:
        scores, patient_ids = calculatePRS(pvals, snp_matrix, snp_index_manager,stats_complete[REFID], num_pat, args)
        updateLog("Total number of SNPs for use in PRS calculations:", str(getEntryCount(geno_ids,False))) 
        if args.serialize:
            import pickle
            pickle.dump(scores, open("scores.dat", "wb"))
    else:
        import pickle
        scores = pickle.load(open(args.debug_pickle, 'rb'))
        patient_ids = ["empty"] * num_pat

    writeScores(scores, patient_ids, args)
    print("Scores written out to ", args.output + ".*.scores.tsv")
    stop = time.time()
    updateLog("Total runtime: "+ str(stop - start), True)
    #TODO: Add cleanup functionality, remove the massive files you've made.
