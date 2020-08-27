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


def filterSumStatSNPs(ss, ss_type, keep_ambig_filter):
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
    #Set this to 0.3 from 0.4 on MArch 16
    ambig_filter = 0.3
    if keep_ambig_filter:
        ambig_filter = 100
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
        #NEAL has unusual case where occassionally ALT != minor allele. We check for that here:
        #1/27- this doesn't matter at all. The beta is in terms of ALT all the time. So don't do anything.
        #ss_t[REF] = np.where(ss_t["minor_allele"] == ss_t["f"],ss_t["s"],ss_t["f"])
        #ss_t[ALT] = np.where(ss_t["minor_allele"] == ss_t["s"],ss_t["s"],ss_t["f"])  
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
            print(str(sizes[t]) + " variants submitted in summary statistics file.")
        else:
            print(comment[t], str(sizes[prev]-sizes[t]))
            prev = t
    #Write out to file #ID, ref, alt, beta, pvalue
    #Check to see if any changes even made
    if ss_t.shape[0] == sizes["start"]:
        print("Summary stats alreadya appear to be filtered. No changes were made")
    ss_t.to_csv(out_name, sep = '\t', index = False) 
    return out_name


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Quick script to filter a summary stats file")
    parser.add_argument("-ss", "--sum_stats", help = "Path to summary stats file. Be sure to specify format if its not DEFAULT. Assumes tab delimited", required = True)
    parser.add_argument("--ss_format", help = "Format of the summary statistics", default = "SAIGE", choices = ["DEFAULT", "SAIGE", "NEAL", "HELIX"])    
    parser.add_argument("--maf", default = 0.01, type = float, help = "Specify what MAF cutoff to use.")
    parser.add_argument("-o", "--output", default = "./prs_" + str(date.today()), help = "Specify where you would like the output file to be written and its prefix.")
    parser.add_argument("--reset_ids", action = "store_true", help = "By default, program assumes genotype ids are in the form chr:loc:ref:alt. Specify this argument if this is NOT the case to.")
    parser.add_argument("--no_ambiguous_snps", help = "Specify this if you wish to remove all ambiguous SNPs whatsover, regardless of MAF. Default is to keep those with MAF < 0.3", action = "store_true", default = False)
    args = parser.parse_args()
    if args.reset_ids:
        VAR_IN_ID = False #Variant information is not included in genotype id
        print("Filtering ambiguous SNPs and indels from SS data...")
    args.sum_stats = filterSumStatSNPs(args.sum_stats, args.ss_format, args.no_ambiguous_snps)
    args.ss_format = "DEFAULT" #ID REF ALT BETA PVAL
    print("Ambiguous SNPs and indels removed. Filtered summary statistics written out to", args.sum_stats, True)
