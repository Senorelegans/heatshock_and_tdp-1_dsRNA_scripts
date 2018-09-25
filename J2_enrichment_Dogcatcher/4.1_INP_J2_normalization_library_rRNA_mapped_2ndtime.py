#!python3.5
# This script will normalize INP and J2 for Deseq2

#load pandas and numpy and regrex
import pandas as pd
import numpy as np
import os
import sys
import re
import argparse
import csv
import numpy as np
import os
import shlex
import shutil
import subprocess
import sys
import string
import glob

#this is finished avg treatment input
#f1 = sys.argv[1]

# Running in parallel
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Part 3.3: Normalizing Calculating run in/on for sense and antisense strands.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--control_TOTAL_BAM_list', nargs="*", action= 'store', metavar='control_TOTAL_BAM_list', help='List of files for your control samples')
    parser.add_argument('--treatment_TOTAL_BAM_list', nargs="*", action= 'store', metavar='treatment_TOTAL_BAM_list', help='List of files for treated samples')
    parser.add_argument('--control_J2_BAM_list', nargs="*", action= 'store', metavar='control_J2_BAM_list', help='List of files for your control samples')
    parser.add_argument('--treatment_J2_BAM_list', nargs="*", action= 'store', metavar='treatment_J2_BAM_list', help='List of files for treated samples')
    parser.add_argument('--input_prefix', action= 'store', metavar='input_prefix', default = "Dogcatcher_out/", help='input prefix. MAKE THIS PART II OUTPUT_PREFIX : Default: Dogcatcher_out/')

    parser.add_argument('--control_TOTAL_STAROUT_list', nargs="*", action= 'store', metavar='control_TOTAL_STAROUT_list', help='List of files for your control samples')
    parser.add_argument('--treatment_TOTAL_STAROUT_list', nargs="*", action= 'store', metavar='treatment_TOTAL_STAROUT_list', help='List of files for treated samples')
    parser.add_argument('--control_J2_STAROUT_list', nargs="*", action= 'store', metavar='control_J2_STAROUT_list', help='List of files for your control samples')
    parser.add_argument('--treatment_J2_STAROUT_list', nargs="*", action= 'store', metavar='treatment_J2_STAROUT_list', help='List of files for treated samples')

    args = parser.parse_args()

    control_TOTAL_BAM_list = args.control_TOTAL_BAM_list
    treatment_TOTAL_BAM_list = args.treatment_TOTAL_BAM_list
    control_J2_BAM_list = args.control_J2_BAM_list
    treatment_J2_BAM_list = args.treatment_J2_BAM_list
    
    control_TOTAL_STAROUT_list = args.control_TOTAL_STAROUT_list
    treatment_TOTAL_STAROUT_list = args.treatment_TOTAL_STAROUT_list
    control_J2_STAROUT_list = args.control_J2_STAROUT_list
    treatment_J2_STAROUT_list = args.treatment_J2_STAROUT_list



    input_prefix = args.input_prefix


    ###Generate col data file
    control_TOTAL_BAM_list = [ "X." + bam.replace("/",".").replace("-",".").replace("..",".") for bam in control_TOTAL_BAM_list ]
    treatment_TOTAL_BAM_list = [ "X." + bam.replace("/",".").replace("-",".").replace("..",".") for bam in treatment_TOTAL_BAM_list ]
    control_J2_BAM_list = [ "X." + bam.replace("/",".").replace("-",".").replace("..",".") for bam in control_J2_BAM_list ]
    treatment_J2_BAM_list = [ "X." + bam.replace("/",".").replace("-",".").replace("..",".") for bam in treatment_J2_BAM_list ]

    #output_prefix + "/" + "Rsubread_initial.R", "w"



    TOTAL_LIST = control_TOTAL_BAM_list + treatment_TOTAL_BAM_list
    J2_LIST=control_J2_BAM_list + treatment_J2_BAM_list



    ###Functions for normalization
    def clean_add_one_to_df(f1):
        """This function will delete the length column, add one, and return the df, and header (list of columns)"""
        df=pd.read_csv(f1, index_col="GeneID", sep='\t')
        del df["Length"]            # get rid of length column
        df.columns = control_TOTAL_BAM_list + treatment_TOTAL_BAM_list + control_J2_BAM_list + treatment_J2_BAM_list
        for col in df.columns:
            df[col] = df[col] + 1    #Add one to every point in df. So we don't divide by zero later
        return (df)

    def get_star_output(f1):
        """This function will get #Number of input reads - (Uniquely mapped reads number + Number of reads mapped to multiple loci)
        for each file from a star output. good for rRNA normalization"""
        with open(f1,"r") as infile:
            input_reads_list = []
            unique_reads_list = []
            multi_reads_list = []
            for line in infile:
                if "Number of input reads" in line:
                    scl = [m.start() for m in re.finditer(r"\t",line)]  #semicolon count list
                    input_reads = line[ scl[0]:].strip("\t")
                    input_reads_list.append(int(input_reads))
                if "Uniquely mapped reads number" in line:
                    scl = [m.start() for m in re.finditer(r"\t",line)]  #semicolon count list
                    unique_reads = line[ scl[0]:].strip("\t")
                    unique_reads_list.append(int(unique_reads))
                if "Number of reads mapped to multiple loci" in line:
                    scl = [m.start() for m in re.finditer(r"\t",line)]  #semicolon count list
                    multi_reads = line[ scl[0]:].strip("\t")
                    multi_reads_list.append(int(multi_reads))
            total_mapped = unique_reads_list[0] + multi_reads_list[0]
            total_no_mapped = input_reads_list[0] - total_mapped


            return total_no_mapped
    # get_star_output(control_TOTAL_STAROUT_list[0])



    #1. Collect mapping reads from star out
    input_total=0
    x = 0
    input_nomapped_list = []
    for file in control_TOTAL_STAROUT_list + treatment_TOTAL_STAROUT_list:
        file_total = get_star_output(file)
        input_nomapped_list.append(file_total)
        input_total = input_total + file_total
        x = x + 1
    inputavg = input_total / x
    x = 0
    j2_total = 0
    j2_nomapped_list = []
    for file in control_J2_STAROUT_list + treatment_J2_STAROUT_list:
        file_total = get_star_output(file)
        j2_nomapped_list.append(file_total)
        j2_total = j2_total + file_total
        x = x + 1
    j2avg = j2_total / x


    input_nomapped_ratio_list = []
    for x in input_nomapped_list:
        ratio = x/inputavg
        input_nomapped_ratio_list.append(ratio)

    j2_nomapped_ratio_list = []
    for x in j2_nomapped_list:
        ratio = x/j2avg
        j2_nomapped_ratio_list.append(ratio)



    ##############

    filelist = [input_prefix + "/min_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_Rsubread_sense.txt", \
    input_prefix + "/min_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_sense.txt", \
    input_prefix + "/plu_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_Rsubread_sense.txt", \
    input_prefix + "/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_sense.txt", \
    input_prefix + "/min_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_Rsubread_antisense.txt", \
    input_prefix + "/min_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_antisense.txt", \
    input_prefix + "/plu_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_Rsubread_antisense.txt", \
    input_prefix + "/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_antisense.txt", ]


    
    #2. Apply normalization plu_sense_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_Rsubread_senseunique.txt"
    for f1 in filelist:
        fN = f1[:-4]
        f1out = fN + "_normalized_rRNA_star.txt"

        df = clean_add_one_to_df(f1)
        x = 0
        for col in TOTAL_LIST:
            df[col] = df[col]/input_nomapped_ratio_list[x]
            df[col] = df[col].apply(round)
            x=x+1

        x=0
        for col in J2_LIST:
            df[col] = df[col]/j2_nomapped_ratio_list[x]
            df[col] = df[col].apply(round)
            x=x+1
        df.to_csv(f1out, sep="\t")
