#!/usr/bin/env python
import argparse
import csv
import numpy as np
import os
import shlex
import shutil
import subprocess
import pandas as pd
import string
import glob



def make_output_folders_from_list(folder_list):
    for folder in folder_list:
        directory = output_prefix + folder
        if not os.path.exists(directory):
            os.makedirs(directory)


def get_chrm_dic(df_gtf):
    """This function will take a gtf and return dictionary of different chrm. And clean the gtf for parsing"""
    df_gtf = df_gtf.copy()
    names=df_gtf['chr'].unique().tolist()
    d_gtf_chr = {chrom : df_gtf.loc[df_gtf.chr==chrom] for chrom in names}
    return d_gtf_chr, names






# Running in parallel
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Part II: Calculating run in/on for sense and antisense strands.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--Dogcatcher_plu_strand_list', nargs="*", action= 'store', metavar='Dogcatcher_plu_strand_list', help='List of files for your control samples')
    parser.add_argument('--Dogcatcher_min_strand_list', nargs="*", action= 'store', metavar='Dogcatcher_min_strand_list', help='List of files for treated samples')
    parser.add_argument('--input_R_template_file', action= 'store', metavar='input_R_template_file', help='input_R_template_file. Example "R_template/R_template_new.R" ')
    parser.add_argument('--output_prefix', action= 'store', metavar='output_prefix', default = "Dogcatcher_out/", help='Output prefix. KEEP THIS THE SAME FROM PART II: Default: Dogcatcher_out/')
    parser.add_argument('--input_prefix', action= 'store', metavar='input_prefix', default = "Dogcatcher_out/", help='input prefix. MAKE THIS PART II OUTPUT_PREFIX : Default: Dogcatcher_out/')
    parser.add_argument('--filter', action= 'store', metavar='filter', default= "longest", help="Length to filter, longest or shortest")

    args = parser.parse_args()
    Dogcatcher_plu_strand_list = args.Dogcatcher_plu_strand_list
    Dogcatcher_min_strand_list = args.Dogcatcher_min_strand_list

    filter = args.filter

    output_prefix = args.output_prefix
    input_prefix = args.input_prefix


    print("Starting 2.5_Dogcatcher_filter.py with ", filter, "contig..." )
    def FilterFunction(dfDOG,dfADOG,filter,gtf_plu_or_minDOG,gtf_plu_or_minADOG):
        if filter == "longest":
            dfDOG = dfDOG.sort_values('DOG_length', ascending=False).drop_duplicates("gene_id_name")
            dfADOG = dfADOG.sort_values('DOG_length', ascending=False).drop_duplicates("gene_id_name")

        if filter == "shortest":
            dfDOG = dfDOG.sort_values('DOG_length', ascending=True).drop_duplicates("gene_id_name")
            dfADOG = dfADOG.sort_values('DOG_length', ascending=True).drop_duplicates("gene_id_name")



        dfDOG.fillna(value=0,inplace=True)
        dfADOG.fillna(value=0,inplace=True)
        DOGdic,DOGchr = get_chrm_dic(dfDOG)


        def checkForOverlap(chr, ADOGstart, ADOGend):

                try:
                    dfdog = DOGdic[chr]

                    df_left_overlap = dfDOG[(dfDOG['DOG_start_local'] < ADOGend) & ( ADOGend < dfDOG['DOG_end_local']) ]
                    df_right_overlap = dfDOG[(dfDOG['DOG_start_local'] < ADOGstart ) & ( ADOGstart < dfDOG['DOG_end_local']) ]
                    df_DOG_inside = dfDOG[( ADOGstart < dfDOG['DOG_start_local']) & ( dfDOG['DOG_end_local'] < ADOGend) ]
                    x = 0
                    x = x + len(df_left_overlap)
                    x = x + len(df_right_overlap)
                    x = x + len(df_DOG_inside)

                    if x == 0:
                        return False
                    else:

                        return True
                except:
                    return False


        dfADOG['ADOG_overlap_DOG'] = dfADOG.apply(lambda row: checkForOverlap(row["chr"], row['DOG_start_local'], row['DOG_end_local']), axis=1)
        dfADOG = dfADOG[dfADOG["ADOG_overlap_DOG"]==False]

        DOG_out_name_bed = output_prefix + "/DOG_contigs/bed/DOG_contigs_" + gtf_plu_or_minDOG + "_gtf_DOG.bed"
        DOG_out_name_csv = output_prefix + "/DOG_contigs/csv/DOG_contigs_" + gtf_plu_or_minDOG + "_gtf_DOG.csv"
        DOG_out_name_gtf = output_prefix + "/DOG_contigs/gtf/DOG_contigs_" + gtf_plu_or_minDOG + "_gtf_DOG.gtf"
        ADOG_out_name_bed = output_prefix + "/ADOG_contigs/bed/ADOG_contigs_" + gtf_plu_or_minADOG + "_gtf_ADOG.bed"
        ADOG_out_name_csv = output_prefix + "/ADOG_contigs/csv/ADOG_contigs_" + gtf_plu_or_minADOG + "_gtf_ADOG.csv"
        ADOG_out_name_gtf = output_prefix + "/ADOG_contigs/gtf/ADOG_contigs_" + gtf_plu_or_minADOG + "_gtf_ADOG.gtf"

        dfDOG_all_run_on_bed = dfDOG[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        dfDOG_all_run_on_bed.to_csv(DOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)
        dfDOG.to_csv(DOG_out_name_csv, sep="\t", index=None, quoting=csv.QUOTE_NONE)
        dfDOG = dfDOG[["chr","source","type", "DOG_start_local", "DOG_end_local","dot","strand","dot2","gene_id"]]
        dfDOG.to_csv(DOG_out_name_gtf, sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)


        dfADOG_all_run_on_bed = dfADOG[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]

        dfADOG_all_run_on_bed.to_csv(ADOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)
        dfADOG.to_csv(ADOG_out_name_csv, sep="\t", index=None, quoting=csv.QUOTE_NONE)
        dfADOG = dfADOG[["chr","source","type", "DOG_start_local", "DOG_end_local","dot","strand","dot2","gene_id"]]
        dfADOG.to_csv(ADOG_out_name_gtf, sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)


    print("Finding contigs and putting into : " + output_prefix)
    make_output_folders_from_list(["DOG_contigs/bed","DOG_contigs/gtf","DOG_contigs/csv","ADOG_contigs/bed","ADOG_contigs/gtf","ADOG_contigs/csv"])


    df_all_files = pd.DataFrame()
    for file in Dogcatcher_plu_strand_list:
        samplename = file[:-13]
        df = pd.read_csv(input_prefix + "DOG/csv/" + file[:-13] + "_plu_gtf_DOG.csv", sep='\t')
        df_all_files = df_all_files.append(df)

    pluDOG = df_all_files

    df_all_files = pd.DataFrame()
    for file in Dogcatcher_min_strand_list:
        samplename = file[:-13]
        df = pd.read_csv(input_prefix + "DOG/csv/" + file[:-13] + "_min_gtf_DOG.csv", sep='\t')
        df_all_files = df_all_files.append(df)

    minDOG = df_all_files

    df_all_files = pd.DataFrame()
    for file in Dogcatcher_plu_strand_list:
        samplename = file[:-13]
        df = pd.read_csv(input_prefix + "ADOG/csv/" + file[:-13] + "_plu_gtf_ADOG.csv", sep='\t')
        df_all_files = df_all_files.append(df)
    pluADOG = df_all_files

    df_all_files = pd.DataFrame()
    for file in Dogcatcher_min_strand_list:
        samplename = file[:-13]
        df = pd.read_csv(input_prefix + "ADOG/csv/" + file[:-13] + "_min_gtf_ADOG.csv", sep='\t')
        df_all_files = df_all_files.append(df)
    minADOG = df_all_files


    FilterFunction(pluDOG,minADOG,filter, gtf_plu_or_minDOG="plu",gtf_plu_or_minADOG="min")
    FilterFunction(minDOG,pluADOG,filter, gtf_plu_or_minDOG="min",gtf_plu_or_minADOG="plu")



