#!/usr/bin/env python
import argparse
import csv
import numpy as np
import os
import shlex
import shutil
import subprocess
import sys
import pandas as pd
import string
import glob
import multiprocessing as mp
import time

sys.path.append(os.getcwd())


def parse_arguments():
    parser = argparse.ArgumentParser(description='Part II: Calculating run in/on for sense and antisense strands.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--annotation_file_with_path', action= 'store', metavar='annotation_file_with_path', help='ensembl or wormbase annotation file for your organism.  The file may be downloaded from http://www.ensembl.org/info/data/ftp/index.html for ensembl or ftp://ftp.wormbase.org/pub/wormbase/releases/ for wormbase')
    parser.add_argument('--BedGraph_input_plu_strand', action= 'store', metavar='BedGraph_input_plu_strand', help='File that contains the plu strand BedGraph')
    parser.add_argument('--BedGraph_input_min_strand', action= 'store', metavar='BedGraph_input_min_strand', help='File that contains the min strand BedGraph')
    parser.add_argument('--output_prefix', action= 'store', metavar='output_prefix', help='Output prefix. Default: Dogcatcher_out/')
    parser.add_argument('--cpus', action= 'store', dest='cpus', metavar='cpus', default= 1, type=int, help='Enter available cpus per node.  The more cpus the faster Dogcatcher performs. Dogcatcher is designed to only work on one node. Default: "1"')
    parser.add_argument('--window_size', action= 'store', dest='window_size', metavar='window_size', default= 100, type=int, help='Sliding window size. Default: "100"')
    parser.add_argument('--coverage_percentage', action= 'store', dest='coverage_percentage', metavar='coverage_percentage', default= 90, type=float, help='Coverage percentage over sliding window. Default: "0.95"')
    parser.add_argument('--is_bedGraph', action= 'store', dest='is_bedGraph', metavar='is_bedGraph', default= False, help='Is the annotation file a bedGraph file.  This is also a compatible format.  The file needs to be a tab seperated bed with optional fields.  Ex. format name\tchr\tstart\tend. Default FALSE change to TRUE')
    parser.add_argument('--plu_DOG', action= 'store', dest='plu_DOG', metavar='plu_DOG', default= True, help='Default True change to False if there are multiprocessing issues and you want to process one condition at a time')
    parser.add_argument('--min_DOG', action= 'store', dest='min_DOG', metavar='min_DOG', default= True, help='Default True change to False if there are multiprocessing issues and you want to process one condition at a time')
    parser.add_argument('--plu_ADOG', action= 'store', dest='plu_ADOG', metavar='plu_ADOG', default= True, help='Default True change to False if there are multiprocessing issues and you want to process one condition at a time')
    parser.add_argument('--min_ADOG', action= 'store', dest='min_ADOG', metavar='min_ADOG', default= True, help='Default True change to False if there are multiprocessing issues and you want to process one condition at a time')

    parser.add_argument('--get_biotypes', action= 'store', dest='get_biotypes', metavar='get_biotypes', default= True, help='Default True change to False if you dont care about biotype overlap')

    args = parser.parse_args()
    return args

args = parse_arguments()

################################################ parameters
BedGraph_input_plu_strand = args.BedGraph_input_plu_strand
BedGraph_input_min_strand = args.BedGraph_input_min_strand
annotation_file_with_path = args.annotation_file_with_path
annotation_file = os.path.basename(annotation_file_with_path)
annotation_path = os.path.dirname(annotation_file_with_path) + "/"
get_biotypes = args.get_biotypes

cpus = args.cpus
is_bedGraph = args.is_bedGraph
window_size = args.window_size
coverage_percentage = float(args.coverage_percentage)*(0.01)
output_prefix = args.output_prefix

samplename = BedGraph_input_plu_strand[:-13]
samplename_nopath = os.path.basename(BedGraph_input_plu_strand[:-13])

print("Starting Part: II")
################################################################################

################################################################################
################################################################################
print("Checking for flattened gtf or BedGraph")
df_gtf_plu_NO_inside_genes = annotation_file_with_path[:-4] + "_plu_NO_inside_genes_same_strand_with_overlap.txt"
df_gtf_min_NO_inside_genes = annotation_file_with_path[:-4] + "_min_NO_inside_genes_same_strand_with_overlap.txt"
gtf_plu_exists = os.path.isfile(df_gtf_plu_NO_inside_genes)
gtf_min_exists = os.path.isfile(df_gtf_min_NO_inside_genes)

if is_bedGraph == True:
    print("Will implement later...")

#Check if the input is a gtf
if gtf_plu_exists == False:
     gtf_min_exists == False #Turn both of the checks to false...

if gtf_min_exists == False:
    print("Error loading gtf. Please run part I again and make sure all gtf files are in the correct location.")

if gtf_plu_exists and gtf_min_exists == True:
    print("Detected gtf files...")
    my_col = ["chr","start","end","gene_id"]

    df_gtf_all_genes = pd.read_csv(annotation_file_with_path[:-4] + "_ALL_GENES_NO_PATCHES.txt", sep="\t", dtype={"chr": str}, low_memory=False)

    df_gtf_plu_bed = pd.read_csv(annotation_file_with_path[:-4] + "_plu_NO_inside_genes_same_strand.BedGraph", sep="\t",header=None, comment="#", names = my_col, dtype={"chr" : str, "start" : int,"end" : int,"gene_id" : str}, low_memory=False)
    df_gtf_min_bed = pd.read_csv(annotation_file_with_path[:-4] + "_min_NO_inside_genes_same_strand.BedGraph", sep="\t",header=None, comment="#", names = my_col, dtype={"chr" : str, "start" : int,"end" : int,"gene_id" : str}, low_memory=False)

    df_gtf_plu_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_plu_inside_genes_same_strand.txt", sep="\t", dtype={"chr": str}, low_memory=False)
    df_gtf_min_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_min_inside_genes_same_strand.txt", sep="\t", dtype={"chr": str}, low_memory=False)

    df_gtf_plu_NO_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_plu_NO_inside_genes_same_strand_with_overlap.txt", sep="\t", dtype={"chr": str}, low_memory=False)
    df_gtf_min_NO_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_min_NO_inside_genes_same_strand_with_overlap.txt", sep="\t", dtype={"chr": str}, low_memory=False)

    df_gtf_plu_inside_genes_opposite_strand = pd.read_csv(annotation_file_with_path[:-4] + "_plu_inside_genes_opposite_strand.txt", sep="\t", dtype={"chr": str}, low_memory=False)
    df_gtf_min_inside_genes_opposite_strand = pd.read_csv(annotation_file_with_path[:-4] + "_min_inside_genes_opposite_strand.txt", sep="\t", dtype={"chr": str}, low_memory=False)

    df_gtf_plu_ALL_with_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_plu_ALL_with_inside_genes.txt", sep="\t", dtype={"chr": str}, low_memory=False)
    df_gtf_min_ALL_with_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_min_ALL_with_inside_genes.txt", sep="\t", dtype={"chr": str}, low_memory=False)

    #print("Detected biotypes are... ")
    chr_names = df_gtf_all_genes["chr"].unique().tolist()
    names_biotype_plu = df_gtf_plu_ALL_with_inside_genes["gene_biotype_name"].unique().tolist()
    names_biotype_min = df_gtf_min_ALL_with_inside_genes["gene_biotype_name"].unique().tolist()
    names_biotype = names_biotype_plu + names_biotype_min
    names_biotype = list(set(names_biotype))   #Change to set to get uniques then change back again
    print("Detected biotypes are ...")
    print(names_biotype)


#####Make bedgraphs for each chr#####################
def make_temp_bedgraph_parallel(chr):

        directory = samplename + "/Dogcatcher_temporary_BedGraph/"
        if not os.path.exists(directory):
            os.makedirs(directory)


        fout = (directory + os.path.basename(f1)[:-13] + "_" + str(chr) + "_" + str(bedgraph_plu_min) + '.BedGraph')
        print("Created Dogcatcher_temporary_BedGraph for chr: ", chr, "*"*5, "Individual worker is", os.getpid())
        if os.path.isfile(fout) == True:#Clear the file from multiple runs if it exists
            outfile = open(fout, 'w')
        with open(f1, "r", newline="\n") as infile, open(fout, "w", newline="\n") as outfile:
            for line in infile:
                tab1 = line.index("\t")
                line_chr = line[:tab1]
                if line_chr == chr:
                    outfile.write(line)

if os.path.exists(samplename + "/Dogcatcher_temporary_BedGraph") == True:
    print("Detected temporary BedGraph files... If process didn't complete or files are truncated delete Bedgraph_input_path/Dogcatcher_temporary_BedGraph folder to run again...")
if os.path.exists(samplename + "/Dogcatcher_temporary_BedGraph") == False:
    print("Detected no outputprefix/temp file. Writing out temporary BedGraphs for each chromosome into BedGraph input folder...")
    os.makedirs(samplename + "/Dogcatcher_temporary_BedGraph")
    x = 0
    for f1 in [BedGraph_input_plu_strand, BedGraph_input_min_strand]:
        if x == 0:
            bedgraph_plu_min = "plu"

        if x == 1:
            bedgraph_plu_min = "min"
        #Start multiprocess with temp bedgraph
        pool    = mp.Pool(processes=cpus)
        multiple_results = [pool.apply_async(make_temp_bedgraph_parallel, args=(chr,)) for chr in chr_names]
        pool.close()
        pool.join()  # block at this line until all processes are done
        print("All multiprocess workers done for *** ", f1)
        x = x + 1

################################################################################
################################################################################
print('Preparing for analysis using Dogcatcher... on file: ' + samplename)


###############################
#Functions to run Dogcatcher
def get_initial_input_files(BedGraph_input_plu_strand, bedgraph_plu_or_min):
    """This function will make a list of all of the bedgraph files in input folder"""
    filelist = [os.path.basename(x) for x in glob.glob(samplename + "/*" + bedgraph_plu_or_min + '.BedGraph')]
    filelist_with_path = glob.glob(samplename + "*" + bedgraph_plu_or_min + '.BedGraph')
    return filelist, filelist_with_path

def get_chr_input_files(BedGraph_input_plu_strand, bedgraph_plu_or_min):
    """This function will make a list of all of the chr bedgraph files in temp folder"""
    filelist = [os.path.basename(x) for x in glob.glob(samplename + "/Dogcatcher_temporary_BedGraph/*" + bedgraph_plu_or_min + '.BedGraph')]
    filelist_with_path = glob.glob(samplename + "/Dogcatcher_temporary_BedGraph/*" + bedgraph_plu_or_min + '.BedGraph')
    return filelist, filelist_with_path

def get_cov_perc(df, chr, start, end):
    """this function will take start and end and return the percent coverage over the area.
    for a BedGraph file for only one chromosome"""
    rows_df = df[ (df["end"] > start) & (df["start"] < end) ]
    if len(rows_df) == 0:
        return 0.0
    else:
        start_no_cov_length = 0
        if start > rows_df["start"].iloc[0]:
            start_no_cov_length = start - rows_df["start"].iloc[0]

        end_no_cov_length = 0
        if end < rows_df["end"].iloc[-1]:
            end_no_cov_length = rows_df["end"].iloc[-1] - end

        cov_sum = rows_df["length"].sum() - (start_no_cov_length + end_no_cov_length)    #add new start and end distance and all rows inbetween
        cov_length = end - start
        if cov_length == 0:
            cov_per = 0.0
            return cov_per
        else:
            cov_per = cov_sum / cov_length
            return cov_per

def header_int_to_char(df):
    """this function will change all numbers to letters for pandas"""
    import string
    x = 0  #rename int headers to letters
    for num in (list(df)):
        df.rename(columns={int(num) : str(string.ascii_uppercase[x])}, inplace=True)
        x = x + 1

def get_chrm_dic(df_gtf):
    """This function will take a gtf and return dictionary of different chrm. And clean the gtf for parsing"""
    df_gtf = df_gtf.copy()
    df_gtf["length"] = df_gtf["end"] - df_gtf["start"]
    #Make df of chrom
    names=df_gtf['chr'].unique().tolist()
    d_gtf_chr = {chrom : df_gtf.loc[df_gtf.chr==chrom] for chrom in names}
    return d_gtf_chr

def get_bg_df(fin):
    my_col = ["chr", "start", "end","count"]
    df_bed = pd.read_csv(fin, sep="\t", header=None, names=my_col, comment="#", dtype={"chr": str, "start": int, "end": int, "count": int}, low_memory=False)
    df_bed["length"] = df_bed["end"] - df_bed["start"]
    df_bed = df_bed[ df_bed["count"] > 1 ]
    return df_bed

def get_chrm_dic_w_biotype(df_gtf, chr_names, plu_or_min, gene_biotype):

    """This function will take a gtf and return strand specific dictionary of different chrm"""
    #Grab all genes that say protein coding
    names_biotype=df_gtf["gene_biotype_name"].unique().tolist()
    d_prot = {prot : df_gtf.loc[df_gtf.gene_biotype_name==prot] for prot in names_biotype}
    df_prot = d_prot[gene_biotype]  #Input for gene_biotype

    try:
        #Make positive and negative strand dictionary of dataframe
        names=df_prot['strand'].unique().tolist()
        d_strand = {strand : df_prot.loc[df_prot.strand==strand] for strand in names}
        df_plu_or_min = d_strand[plu_or_min]
        try:
            d_gtf_chr = {chrom : df_plu_or_min.loc[df_plu_or_min.chr==chrom] for chrom in chr_names}
            return d_gtf_chr
        except:
            d_gtf_chr = "NA"
            return d_gtf_chr
            print("No chr for this biotype")
    except:
            d_gtf_chr = "NA"
            return d_gtf_chr
            print("no dictionary for this strand")

def get_any_gene_overlap_to_list(df_gtf, chr, runon_start, runon_end):
        try:

            rows_df = df_gtf[ (runon_start < df_gtf["end"]) & (df_gtf["start"] < runon_end) ]
            if len(rows_df) == 0:
                return "NA"
            else:
                return rows_df["gene_id_name"].tolist()
        except:
            print("")

def checkForOverlap(df_gtf, runon_start, runon_end):
        try:
            rows_df = df_gtf[ (runon_start < df_gtf["end"]) & (df_gtf["start"] < runon_end) ]
            if len(rows_df) == 0:
                return False
            else:
                return True
        except:
            print("")


################################################################################
################################################################################
################################################################################
################################################################################
#RUNNING THE PROGRAM

print("Creating output folders...")
def make_output_folders_from_list(folder_list):
    for folder in folder_list:
        directory = output_prefix + folder
        if not os.path.exists(directory):
            os.makedirs(directory)

make_output_folders_from_list(["DOG/bed","DOG/gtf","DOG/csv","ADOG/bed","ADOG/gtf","ADOG/csv","DOG/Fstitch/"])

print("Generating chrm dictionaries of plu and minus strand gtf...")
d_gtf_chr = get_chrm_dic(df_gtf_all_genes)
d_gtf_plu_chr = get_chrm_dic(df_gtf_plu_NO_inside_genes)
d_gtf_min_chr = get_chrm_dic(df_gtf_min_NO_inside_genes)



################################################################################
def get_right_DOG(df, df_bed_chr, window_size, coverage_percentage):
        df5p = df
        df = df.copy()
        df = df.reset_index()
        df.fillna(value=0,inplace=True)    #fill with na
        df_left_overlap = df[ df["start"] < df["end"].shift(1) ]
        df_right_overlap = df[ df["start"].shift(-1) < df["end"] ].copy()
        df_both_end_overlap = df[ ( df["start"] < df["end"].shift(1) ) & ( df["start"].shift(-1) < df["end"] )]
        df["inside_gene"] = ( df["start"].shift(1) < df["start"] ) & ( df["end"] < df["end"].shift(1) )
        df["right_int_length"] = df["start"].shift(-1) - df["end"]       #Calculate 3p int after you get rid of inside genes
        df["right_overlap"] = df["start"].shift(-1) < df["end"]       #Get rid of genes that have another gene in 3' end.
        df = df[ ~df.right_overlap ]       #drop the values
        del df["inside_gene"]
        del df["right_overlap"]

        df["DOG_length"] = 0        #Initialize run on length column
        df_final = pd.DataFrame()      #Initialize final dataframe

        #For everything less than window size, give full length if over 95 cov. else 0
        try:
            df_LT_w = df [ df['right_int_length'] <= window_size ].copy()    #Less than window size
            #print("Number of genes under intergenic ", window_size, "window size = ", len(df_LT_w))
            df_LT_w['cov_perc'] = df_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + row["right_int_length"]), axis=1)
            df_LT_w["DOG_length"] = np.where(df_LT_w['cov_perc'] > coverage_percentage, df_LT_w['right_int_length'], 0 )
            df_final = df_final.append(df_LT_w)
        except:
            print("")

        #For genes over window size
        df.fillna(value=100000,inplace=True) #Give the last gene a 100kb dummy length (Or put chrom sizes in later)
        df = df [ df['right_int_length'] > window_size ]
        #print("Number of genes over intergenic ", window_size, "window size = ", len(df))
        counter = window_size
        while len(df) > 0:
            df['cov_perc'] = df.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + window_size), axis=1)
            df_under_95 = df [ df['cov_perc'] < coverage_percentage  ]
            if len(df_under_95) > 0:
                df_under_95 = df [ df['cov_perc'] < coverage_percentage  ].copy()
                df_under_95["DOG_length"] = df["DOG_length"] + window_size - counter
                df_final = df_final.append(df_under_95)
            df = df [ df['cov_perc'] > coverage_percentage  ]
            window_size = window_size + counter
        df_final = df_final.append(df_right_overlap)   #Add back genes in that had right overlap
        df = df_final                                 #Change name back to df_final for easier coding
        df.fillna(value=0,inplace=True)
        df = df.sort_values("start")



        df["TYPE"] = "DOG"
        df["DOG_start_local"] = df["end"]                                      #start coordinate of run on
        df["DOG_start_meta"] = df["end"]
        df["DOG_end_meta"] = df["DOG_start_local"] + df["DOG_length"]   #end coordinate of run on
        df["DOG_into_downstream_gene"] = ( df["DOG_length"] >= df['right_int_length'] ) & ( df["DOG_length"] > 0 )
        df["DOG_end_local"] = np.where(df["DOG_into_downstream_gene"] == True , df['start'].shift(-1), df["DOG_end_meta"] )
        df["DOG_from_upstream_gene"] = np.where(df["DOG_into_downstream_gene"].shift(1) == True , True, False)

        #df = df[ ~df.right_overlap ].copy()       #drop the values
        # #print("Number of genes without right overlap", len(df))
        # #Drop some of the columns not needed




        ###############################################################
        ###############################################################
        #####################Process 5p
        df5p = df5p.reset_index()

        df5p_left_overlap = df5p[df5p["start"] < df5p["end"].shift(1)]
        df5p['left_int_length'] = (df5p["end"].shift(1) - df5p["start"]) * (-1)  # Calculate 3p int after you get rid of
        df5p["left_overlap"] = df5p["start"] < df5p["end"].shift(1)  # Get rid of genes that have another gene in 3'

        df5p = df5p[ ~df5p.left_overlap ]       #drop the values
        del df5p["left_overlap"]
        df5p['DOG_length'] = 0  # Initialize run on length column
        df5p_final = pd.DataFrame()  # Initialize final dataframe
        try:
            df5p_LT_w = df5p [ df5p['left_int_length'] <= window_size ].copy()    #Less than window size
            #print("Number of genes under intergenic ", window_size, "window size = ", len(df5p_LT_w))
            df5p_LT_w['cov_perc'] = df5p_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - row['left_int_length'], row["start"]), axis=1)
            df5p_LT_w['DOG_length'] = np.where(df5p_LT_w['cov_perc'] > coverage_percentage , df5p_LT_w['left_int_length'], 0 )
            df5p_final = df5p_final.append(df5p_LT_w)
        except:
            print("")
        df5p.fillna(value=0, inplace=True)  # Give the first gene a 0 dummy length (Or put chrom sizes in later)
        df5p = df5p[df5p['left_int_length'] > window_size]
        counter = window_size
        while len(df5p) > 0:
            df5p['cov_perc'] = df5p.apply(
                lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - window_size, row['start']),
                axis=1)
            df5p_under_95 = df5p[df5p['cov_perc'] < coverage_percentage]
            if len(df5p_under_95) > 0:
                df5p_under_95 = df5p[df5p['cov_perc'] < coverage_percentage].copy()
                df5p_under_95['DOG_length'] = df5p['DOG_length'] + window_size - counter
                df5p_final = df5p_final.append(df5p_under_95)
            df5p = df5p[df5p['cov_perc'] > coverage_percentage]
            window_size = window_size + counter
        df5p_final = df5p_final.append(df5p_left_overlap).copy()  # Add back genes in that had right overlap
        df5p = df5p_final  # Change name back to df5p_final for easier coding
        df5p.fillna(value=0, inplace=True)
        df5p = df5p.sort_values("start")

        df5p["TYPE"] = "POG"
        df5p["DOG_end_local"] = df5p["start"]  # start coordinate of run on
        df5p["DOG_end_meta"] = df5p["start"]  # start coordinate of run on
        df5p["DOG_start_meta"] = df5p["DOG_end_meta"] - df5p['DOG_length']  # end coordinate of run on
        df5p["POG_into_upstream_gene"] = (df5p['DOG_length'] >= df5p['left_int_length']) & (df5p['DOG_length'] > 0)
        df5p["DOG_start_local"] = np.where(df5p["POG_into_upstream_gene"] == True, df5p['end'].shift(1),df5p["DOG_start_meta"])  # run on length when stops at next gene
        df5p["DOG_from_downstream_gene"] = np.where(df5p["POG_into_upstream_gene"].shift(1) == True, True,False)  # run on from upstream gene

        df5p.fillna(value=0, inplace=True)
        df = pd.concat([df,df5p])
        df = df[df["POG_into_upstream_gene"] != True ]
        return df

def running_plu_DOG(f1):
    file_name = os.path.basename(f1)
    print("Working on file...", file_name, "*"*5, "Individual worker is", os.getpid())
    df_bed_chr = get_bg_df(f1)
    df_final_all_chrm_all_files = pd.DataFrame()
    for chrom in df_bed_chr["chr"].unique().tolist(): #Do in a loop to get all chrm for plu genes
        try: #Use a try incase the gtf doesn't have that chr from the chr bed file
            df = get_right_DOG(d_gtf_plu_chr[chrom], df_bed_chr, window_size, coverage_percentage)
            df["sample_name"] = file_name[:-(14+len(str(chrom)))]
            df_final_all_chrm_all_files = df_final_all_chrm_all_files.append(df)
        except:
            df_final_all_chrm_all_files = df_final_all_chrm_all_files

    return df_final_all_chrm_all_files

def get_left_DOG(df, df_bed_chr, window_size, coverage_percentage):
    """This function will take a dataframe and return a bed file of left run on and gtf of all genes that have run on"""

    df = df.copy()
    df5p = df
    df = df.reset_index()
    # Make dataframes of genes that overlap, and are between other genes
    df_left_overlap = df[df["start"] < df["end"].shift(1)]
    df['left_int_length'] = (df["end"].shift(1) - df["start"]) * (-1)  # Calculate 3p int after you get rid of inside genes. turn to positive value
    df["left_overlap"] = df["start"] < df["end"].shift(1)  # Get rid of genes that have another gene in 3' end. (Keep while calculating 3p int length because 5' end is important)
    df = df[~df.left_overlap].copy()  # drop the values
    # Drop some of the columns not needed
    del df["left_overlap"]


    # Run on Loop (with coverage)
    df['DOG_length'] = 0  # Initialize run on length column
    df_final = pd.DataFrame()  # Initialize final dataframe

    # For everything less than window size, give full length if over 95 cov. else 0
    try:
        df_LT_w = df [ df['left_int_length'] <= window_size ].copy()    #Less than window size
        #print("Number of genes under intergenic ", window_size, "window size = ", len(df_LT_w))
        df_LT_w['cov_perc'] = df_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - row['left_int_length'], row["start"]), axis=1)
        df_LT_w['DOG_length'] = np.where(df_LT_w['cov_perc'] > coverage_percentage , df_LT_w['left_int_length'], 0 )
        df_final = df_final.append(df_LT_w)
    except:
        print("")
    # For genes over window size
    df.fillna(value=0, inplace=True)  # Give the first gene a 0 dummy length (Or put chrom sizes in later)
    df = df[df['left_int_length'] > window_size]
    #print("Number of genes over intergenic ", window_size, "window size = ", len(df))
    counter = window_size
    while len(df) > 0:
        df['cov_perc'] = df.apply(
            lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - window_size, row['start']),
            axis=1)
        df_under_95 = df[df['cov_perc'] < coverage_percentage]
        if len(df_under_95) > 0:
            df_under_95 = df[df['cov_perc'] < coverage_percentage].copy()
            df_under_95['DOG_length'] = df['DOG_length'] + window_size - counter
            df_final = df_final.append(df_under_95)
        df = df[df['cov_perc'] > coverage_percentage]
        window_size = window_size + counter

    df_final = df_final.append(df_left_overlap)  # Add back genes in that had right overlap
    df = df_final  # Change name back to df_final for easier coding
    df.fillna(value=0, inplace=True)


    df["TYPE"] = "DOG"

    df = df.sort_values("start")
    df["DOG_end_local"] = df["start"]  # start coordinate of run on
    df["DOG_end_meta"] = df["start"]  # start coordinate of run on
    df["DOG_start_meta"] = df["DOG_end_meta"] - df['DOG_length']  # end coordinate of run on
    df["DOG_into_downstream_gene"] = (df['DOG_length'] >= df['left_int_length']) & (df['DOG_length'] > 0)
    df["DOG_start_local"] = np.where(df["DOG_into_downstream_gene"] == True, df['end'].shift(1),df["DOG_start_meta"])
    df["DOG_from_upstream_gene"] = np.where(df["DOG_into_downstream_gene"].shift(1) == True, True,False)  # run on from upstream gene
    df.fillna(value=0, inplace=True)

    df_gene = df

    # ###########################################



    #############
    #Neeeeed to put in five prime

    df = df5p
    df = df.reset_index()
    df.fillna(value=0,inplace=True)    #fill with na
    df_left_overlap = df[ df["start"] < df["end"].shift(1) ]
    df_right_overlap = df[ df["start"].shift(-1) < df["end"] ].copy()
    df_both_end_overlap = df[ ( df["start"] < df["end"].shift(1) ) & ( df["start"].shift(-1) < df["end"] )]
    df["inside_gene"] = ( df["start"].shift(1) < df["start"] ) & ( df["end"] < df["end"].shift(1) )
    df["right_int_length"] = df["start"].shift(-1) - df["end"]       #Calculate 3p int after you get rid of inside genes
    df["right_overlap"] = df["start"].shift(-1) < df["end"]       #Get rid of genes that have another gene in 3' end.
    df = df[ ~df.right_overlap ]       #drop the values
    del df["right_overlap"]

    df["DOG_length"] = 0        #Initialize run on length column
    df_final = pd.DataFrame()      #Initialize final dataframe

    #For everything less than window size, give full length if over 95 cov. else 0
    try:
        df_LT_w = df [ df['right_int_length'] <= window_size ].copy()    #Less than window size
        #print("Number of genes under intergenic ", window_size, "window size = ", len(df_LT_w))
        df_LT_w['cov_perc'] = df_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + row["right_int_length"]), axis=1)
        df_LT_w["DOG_length"] = np.where(df_LT_w['cov_perc'] > coverage_percentage, df_LT_w['right_int_length'], 0 )
        df_final = df_final.append(df_LT_w)
    except:
        print("")

    #For genes over window size
    df.fillna(value=100000,inplace=True) #Give the last gene a 100kb dummy length (Or put chrom sizes in later)
    df = df [ df['right_int_length'] > window_size ]
    #print("Number of genes over intergenic ", window_size, "window size = ", len(df))
    counter = window_size
    while len(df) > 0:
        df['cov_perc'] = df.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + window_size), axis=1)
        df_under_95 = df [ df['cov_perc'] < coverage_percentage  ]
        if len(df_under_95) > 0:
            df_under_95 = df [ df['cov_perc'] < coverage_percentage  ].copy()
            df_under_95["DOG_length"] = df["DOG_length"] + window_size - counter
            df_final = df_final.append(df_under_95)
        df = df [ df['cov_perc'] > coverage_percentage  ]
        window_size = window_size + counter
    df_final = df_final.append(df_right_overlap)   #Add back genes in that had right overlap
    df = df_final                                  #Change name back to df_final for easier coding
    df.fillna(value=0,inplace=True)
    df = df.sort_values("start")

    df["TYPE"] = "POG"
    df["DOG_start_local"] = df["end"]
    df["DOG_start_meta"] = df["end"]#start coordinate of run on
    df["DOG_end_meta"] = df["DOG_start_meta"] + df["DOG_length"]   #end coordinate of run on
    df["POG_into_upstream_gene"] = ( df["DOG_length"] >= df['right_int_length'] ) & ( df["DOG_length"] > 0 )

    df["DOG_end_local"] = np.where(df["POG_into_upstream_gene"] == True , df['start'].shift(-1), df["DOG_end_meta"] )

    df["DOG_from_downstream_gene"] = np.where(df["POG_into_upstream_gene"].shift(1) == True, True,False)  # run on from upstream gene
    # upstream gene
    df.fillna(value=0, inplace=True)
    df = pd.concat([df_gene,df])
    df = df[df["POG_into_upstream_gene"] != True ]
    return df

def running_min_DOG(f1):
    file_name = os.path.basename(f1)
    print("Working on file...", file_name, "*"*5, "Individual worker is", os.getpid())
    df_bed_chr = get_bg_df(f1)
    df_final_all_chrm_all_files = pd.DataFrame()
    for chrom in df_bed_chr["chr"].unique().tolist(): #Do in a loop to get all chrm for plu genes
        try:
            df = get_left_DOG(d_gtf_min_chr[chrom], df_bed_chr, window_size, coverage_percentage)
            df["sample_name"] = file_name[:-(14+len(str(chrom)))]
            df_final_all_chrm_all_files = df_final_all_chrm_all_files.append(df)
        except:
            df_final_all_chrm_all_files = df_final_all_chrm_all_files
    return df_final_all_chrm_all_files


#runin
def get_right_ADOG(df, df_bed_chr, window_size, coverage_percentage, d_gtf_chr):
        df5p = df
        df = df.copy()
        df = df.reset_index()
        df.fillna(value=0,inplace=True)    #fill with na
        df_left_overlap = df[ df["start"] < df["end"].shift(1) ]
        df_right_overlap = df[ df["start"].shift(-1) < df["end"] ].copy()
        df_both_end_overlap = df[ ( df["start"] < df["end"].shift(1) ) & ( df["start"].shift(-1) < df["end"] )]
        df["inside_gene"] = ( df["start"].shift(1) < df["start"] ) & ( df["end"] < df["end"].shift(1) )
        df["right_int_length"] = df["start"].shift(-1) - df["end"]       #Calculate 3p int after you get rid of inside genes
        df["right_overlap"] = df["start"].shift(-1) < df["end"]       #Get rid of genes that have another gene in 3' end.
        df = df[ ~df.right_overlap ]       #drop the values
        del df["inside_gene"]
        del df["right_overlap"]

        df["DOG_length"] = 0        #Initialize run on length column
        df_final = pd.DataFrame()      #Initialize final dataframe

        #For everything less than window size, give full length if over 95 cov. else 0
        try:
            df_LT_w = df [ df['right_int_length'] <= window_size ].copy()    #Less than window size
            #print("Number of genes under intergenic ", window_size, "window size = ", len(df_LT_w))
            df_LT_w['cov_perc'] = df_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + row["right_int_length"]), axis=1)
            df_LT_w["DOG_length"] = np.where(df_LT_w['cov_perc'] > coverage_percentage, df_LT_w['right_int_length'], 0 )
            df_final = df_final.append(df_LT_w)
        except:
            print("")

        #For genes over window size
        df.fillna(value=100000,inplace=True) #Give the last gene a 100kb dummy length (Or put chrom sizes in later)
        df = df [ df['right_int_length'] > window_size ]
        #print("Number of genes over intergenic ", window_size, "window size = ", len(df))
        counter = window_size
        while len(df) > 0:
            df['cov_perc'] = df.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + window_size), axis=1)
            df_under_95 = df [ df['cov_perc'] < coverage_percentage  ]
            if len(df_under_95) > 0:
                df_under_95 = df [ df['cov_perc'] < coverage_percentage  ].copy()
                df_under_95["DOG_length"] = df["DOG_length"] + window_size - counter
                df_final = df_final.append(df_under_95)
            df = df [ df['cov_perc'] > coverage_percentage  ]
            window_size = window_size + counter
        df_final = df_final.append(df_right_overlap)   #Add back genes in that had right overlap
        df = df_final                                 #Change name back to df_final for easier coding
        df.fillna(value=0,inplace=True)
        df = df.sort_values("start")



        df["TYPE"] = "ADOG"
        df["DOG_start_local"] = df["end"]                                      #start coordinate of run on
        df["DOG_start_meta"] = df["end"]
        df["DOG_end_meta"] = df["DOG_start_local"] + df["DOG_length"]   #end coordinate of run on
        df["DOG_into_downstream_gene"] = ( df["DOG_length"] >= df['right_int_length'] ) & ( df["DOG_length"] > 0 )
        df["DOG_end_local"] = np.where(df["DOG_into_downstream_gene"] == True , df['start'].shift(-1), df["DOG_end_meta"] )
        df["DOG_from_upstream_gene"] = np.where(df["DOG_into_downstream_gene"].shift(1) == True , True, False)


        ###############################################################
        ###############################################################
        #####################Process 5p
        df5p = df5p.reset_index()

        df5p_left_overlap = df5p[df5p["start"] < df5p["end"].shift(1)]
        df5p['left_int_length'] = (df5p["end"].shift(1) - df5p["start"]) * (-1)  # Calculate 3p int after you get rid of
        df5p["left_overlap"] = df5p["start"] < df5p["end"].shift(1)  # Get rid of genes that have another gene in 3'

        df5p = df5p[ ~df5p.left_overlap ]       #drop the values
        del df5p["left_overlap"]
        
        df5p['DOG_length'] = 0  # Initialize run on length column
        df5p_final = pd.DataFrame()  # Initialize final dataframe

        try:
            df5p_LT_w = df5p [ df5p['left_int_length'] <= window_size ].copy()    #Less than window size
            #print("Number of genes under intergenic ", window_size, "window size = ", len(df5p_LT_w))
            df5p_LT_w['cov_perc'] = df5p_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - row['left_int_length'], row["start"]), axis=1)
            df5p_LT_w['DOG_length'] = np.where(df5p_LT_w['cov_perc'] > coverage_percentage , df5p_LT_w['left_int_length'], 0 )
            df5p_final = df5p_final.append(df5p_LT_w)
        except:
            print("")
        df5p.fillna(value=0, inplace=True)  # Give the first gene a 0 dummy length (Or put chrom sizes in later)
        df5p = df5p[df5p['left_int_length'] > window_size]
        counter = window_size
        while len(df5p) > 0:
            df5p['cov_perc'] = df5p.apply(
                lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - window_size, row['start']),
                axis=1)
            df5p_under_95 = df5p[df5p['cov_perc'] < coverage_percentage]
            if len(df5p_under_95) > 0:
                df5p_under_95 = df5p[df5p['cov_perc'] < coverage_percentage].copy()
                df5p_under_95['DOG_length'] = df5p['DOG_length'] + window_size - counter
                df5p_final = df5p_final.append(df5p_under_95)
            df5p = df5p[df5p['cov_perc'] > coverage_percentage]
            window_size = window_size + counter
        df5p_final = df5p_final.append(df5p_left_overlap)  # Add back genes in that had right overlap
        df5p = df5p_final  # Change name back to df5p_final for easier coding
        df5p.fillna(value=0, inplace=True)
        df5p = df5p.sort_values("start")



        df5p["TYPE"] = "APOG"
        df5p["DOG_end_local"] = df5p["start"]  # start coordinate of run on
        df5p["DOG_end_meta"] = df5p["start"]  # start coordinate of run on
        df5p["DOG_start_meta"] = df5p["DOG_end_meta"] - df5p['DOG_length']  # end coordinate of run on
        df5p["POG_into_upstream_gene"] = (df5p['DOG_length'] >= df5p['left_int_length']) & (df5p['DOG_length'] > 0)
        df5p["DOG_start_local"] = np.where(df5p["POG_into_upstream_gene"] == True, df5p['end'].shift(1),df5p["DOG_start_meta"])  # run on length when stops at next gene
        df5p["DOG_from_downstream_gene"] = np.where(df5p["POG_into_upstream_gene"].shift(1) == True, True,False)  # run on from upstream gene

        df5p.fillna(value=0, inplace=True)
        df = pd.concat([df,df5p])
        df = df[df["POG_into_upstream_gene"] != True ]

        df["AnyGeneOverlap"] = df.apply(lambda row: checkForOverlap(d_gtf_chr, row["DOG_start_meta"], row["DOG_end_meta"]), axis=1)
        df = df[ df["AnyGeneOverlap"] == False]

        return df

def running_plu_ADOG(f1):
    file_name = os.path.basename(f1)
    df_final_all_chrm_all_files = pd.DataFrame()

    print("Working on file...", file_name, "*"*5, "Individual worker is", os.getpid())
    df_bed_chr = get_bg_df(f1)
    df_final_all_chrm_all_files = pd.DataFrame()
    for chrom in df_bed_chr["chr"].unique().tolist(): #Do in a loop to get all chrm for plu genes
        try:
            df = get_right_ADOG(d_gtf_plu_chr[chrom], df_bed_chr, window_size, coverage_percentage,d_gtf_chr[chrom])
            df["sample_name"] = file_name[:-(14+len(str(chrom)))]
            df_final_all_chrm_all_files = df_final_all_chrm_all_files.append(df)
        except:
            df_final_all_chrm_all_files = df_final_all_chrm_all_files
    return df_final_all_chrm_all_files

def get_left_ADOG(df, df_bed_chr, window_size, coverage_percentage,d_gtf_chr):
        """This function will take a dataframe and return a bed file of left run on and gtf of all genes that have run on"""
        df = df.copy()
        df5p = df
        df = df.reset_index()
        # Make dataframes of genes that overlap, and are between other genes
        df_left_overlap = df[df["start"] < df["end"].shift(1)]
        df['left_int_length'] = (df["end"].shift(1) - df["start"]) * (-1)  # Calculate 3p int after you get rid of inside genes. turn to positive value
        df["left_overlap"] = df["start"] < df["end"].shift(1)  # Get rid of genes that have another gene in 3' end. (Keep while calculating 3p int length because 5' end is important)

        # Drop some of the columns not needed

        df = df[ ~df.left_overlap ]       #drop the values
        del df["left_overlap"]

        # Run on Loop (with coverage)
        df['DOG_length'] = 0  # Initialize run on length column
        df_final = pd.DataFrame()  # Initialize final dataframe

        # For everything less than window size, give full length if over 95 cov. else 0
        try:
            df_LT_w = df [ df['left_int_length'] <= window_size ].copy()    #Less than window size
            #print("Number of genes under intergenic ", window_size, "window size = ", len(df_LT_w))
            df_LT_w['cov_perc'] = df_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - row['left_int_length'], row["start"]), axis=1)
            df_LT_w['DOG_length'] = np.where(df_LT_w['cov_perc'] > coverage_percentage , df_LT_w['left_int_length'], 0 )
            df_final = df_final.append(df_LT_w)
        except:
            print("")
        # For genes over window size
        df.fillna(value=0, inplace=True)  # Give the first gene a 0 dummy length (Or put chrom sizes in later)
        df = df[df['left_int_length'] > window_size]
        #print("Number of genes over intergenic ", window_size, "window size = ", len(df))
        counter = window_size
        while len(df) > 0:
            df['cov_perc'] = df.apply(
                lambda row: get_cov_perc(df_bed_chr, row['chr'], row['start'] - window_size, row['start']),
                axis=1)
            df_under_95 = df[df['cov_perc'] < coverage_percentage]
            if len(df_under_95) > 0:
                df_under_95 = df[df['cov_perc'] < coverage_percentage].copy()
                df_under_95['DOG_length'] = df['DOG_length'] + window_size - counter
                df_final = df_final.append(df_under_95)
            df = df[df['cov_perc'] > coverage_percentage]
            window_size = window_size + counter

        df_final = df_final.append(df_left_overlap)   #Add back genes in that had right overlap
        df = df_final.copy()  # Change name back to df_final for easier coding
        df.fillna(value=0, inplace=True)


        df["TYPE"] = "ADOG"

        df = df.sort_values("start")
        df["DOG_end_local"] = df["start"]  # start coordinate of run on
        df["DOG_end_meta"] = df["start"]  # start coordinate of run on
        df["DOG_start_meta"] = df["DOG_end_meta"] - df['DOG_length']  # end coordinate of run on
        df["DOG_into_downstream_gene"] = (df['DOG_length'] >= df['left_int_length']) & (df['DOG_length'] > 0)
        df["DOG_start_local"] = np.where(df["DOG_into_downstream_gene"] == True, df['end'].shift(1),df["DOG_start_meta"])
        df["DOG_from_upstream_gene"] = np.where(df["DOG_into_downstream_gene"].shift(1) == True, True,False)  # run on from upstream gene
        df.fillna(value=0, inplace=True)

        df_gene = df

        # ###########################################



        #############
        #Neeeeed to put in five prime

        df = df5p
        df = df.reset_index()
        df.fillna(value=0,inplace=True)    #fill with na
        df_left_overlap = df[ df["start"] < df["end"].shift(1) ]
        df_right_overlap = df[ df["start"].shift(-1) < df["end"] ].copy()
        df_both_end_overlap = df[ ( df["start"] < df["end"].shift(1) ) & ( df["start"].shift(-1) < df["end"] )]
        df["inside_gene"] = ( df["start"].shift(1) < df["start"] ) & ( df["end"] < df["end"].shift(1) )
        df["right_int_length"] = df["start"].shift(-1) - df["end"]       #Calculate 3p int after you get rid of inside genes
        df["right_overlap"] = df["start"].shift(-1) < df["end"]       #Get rid of genes that have another gene in 3' end.
        df = df[ ~df.right_overlap ]       #drop the values
        del df["inside_gene"]
        del df["right_overlap"]


        df["DOG_length"] = 0        #Initialize run on length column
        df_final = pd.DataFrame()      #Initialize final dataframe

        #For everything less than window size, give full length if over 95 cov. else 0
        try:
            df_LT_w = df [ df['right_int_length'] <= window_size ].copy()    #Less than window size
            #print("Number of genes under intergenic ", window_size, "window size = ", len(df_LT_w))
            df_LT_w['cov_perc'] = df_LT_w.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + row["right_int_length"]), axis=1)
            df_LT_w["DOG_length"] = np.where(df_LT_w['cov_perc'] > coverage_percentage, df_LT_w['right_int_length'], 0 )
            df_final = df_final.append(df_LT_w)
        except:
            print("")

        #For genes over window size
        df.fillna(value=100000,inplace=True) #Give the last gene a 100kb dummy length (Or put chrom sizes in later)
        df = df [ df['right_int_length'] > window_size ]
        #print("Number of genes over intergenic ", window_size, "window size = ", len(df))
        counter = window_size
        while len(df) > 0:
            df['cov_perc'] = df.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + window_size), axis=1)
            df_under_95 = df [ df['cov_perc'] < coverage_percentage  ]
            if len(df_under_95) > 0:
                df_under_95 = df [ df['cov_perc'] < coverage_percentage  ].copy()
                df_under_95["DOG_length"] = df["DOG_length"] + window_size - counter
                df_final = df_final.append(df_under_95)
            df = df [ df['cov_perc'] > coverage_percentage  ]
            window_size = window_size + counter
        df_final = df_final.append(df_right_overlap)   #Add back genes in that had right overlap
        df = df_final                                  #Change name back to df_final for easier coding
        df.fillna(value=0,inplace=True)
        df = df.sort_values("start")

        df["TYPE"] = "APOG"
        df["DOG_start_local"] = df["end"]
        df["DOG_start_meta"] = df["end"]#start coordinate of run on
        df["DOG_end_meta"] = df["DOG_start_meta"] + df["DOG_length"]   #end coordinate of run on
        df["POG_into_upstream_gene"] = ( df["DOG_length"] >= df['right_int_length'] ) & ( df["DOG_length"] > 0 )

        df["DOG_end_local"] = np.where(df["POG_into_upstream_gene"] == True , df['start'].shift(-1), df["DOG_end_meta"] )

        df["DOG_from_downstream_gene"] = np.where(df["POG_into_upstream_gene"].shift(1) == True, True,False)  # run on from upstream gene

        # upstream gene
        df.fillna(value=0, inplace=True)
        df = pd.concat([df_gene,df])
        df = df[df["POG_into_upstream_gene"] != True ]



        df["AnyGeneOverlap"] = df.apply(lambda row: checkForOverlap(d_gtf_chr, row["DOG_start_meta"], row["DOG_end_meta"]), axis=1)
        df = df[ df["AnyGeneOverlap"] == False]


        return df

def running_min_ADOG(f1):
    file_name = os.path.basename(f1)
    print("Working on file...", file_name, "*"*5, "Individual worker is", os.getpid())
    df_bed_chr = get_bg_df(f1)
    df_final_all_chrm_all_files = pd.DataFrame()
    for chrom in df_bed_chr["chr"].unique().tolist(): #Do in a loop to get all chrm for plu genes
        try:
            df = get_left_ADOG(d_gtf_min_chr[chrom], df_bed_chr, window_size, coverage_percentage,d_gtf_chr[chrom])
            df["sample_name"] = file_name[:-(14+len(str(chrom)))]
            df_final_all_chrm_all_files = df_final_all_chrm_all_files.append(df)
        except:
            df_final_all_chrm_all_files = df_final_all_chrm_all_files
    return df_final_all_chrm_all_files


# Running in parallel
if __name__ == '__main__':
    args = parse_arguments()
    plu_DOG = args.plu_DOG
    min_DOG = args.min_DOG
    plu_ADOG = args.plu_ADOG
    min_ADOG = args.min_ADOG

    last_time = time.time()

    if plu_DOG == True:

        """This will run all of the files on the plus strand to calculate run on"""
        print("Processing plu strand DOG", "*-*"*40)
        gtf_plu_or_min = "plu"
        filelist, filelist_with_path = get_chr_input_files(BedGraph_input_plu_strand, gtf_plu_or_min)

        pool    = mp.Pool(processes=cpus)
        multiple_results = [pool.apply_async(running_plu_DOG, args=(f1,)) for f1 in filelist_with_path]
        pool.close()
        pool.join()  # block at this line until all processes are done
        print("All multiprocess workers done... ")
        df_all_files_list = ([ res.get() for res in multiple_results])
        #END PARALLEL PROCESSING

        df_all_files = pd.DataFrame()
        for df in df_all_files_list:
            df_all_files = df_all_files.append(df)

        #Clean up output
        df = df_all_files #Change to df for easier coding

        # print(df)
        df = df [ df['DOG_length'] > window_size]           #PLU Get only sections where the run on length is bigger than the base window
        df = df.copy()
        del df["cov_perc"]
        del df["index"]

        def get_plu_biotypes_parallel(chr):
            print("Working on biotypes for chr: ", chr, "*"*5, "Individual worker is", os.getpid())
            for biotype in names_biotype:
                #print("working on run on overlap for biotype:" + biotype)
                d_same_strand_gtf_chr = get_chrm_dic_w_biotype(df_gtf_all_genes, chr_names, "+", biotype)
                d_opposite_strand_gtf_chr = get_chrm_dic_w_biotype(df_gtf_all_genes, chr_names, "-", biotype)
                if d_same_strand_gtf_chr != "NA":
                    #print("working on runon_overlap_gene_same_strand...")
                    df['DOG_overlap_same_strand_local:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_same_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_meta"], row["DOG_end_local"]), axis=1)
                    df['DOG_overlap_same_strand_meta:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_same_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_meta"], row["DOG_end_meta"]), axis=1)
                if d_opposite_strand_gtf_chr != "NA":
                    #print("working on runon_overlap_gene_opposite_strand...")
                    df['DOG_overlap_opposite_strand_local:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_opposite_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_local"], row["DOG_end_local"]), axis=1)
                    df['DOG_overlap_opposite_strand_meta:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_opposite_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_meta"], row["DOG_end_meta"]), axis=1)
            return df


        if get_biotypes == True:
            #Start multiprocess with biotypes
            chr_names = df["chr"].unique().tolist()

            pool    = mp.Pool(processes=cpus)
            multiple_results = [pool.apply_async(get_plu_biotypes_parallel, args=(chr,)) for chr in chr_names]
            pool.close()
            pool.join()  # block at this line until all processes are done
            print("All multiprocess workers done... ")
            df_all_files_list = ([ res.get() for res in multiple_results])
            #END PARALLEL PROCESSING
            df_all_files = pd.DataFrame()
            for df in df_all_files_list:
                df_all_files = df_all_files.append(df)
            df = df_all_files
            df.fillna(value="0",inplace=True)



        print("Done with biotypes.........")

        df["DOG_start_local"] = df["DOG_start_local"].astype(int)
        df["DOG_start_meta"] = df["DOG_start_local"].astype(int)
        df["DOG_end_local"] = df["DOG_end_local"].astype(int)
        df["DOG_end_meta"] = df["DOG_end_meta"].astype(int)
        df["DOG_length"] = df["DOG_length"].astype(int)


        #Training samples for Fstitch
        df_2 = df[ ["chr", "DOG_start_local", "DOG_end_local"] ].copy()
        df_2["number"] = "1"
        df_2.to_csv(output_prefix + "DOG/Fstitch/" + samplename_nopath + "_"+ gtf_plu_or_min + "_gtf_FSTITCH_train_DOG.bed", sep="\t", index=None, quoting=csv.QUOTE_NONE, header=None)

        DOG_out_name_bed = output_prefix + "/DOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_DOG.bed"
        DOG_out_name_csv = output_prefix + "/DOG/csv/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_DOG.csv"
        DOG_out_name_gtf = output_prefix + "/DOG/gtf/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_DOG.gtf"

        #Bedgraph file

        df_bed = df[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df_bed.to_csv(DOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)

        df.to_csv(DOG_out_name_csv, sep="\t", index=None, quoting=csv.QUOTE_NONE)

        df_DOG_out_name_gtf = df[["chr","source","type","DOG_start_local", "DOG_end_local","dot","strand","dot2","gene_id"]]
        df_DOG_out_name_gtf.to_csv(DOG_out_name_gtf, sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)

        ################### 5p
        df5p = df[df["TYPE"]=="POG"]
        DOG_out_name_bed = output_prefix + "/DOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_POG.bed"
        df5p_bed = df5p[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df5p_bed.to_csv(DOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)


    if min_DOG == True:

        print("Processing min strand DOG", "*-*"*10)
        gtf_plu_or_min = "min"
        filelist, filelist_with_path = get_chr_input_files(BedGraph_input_plu_strand, gtf_plu_or_min)

        x = 0
        pool    = mp.Pool(processes=cpus)
        multiple_results = [pool.apply_async(running_min_DOG, args=(f1,)) for f1 in filelist_with_path]
        pool.close()
        pool.join()  # block at this line until all processes are done
        print("All multiprocess workers done... ")
        df_all_files_list = ([ res.get() for res in multiple_results])
        #END PARALLEL PROCESSING

        df_all_files = pd.DataFrame()
        for df in df_all_files_list:
            df_all_files = df_all_files.append(df)

        #Clean up output
        df = df_all_files #Change to df for easier coding


        ################################################################################

        # Clean up output
        df = df_all_files #Change to df for easier coding
        df.fillna(value=0,inplace=True)

        df = df [ df['DOG_length'] > window_size]
        df = df.copy()
        del df["cov_perc"]
        del df["index"]


        def get_min_biotypes_parallel(chr):
            print("Working on biotypes for chr: ", chr, "*"*5, "Individual worker is", os.getpid())
            for biotype in names_biotype:
                d_same_strand_gtf_chr = get_chrm_dic_w_biotype(df_gtf_all_genes, chr_names, "-", biotype)
                d_opposite_strand_gtf_chr = get_chrm_dic_w_biotype(df_gtf_all_genes, chr_names, "+", biotype)

                if d_same_strand_gtf_chr != "NA":
                    # print("working on runon_overlap_gene_same_strand...")
                    df['DOG_overlap_same_strand_local:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_same_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_local"], row["DOG_end_local"]), axis=1)
                    df['DOG_overlap_same_strand_meta:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_same_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_meta"], row["DOG_end_meta"]), axis=1)
                if d_opposite_strand_gtf_chr != "NA":
                    #print("working on runon_overlap_gene_opposite_strand...")
                    df['DOG_overlap_opposite_strand_local:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_opposite_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_local"], row["DOG_end_local"]), axis=1)
                    df['DOG_overlap_opposite_strand_meta:' + biotype] = df.apply(lambda row: get_any_gene_overlap_to_list(d_opposite_strand_gtf_chr[row['chr']], row["chr"], row["DOG_start_meta"], row["DOG_end_meta"]), axis=1)
            return df

        if get_biotypes == True:
            #Start multiprocess with biotypes
            chr_names = df["chr"].unique().tolist()

            pool    = mp.Pool(processes=cpus)
            multiple_results = [pool.apply_async(get_min_biotypes_parallel, args=(chr,)) for chr in chr_names]
            pool.close()
            pool.join()  # block at this line until all processes are done
            print("All multiprocess workers done... ")
            df_all_files_list = ([ res.get() for res in multiple_results])
            #END PARALLEL PROCESSING
            df_all_files = pd.DataFrame()
            for df in df_all_files_list:
                df_all_files = df_all_files.append(df)
            df = df_all_files


        df.fillna(value=0,inplace=True)



        df["DOG_start_local"] = df["DOG_start_local"].astype(int)
        df["DOG_start_meta"] = df["DOG_start_local"].astype(int)
        df["DOG_end_local"] = df["DOG_end_local"].astype(int)
        df["DOG_end_meta"] = df["DOG_end_meta"].astype(int)
        df["DOG_length"] = df["DOG_length"].astype(int)

        df5p = df[df["TYPE"] == "POG"]
        df = df[df["TYPE"] == "DOG"].copy()

        DOG_out_name_bed = output_prefix + "/DOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_DOG.bed"
        DOG_out_name_csv = output_prefix + "/DOG/csv/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_DOG.csv"
        DOG_out_name_gtf = output_prefix + "/DOG/gtf/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_DOG.gtf"

        df_bed_DOG = df[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df_bed_DOG.to_csv(DOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)

        #Write out total dataframe
        df.to_csv(DOG_out_name_csv, sep="\t", index=None, quoting=csv.QUOTE_NONE)

        df_gtf_DOG = df[ ["chr","source","type","DOG_start_local", "DOG_end_local","dot","strand","dot2","gene_id"]]
        df_gtf_DOG.to_csv(DOG_out_name_gtf, sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)


        #########5p
        DOG_out_name_bed = output_prefix + "/DOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_POG.bed"

        #df5p_all_DOG_longest = df5p.sort_values('DOG_length', ascending=False).drop_duplicates("gene_id_name")
        df5p_bed = df5p[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df5p_bed.to_csv(DOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)


    ####################################
    #RUNIN... ADOGS AND APOGS

    if plu_ADOG == True:
        """This will run all of the files on the plus gtf strand to calculate run in"""
        print("Processing plu strand ADOG", "*-*"*10)
        gtf_plu_or_min = "plu"
        
        filelist, filelist_with_path = get_chr_input_files(BedGraph_input_plu_strand, "min")
        pool    = mp.Pool(processes=cpus)

        multiple_results = [pool.apply_async(running_plu_ADOG, args=(f1,)) for f1 in filelist_with_path]
        pool.close()
        pool.join()  # block at this line until all processes are done
        print("All multiprocess workers done... ")
        df_all_files_list = ([ res.get() for res in multiple_results])
        #END PARALLEL PROCESSING

        df_all_files = pd.DataFrame()
        for df in df_all_files_list:
            df_all_files = df_all_files.append(df)

        #Clean up output
        df = df_all_files #Change to df for easier coding


        df = df [ df['DOG_length'] > window_size].copy()

        del df["cov_perc"]
        del df["index"]
        #Bedgraph file

        df["DOG_start_local"] = df["DOG_start_local"].astype(int)
        df["DOG_start_meta"] = df["DOG_start_local"].astype(int)
        df["DOG_end_local"] = df["DOG_end_local"].astype(int)
        df["DOG_end_meta"] = df["DOG_end_meta"].astype(int)
        df["DOG_length"] = df["DOG_length"].astype(int)



        ADOG_out_name_bed = output_prefix + "/ADOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_ADOG.bed"
        ADOG_out_name_csv = output_prefix + "/ADOG/csv/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_ADOG.csv"
        ADOG_out_name_gtf = output_prefix + "/ADOG/gtf/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_ADOG.gtf"

        df_bed = df[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df_bed.to_csv(ADOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)
        df.to_csv(ADOG_out_name_csv, sep="\t", index=None, quoting=csv.QUOTE_NONE)
        df_ADOG_out_name_gtf = df[["chr","source","type","DOG_start_local", "DOG_end_local","dot","strand","dot2","gene_id"]]
        df_ADOG_out_name_gtf.to_csv(ADOG_out_name_gtf, sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)
        print("Done with plu strand ADOG", "*-*"*20)


        ################### 5p

        df5p = df[df["TYPE"]=="APOG"]
        APOG_out_name_bed = output_prefix + "/ADOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_APOG.bed"

        df_bed = df5p[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df_bed.to_csv(APOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)

    if min_ADOG == True:
        """This will run all of the files on the min gtf strand to calculate run in"""
        print("Processing min strand ADOG", "*-*"*10)


        gtf_plu_or_min = "min"
        filelist, filelist_with_path = get_chr_input_files(BedGraph_input_plu_strand, "plu")
        pool    = mp.Pool(processes=cpus)

        multiple_results = [pool.apply_async(running_min_ADOG, args=(f1,)) for f1 in filelist_with_path]
        pool.close()
        pool.join()  # block at this line until all processes are done
        print("All multiprocess workers done... ")
        df_all_files_list = ([ res.get() for res in multiple_results])
        #END PARALLEL PROCESSING

        df_all_files = pd.DataFrame()
        for df in df_all_files_list:
            df_all_files = df_all_files.append(df)

        #Clean up output
        df = df_all_files #Change to df for easier coding
        df = df [ df['DOG_length'] > window_size].copy()

        del df["cov_perc"]
        del df["index"]
        #Bedgraph file

        df["DOG_start_local"] = df["DOG_start_local"].astype(int)
        df["DOG_start_meta"] = df["DOG_start_local"].astype(int)
        df["DOG_end_local"] = df["DOG_end_local"].astype(int)
        df["DOG_end_meta"] = df["DOG_end_meta"].astype(int)
        df["DOG_length"] = df["DOG_length"].astype(int)

        df5p = df[df["TYPE"]=="APOG"]

        ADOG_out_name_bed = output_prefix + "/ADOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_ADOG.bed"
        ADOG_out_name_csv = output_prefix + "/ADOG/csv/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_ADOG.csv"
        ADOG_out_name_gtf = output_prefix + "/ADOG/gtf/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_ADOG.gtf"

        df_bed = df[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df_bed.to_csv(ADOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)
        df.to_csv(ADOG_out_name_csv, sep="\t", index=None, quoting=csv.QUOTE_NONE)
        df_ADOG_out_name_gtf = df[["chr","source","type","DOG_start_local", "DOG_end_local","dot","strand","dot2","gene_id"]]
        df_ADOG_out_name_gtf.to_csv(ADOG_out_name_gtf, sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)


        ################### 5p


        APOG_out_name_bed = output_prefix + "/ADOG/bed/" + samplename_nopath + "_" + gtf_plu_or_min + "_gtf_APOG.bed"
        df_bed = df5p[ ["chr", "DOG_start_local", "DOG_end_local", "gene_id_name"] ]
        df_bed.to_csv(APOG_out_name_bed, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)
        print("done with min strand ADOG", "*-*"*20)


    print("*-*"*40)
    print('Finished Part II... Time took {} seconds'.format(time.time()-last_time))



