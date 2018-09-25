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




# Running in parallel
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Part III: Creating gtf of ADOG (Anti-Sense) and DOG (Sense) contigs with Non-Significant genes.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--annotation_file_with_path', action= 'store', metavar='annotation_file_with_path', help='ensembl or wormbase annotation file for your organism.  The file may be downloaded from http://www.ensembl.org/info/data/ftp/index.html for ensembl or ftp://ftp.wormbase.org/pub/wormbase/releases/ for wormbase')
    parser.add_argument('--input_DESeq2_sense_file', action= 'store', metavar='input_DESeq2_sense_file', help='Raw Sense output file from DSeq2 (Must have no column header for first column). Example 0.0_input_files/')
    parser.add_argument('--input_DESeq2_antisense_file', action= 'store', metavar='input_DESeq2_antisense_file', help='Raw Anti-Sense output file from DSeq2 (Must have no column header for first column). Example 0.0_input_files/')
    parser.add_argument('--output_prefix', action= 'store', metavar='output_prefix', default = "Dogcatcher_out/", help='Output prefix. KEEP THIS THE SAME FROM PART II: Default: Dogcatcher_out/')
    parser.add_argument('--input_prefix', action= 'store', metavar='input_prefix', default = "Dogcatcher_out/", help='input prefix. MAKE THIS PART II OUTPUT_PREFIX : Default: Dogcatcher_out/')
    parser.add_argument('--input_prefix_initial', action= 'store', metavar='input_prefix_inital', default = "Dogcatcher_out/", help='input prefix. MAKE THIS PART II OUTPUT_PREFIX : Default: Dogcatcher_out/')
    parser.add_argument('--padj', action= 'store', dest='padj', metavar='padj', default= 0.05, type=float, help='DSeq2 padj value cutoff to select non-significant genes. Default: "0.05"')

    args = parser.parse_args()

    ################################################ parameters

    annotation_file_with_path = args.annotation_file_with_path
    annotation_file = os.path.basename(annotation_file_with_path)
    annotation_path = os.path.dirname(annotation_file_with_path)



    padj = args.padj
    output_prefix = args.output_prefix + "/"
    input_prefix = args.input_prefix + "/"

    input_prefix_initial = args.input_prefix_initial + "/"
    input_DESeq2_sense_file =     input_prefix_initial + "DESeq2_sense.csv"
    input_DESeq2_antisense_file = input_prefix_initial + "DESeq2_antisense.csv"



    #Input files from last script
    gtf_col_names = ["chr","source","type","start","end","dot","strand","dot2","gene_id"]
    #df_gtf_all = pd.read_csv(annotation_file_with_path[:-4] + "_ALL_GENES_NO_PATCHES.txt", sep="\t")

    df_gtf_plu_no_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_plu_NO_inside_genes_same_strand.gtf", names=gtf_col_names, sep="\t")
    df_gtf_min_no_inside_genes = pd.read_csv(annotation_file_with_path[:-4] + "_min_NO_inside_genes_same_strand.gtf", names=gtf_col_names, sep="\t")


    df_plu_DOG_csv = pd.read_csv(input_prefix + "DOG_contigs/csv/DOG_contigs_plu_gtf_DOG.csv", sep="\t")
    df_min_DOG_csv = pd.read_csv(input_prefix + "DOG_contigs/csv/DOG_contigs_min_gtf_DOG.csv", sep="\t")
    df_plu_ADOG_csv = pd.read_csv(input_prefix + "ADOG_contigs/csv/ADOG_contigs_plu_gtf_ADOG.csv", sep="\t")
    df_min_ADOG_csv = pd.read_csv(input_prefix + "ADOG_contigs/csv/ADOG_contigs_min_gtf_ADOG.csv", sep="\t")


    df_plu_DOG_gtf = pd.read_csv(input_prefix + "DOG_contigs/gtf/DOG_contigs_plu_gtf_DOG.gtf",    names=gtf_col_names, sep="\t")
    df_min_DOG_gtf = pd.read_csv(input_prefix + "DOG_contigs/gtf/DOG_contigs_min_gtf_DOG.gtf",    names=gtf_col_names, sep="\t")
    df_plu_ADOG_gtf = pd.read_csv(input_prefix + "ADOG_contigs/gtf/ADOG_contigs_plu_gtf_ADOG.gtf", names=gtf_col_names, sep="\t")
    df_min_ADOG_gtf = pd.read_csv(input_prefix + "ADOG_contigs/gtf/ADOG_contigs_min_gtf_ADOG.gtf", names=gtf_col_names, sep="\t")


    print("Starting Part III: DESeq2 and Rsubread...")
    print("Creating output folders...")
    directory = output_prefix
    if not os.path.exists(directory):
        os.makedirs(directory)

    ##############################################
    ##############################################
    ### Create gtf of run on/in contigs with      non-significant genes (without run on/in genes)
    #load original gtf

    def header_int_to_char(df):
        """this function will change all numbers to letters for pandas"""
        import string
        x = 0  #rename int headers to letters
        for num in (list(df)):
            df.rename(columns={int(num) : str(string.ascii_uppercase[x])}, inplace=True)
            x = x + 1

    def get_gene_id_name_column(df_gtf):
        """This function will take a gtf and unpack the gene id and biotypes column"""
        #Expand the gene id column, then merge back to get gene_id name in new column
        df_gtf["type"] = "exon"
        df_gene_id = df_gtf["gene_id"].str.split(';',expand=True)
        header_int_to_char(df_gene_id)
        for column_name in list(df_gene_id): #Find column that has gene_biotype in it
            if df_gene_id[column_name].str.contains("gene_biotype").any() == True:
                biotype_header = column_name
        df_merge = pd.merge(
            df_gtf,
            df_gene_id,
            how="left",
            left_index = True,
            right_index = True,
        )
        #print(df_merge)
        p = r'"([A-Za-z0-9_\./\\-]*)"'  #Extract gene name. everything in quotes
        df_merge["gene_id_name"] = (df_merge["A"]
                            .str
                            .extract(p, expand=True))
        df_merge = df_merge[["chr","source","type","start","end","dot","strand","dot2","gene_id", "gene_id_name"]]
        return df_merge

    def get_nonsig_dseq2(dseq_input, df_gtf_all, padj):
        """This function will take in a DSeq2 normalized matrix and make a gtf of all non-sig genes"""
        df = pd.read_csv(dseq_input, index_col=None, sep=',')
        df["not_significant"] = np.where(df["padj"] > padj, True, False)
        df_greater_padj = df[df.not_significant]  #Keep all True values. (Not significant)
        del df_greater_padj["not_significant"]
        df_na = df[pd.isnull(df['padj'])]
        del df_na["not_significant"]

        if len(df_na) > 0:
            df_dseq = pd.concat([df_greater_padj, df_na])
        else:
            df_dseq = df_greater_padj
        df_dseq = df_dseq.copy()
        df_dseq.rename(columns={"Unnamed: 0" : "gene_id_name"}, inplace=True)
        df_nonsig_gtf = df_gtf_all[ df_gtf_all.gene_id_name.isin(df_dseq.gene_id_name) ]
        df_nonsig_gtf = df_nonsig_gtf[["chr","source","type","start","end","dot","strand","dot2","gene_id", "gene_id_name"]]
        return df_nonsig_gtf

    def runion_gtf_with_nonsig_gtf(df_nonsig_gtf, df_runion_gtf):
        """Remove genes that have run-in or run-on from non-sig gtf.
        Then concatenate."""
        #Take out runion genic regions from nonsig gtf
        df_dseq_no_runion = df_nonsig_gtf [ ~df_nonsig_gtf.isin(df_runion_gtf.gene_id_name) ]
        #add run on/in intergenic regions to non-sig dseq2 gtf
        df_nonsig_gtf = pd.concat([df_dseq_no_runion, df_runion_gtf])
        del df_nonsig_gtf["gene_id_name"]
        return df_nonsig_gtf

    ############################################################
    #FOR SENSE RUN ON
    #get gene_id_name_column for gtf's
    df_plu_DOG_gtf = get_gene_id_name_column(df_plu_DOG_gtf)
    df_min_DOG_gtf = get_gene_id_name_column(df_min_DOG_gtf)

    df_gtf_plu_no_inside_genes = get_gene_id_name_column(df_gtf_plu_no_inside_genes)
    df_gtf_min_no_inside_genes = get_gene_id_name_column(df_gtf_min_no_inside_genes)

    #Get nonsig genes from deseq2
    df_plu_sense_nonsig_gtf = get_nonsig_dseq2(input_DESeq2_sense_file, df_gtf_plu_no_inside_genes, padj)
    df_min_sense_nonsig_gtf = get_nonsig_dseq2(input_DESeq2_sense_file, df_gtf_min_no_inside_genes, padj)

    #Combine non-sig gtf with runon
    df_plu_sense_nonsig_runon = runion_gtf_with_nonsig_gtf(df_plu_sense_nonsig_gtf, df_plu_DOG_gtf)
    df_min_sense_nonsig_runon = runion_gtf_with_nonsig_gtf(df_min_sense_nonsig_gtf, df_min_DOG_gtf)

    #print(df_plu_sense_nonsig_runon)
    df_plu_sense_nonsig_runon.to_csv(output_prefix + "plu_ALL_SAMPLES_DOG_with_nonsig.gtf", sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)
    df_min_sense_nonsig_runon.to_csv(output_prefix + "min_ALL_SAMPLES_DOG_with_nonsig.gtf", sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)


    ############################################################
    #FOR ANTISENSE RUN IN
    #Get nonsig genes from deseq2
    #get gene_id_name_column for gtf's

    df_plu_ADOG_gtf = get_gene_id_name_column(df_plu_ADOG_gtf)
    df_min_ADOG_gtf = get_gene_id_name_column(df_min_ADOG_gtf)

    df_gtf_plu_no_inside_genes = get_gene_id_name_column(df_gtf_plu_no_inside_genes)
    df_gtf_min_no_inside_genes = get_gene_id_name_column(df_gtf_min_no_inside_genes)

    df_plu_antisense_nonsig_gtf = get_nonsig_dseq2(input_DESeq2_antisense_file, df_gtf_plu_no_inside_genes, padj)
    df_min_antisense_nonsig_gtf = get_nonsig_dseq2(input_DESeq2_antisense_file, df_gtf_min_no_inside_genes, padj)

    #Combine non-sig gtf with runin
    df_plu_sense_nonsig_runin = runion_gtf_with_nonsig_gtf(df_plu_antisense_nonsig_gtf, df_plu_ADOG_gtf)
    df_min_sense_nonsig_runin = runion_gtf_with_nonsig_gtf(df_min_antisense_nonsig_gtf, df_plu_ADOG_gtf)

    df_plu_sense_nonsig_runin.to_csv(output_prefix + "plu_ALL_SAMPLES_ADOG_with_nonsig.gtf", sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)
    df_min_sense_nonsig_runin.to_csv(output_prefix + "min_ALL_SAMPLES_ADOG_with_nonsig.gtf", sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)

    gtf_list = [
    output_prefix + "plu_ALL_SAMPLES_DOG_with_nonsig.gtf", \
    output_prefix + "min_ALL_SAMPLES_DOG_with_nonsig.gtf", \
    output_prefix + "plu_ALL_SAMPLES_ADOG_with_nonsig.gtf", \
    output_prefix + "min_ALL_SAMPLES_ADOG_with_nonsig.gtf", \
    ]



    #################################################################
    #Create the R scripts to run from new gtf's
    def create_multiple_gtf_R_scripts():
        for gtf in gtf_list:
            with open(input_prefix_initial + "Rsubread_initial.R", "r") as infile, open(gtf[:-4] + "Rsubread.R", "w") as outfile:
                for line in infile:
                    line = line.replace(annotation_file_with_path, str(gtf))
                    line = line.replace(input_prefix_initial + "Rsubread_sense.txt", gtf[:-4] + "_Dogcatcher_Rsubread_sense.txt")
                    line = line.replace(input_prefix_initial + "Rsubread_antisense.txt", gtf[:-4] + "_Dogcatcher_Rsubread_antisense.txt")
                    outfile.write(line)

            with open(input_prefix_initial + "DESeq2_initial.R", "r") as infile, open(gtf[:-4] + "DESeq2.R", "w") as outfile:
                for line in infile:
                    line = line.replace(annotation_file_with_path, str(gtf))
                    line = line.replace(input_prefix_initial + "Rsubread_sense_normalized_rRNA_star.txt", gtf[:-4] + "_Dogcatcher_Rsubread_sense_normalized_rRNA_star.txt")
                    line = line.replace(input_prefix_initial + "Rsubread_antisense_normalized_rRNA_star.txt", gtf[:-4] + "_Dogcatcher_Rsubread_antisense_normalized_rRNA_star.txt")
                    line = line.replace(input_prefix_initial + "DESeq2_sense.csv", gtf[:-4] + "_Dogcatcher_DESeq2_sense.csv")
                    line = line.replace(input_prefix_initial + "DESeq2_antisense.csv", gtf[:-4] + "_Dogcatcher_DESeq2_antisense.csv")
                    outfile.write(line)
    create_multiple_gtf_R_scripts()

