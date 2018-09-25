#!/usr/bin/env python
import argparse
import csv
import numpy as np
import os
import sys
import pandas as pd
import string
import re







# Running in parallel
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Part V: Getting sig DOG/ADOG and matching back with overlapping genes.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--annotation_file_with_path', action= 'store', metavar='annotation_file_with_path', help='ensembl or wormbase annotation file for your organism.  The file may be downloaded from http://www.ensembl.org/info/data/ftp/index.html for ensembl or ftp://ftp.wormbase.org/pub/wormbase/releases/ for wormbase')
    parser.add_argument('--input_DOG_DESeq2_plu_sense_file', action= 'store', metavar='input_DOG_DESeq2_plu_sense_file', help='Run on Raw Sense output file from DSeq2 (Must have no column header for first column). Example 0.0_input_files/')
    parser.add_argument('--input_DOG_DESeq2_min_sense_file', action= 'store', metavar='input_DOG_DESeq2_min_sense_file', help='Run on Raw Sense output file from DSeq2 (Must have no column header for first column). Example 0.0_input_files/')
    parser.add_argument('--input_ADOG_DESeq2_plu_antisense_file', action= 'store', metavar='input_ADOG_DESeq2_plu_sense_file', help='Run on Raw Sense output file from DSeq2 (Must have no column header for first column). Example 0.0_input_files/')
    parser.add_argument('--input_ADOG_DESeq2_min_antisense_file', action= 'store', metavar='input_ADOG_DESeq2_min_sense_file', help='Run on Raw Sense output file from DSeq2 (Must have no column header for first column). Example 0.0_input_files/')
    parser.add_argument('--output_prefix', action= 'store', metavar='output_prefix', default = "Dogcatcher_out/", help='Output prefix. KEEP THIS THE SAME FROM PART II: Default: Dogcatcher_out/')
    parser.add_argument('--input_prefix', action= 'store', metavar='input_prefix', default = "Dogcatcher_out/", help='Output prefix. KEEP THIS THE SAME FROM PART II: Default: Dogcatcher_out/')
    parser.add_argument('--padj', action= 'store', dest='padj', metavar='padj', default= 0.05, type=float, help='DSeq2 padj value cutoff to select non-significant genes. Default: "0.05"')
    parser.add_argument('--log2FoldChange', action= 'store', dest='log2FoldChange', metavar='log2FoldChange', default= 0.0, type=float, help='DSeq2 padj value cutoff to select non-significant genes. Default: "0.05"')
    args = parser.parse_args()

    ################################################ parameters

    annotation_file_with_path = args.annotation_file_with_path
    annotation_file = os.path.basename(annotation_file_with_path)
    annotation_path = os.path.dirname(annotation_file_with_path)

    input_DOG_DESeq2_plu_sense_file = args.input_DOG_DESeq2_plu_sense_file
    input_DOG_DESeq2_min_sense_file = args.input_DOG_DESeq2_min_sense_file
    input_ADOG_DESeq2_plu_antisense_file = args.input_ADOG_DESeq2_plu_antisense_file
    input_ADOG_DESeq2_min_antisense_file = args.input_ADOG_DESeq2_min_antisense_file

    padj = args.padj
    log2FoldChange = args.log2FoldChange
    strpadj = str(padj)
    strlog2FoldChange = str(log2FoldChange)
    output_prefix = args.output_prefix + "/"
    input_prefix = args.input_prefix + "/"



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


    print("Starting Part 5: filtering significant DOGS/ADOGS")
    print("Creating output folders...")
    def make_output_folders_from_list(folder_list):
        for folder in folder_list:
            directory = output_prefix + folder
            if not os.path.exists(directory):
                os.makedirs(directory)


    greater_output = "less_than_padj_" + strpadj + "/greater_than_log2FoldChange_" + strlog2FoldChange + "/"
    greater_output_name = "less_than_padj_" + strpadj + "_greater_than_log2FoldChange_" + strlog2FoldChange
    less_output = "less_than_padj_" + strpadj + "/less_than_log2FoldChange_" + strlog2FoldChange + "/"
    less_output_name = "less_than_padj_" + strpadj + "_less_than_log2FoldChange_" + strlog2FoldChange
    make_output_folders_from_list(["/Significant/DOG/" + greater_output,
                                   "/Significant/DOG/" + less_output,
                                   "/Significant/DOG/biotypes/" + greater_output,
                                   "/Significant/DOG/biotypes/" + less_output,
                                   "/Significant/DOG/bed/" + greater_output,
                                   "/Significant/DOG/bed/" + less_output,
                                   "/Significant/ADOG/" + greater_output,
                                   "/Significant/ADOG/" + less_output,
                                   "/Significant/ADOG/biotypes/" + greater_output,
                                   "/Significant/ADOG/biotypes/" + less_output,
                                   "/Significant/ADOG/bed/" + greater_output,
                                   "/Significant/ADOG/bed/" + less_output,
                                   "/All/DOG/",
                                   "/All/DOG/biotypes/",
                                   "/All/DOG/bed/",
                                   "/All/ADOG/",
                                   "/All/ADOG/biotypes/",
                                   "/All/ADOG/bed/",
                                   ])

    #List for combined outputs
    #RUNON
    combined_DOG_list = ["gene_id_name",
    "gene_biotype_name",
    "TYPE",
    "chr",
    "strand",
    "DOG_length",
    "DOG_start_local",
    "DOG_end_local",
    "DOG_start_meta",
    "DOG_end_meta",
    "DOG_from_downstream_gene",
    "DOG_from_upstream_gene",
    "DOG_into_downstream_gene",
    "POG_into_upstream_gene",
    "padj", "baseMean", "pvalue", "log2FoldChange", "lfcSE", "stat",
    "source",
    "length",
    "left_int_length",
    "right_int_length",
    "sample_name",
    "start",
    "end",
    "type"]

    combined_ADOG_list = combined_DOG_list

    #########HELPER FUNCTIONS
    def filter_list_for_string(input_list, str):
        """This function will filter a list for the string of choice"""
        r = re.compile(".*" + str)
        newlist = list((filter(r.match, input_list)))
        return newlist


        #     DOG_overlap_list = filter_list_for_string(cols, "_overlap_")
        #     combined_list = DOG_overlap_list + inside_gene_list
        #     #no_biotype_list = list(set(cols) - set(combined_list) ) # Use sets to find differences in lists
        #     df_biotypes_sorted = df[combined_DOG_list +  DOG_overlap_list + inside_gene_list]
        # if DOG_or_ADOG == "ADOG":
        #     inside_gene_list = filter_list_for_string(cols, "inside_gene_")
        #     df_biotypes_sorted = df[combined_DOG_list +  inside_gene_list]
        # return df_no_biotypes, df_biotypes_sorted

    def get_unpacked_inside_genes(df_gtf, names_biotype):
            """This function will loop through all biotypes,
            get rid of nan values, unpack the list of genes
            return as a single col df"""
            final_inside_genes_list = []
            for biotype in names_biotype:
                inside_genes = df_gtf['DOG_overlap_opposite_strand_local.' + biotype].astype(object).tolist() #Make list of inside genes
                inside_genes = [x for x in inside_genes if str(x) != 'nan']          #Remove nan
                inside_genes = [i.strip('[]') for i in inside_genes]
                inside_genes = [i.strip() for i in inside_genes]
                inside_genes = [w.replace(" ", "",) for w in inside_genes]
                inside_genes = [w.replace("\'", "",) for w in inside_genes]
                final_inside_genes_list.extend(inside_genes)
            #create a dataframe of the inside genes, some columns will have more than one gene
            df = pd.DataFrame({'tags': final_inside_genes_list})
            #Make a new list out of the columns
            all_genes_list = []
            for x in range(len(df)):
                list = df["tags"].iloc[x].split(",")
                all_genes_list.extend(list)
            df = pd.DataFrame({'inside_gene_same_strand_unpacked': all_genes_list})
            return df


    def get_unpacked_csv(dseq_input, df_runion_csv, DOG_or_ADOG):
        """This function will take in a DSeq2 normalized matrix and match this with your runion output"""
        df_dseq = pd.read_csv(dseq_input, index_col=None, sep=',')
        df_dseq.rename(columns={"Unnamed: 0" : "gene_id_name"}, inplace=True) #Rename unnamed column of deseq2
        df = df_dseq.merge(df_runion_csv, left_on="gene_id_name", right_on='gene_id_name', how='inner')
        ########################
        del df["gene_id"]
        del df["dot"]
        del df["dot2"]
        #############################
        #BIOTYPE WRANGLING
        cols = list(df.columns.values)
        #DOG (Will have overlapping biotypes)
        df_final = pd.DataFrame(columns=['gene_id_name',"strand","gene_biotype_name","overlapped_gene","overlapped_biotype","same_or_opposite_strand"])
        if DOG_or_ADOG == "DOG":
            inside_gene_list = filter_list_for_string(cols, "inside_gene_")
            DOG_overlap_same_list = filter_list_for_string(cols, "DOG_overlap_same_strand_local:")
            for col in DOG_overlap_same_list:
                x = 0
                for x in range(len(df)):
                    cell = df[col].iloc[x]
                    try:
                        df[col].iloc[x].astype(float)
                    except:
                        if "[" in cell:
                            cell = cell.replace("[","").replace("]","").replace(" ","").replace("\'","")
                            celllist = cell.split(",")
                            for gene in celllist:
                                df_final = df_final.append({"gene_id_name": df["gene_id_name"].iloc[x],
                                                 "strand": df["strand"].iloc[x],
                                                 "gene_biotype_name": df["gene_biotype_name"].iloc[x],
                                                 "overlapped_gene": gene,
                                                 "overlapped_biotype":col[len("DOG_overlap_same_strand_local:"):],
                                                 "same_or_opposite_strand":"same"}, ignore_index=True)
                    x = x + 1
            DOG_overlap_opposite_list = filter_list_for_string(cols, "DOG_overlap_opposite_strand_local:")
            for col in DOG_overlap_opposite_list:
                x = 0
                for x in range(len(df)):
                    cell = df[col].iloc[x]
                    try:
                        df[col].iloc[x].astype(float)
                    except:
                        if "[" in cell:
                            cell = cell.replace("[","").replace("]","").replace(" ","").replace("\'","")
                            celllist = cell.split(",")
                            for gene in celllist:
                                df_final = df_final.append({"gene_id_name": df["gene_id_name"].iloc[x],
                                                 "gene_biotype_name": df["gene_biotype_name"].iloc[x],
                                                 "strand": df["strand"].iloc[x],
                                                 "overlapped_gene": gene,
                                                 "overlapped_biotype":col[len("DOG_overlap_opposite_strand_local:"):],
                                                 "same_or_opposite_strand":"opposite"}, ignore_index=True)
                    x = x + 1
                    

                    
        return df_final



    ##############################################
    ##############################################
    ### Create gtf of run on/in contigs with      non-significant genes (without run on/in genes)
    #load original gtf

    def get_merged_csv(dseq_input, df_runion_csv, DOG_or_ADOG):
        """This function will take in a DSeq2 normalized matrix and match this with your runion output"""
        df_dseq = pd.read_csv(dseq_input, index_col=None, sep=',')
        df_dseq.rename(columns={"Unnamed: 0" : "gene_id_name"}, inplace=True) #Rename unnamed column of deseq2
        df = df_dseq.merge(df_runion_csv, left_on="gene_id_name", right_on='gene_id_name', how='inner')

        ########################
        del df["gene_id"]
        del df["dot"]
        del df["dot2"]
        #############################
        #BIOTYPE WRANGLING
        cols = list(df.columns.values)
        #DOG (Will have overlapping biotypes)
        df_no_biotypes = df[combined_DOG_list]
        if DOG_or_ADOG == "DOG":
            inside_gene_list = filter_list_for_string(cols, "inside_gene_")
            DOG_overlap_list = filter_list_for_string(cols, "_overlap_")
            combined_list = DOG_overlap_list + inside_gene_list
            #no_biotype_list = list(set(cols) - set(combined_list) ) # Use sets to find differences in lists
            df_biotypes_sorted = df[combined_DOG_list +  DOG_overlap_list + inside_gene_list]
        if DOG_or_ADOG == "ADOG":
            inside_gene_list = filter_list_for_string(cols, "inside_gene_")
            df_biotypes_sorted = df[combined_DOG_list +  inside_gene_list]
        return df_no_biotypes, df_biotypes_sorted

    #############################################################
    #############################################################
    #############################################################
    #############################################################
    #ALL RUNON

        #Lists to output bed files
    bed_list = ["chr", "DOG_start_local","DOG_end_local","gene_id_name"]



    df_no_biotypes_DESeq2_plu_sense_file, df_biotypes_sorted_DESeq2_plu_sense_file = get_merged_csv(input_DOG_DESeq2_plu_sense_file, df_plu_DOG_csv, "DOG")
    df_no_biotypes_DESeq2_min_sense_file, df_biotypes_sorted_DESeq2_min_sense_file = get_merged_csv(input_DOG_DESeq2_min_sense_file, df_min_DOG_csv, "DOG")

                            #Make bed files
    df_no_biotypes_DESeq2_plu_sense_file_bed = df_no_biotypes_DESeq2_plu_sense_file[bed_list]
    df_no_biotypes_DESeq2_min_sense_file_bed = df_no_biotypes_DESeq2_min_sense_file[bed_list]
    df_no_biotypes_DESeq2_plu_sense_file_bed.to_csv(output_prefix + "/All/DOG/bed/" + "plu_sense_DOG.bed", sep="\t", index=None,header=None)
    df_no_biotypes_DESeq2_min_sense_file_bed.to_csv(output_prefix + "/All/DOG/bed/" + "min_sense_DOG.bed", sep="\t", index=None,header=None)

    #Combine_DOG plus and minus strands
    combined_DOG = pd.concat([df_no_biotypes_DESeq2_plu_sense_file, df_no_biotypes_DESeq2_min_sense_file])
    combined_DOG = combined_DOG [combined_DOG_list]
    combined_DOG.to_csv(output_prefix + "/All/DOG/" + "combined_DOG.csv", sep="\t", index=None)
    combined_DOG_biotypes = pd.concat([df_biotypes_sorted_DESeq2_plu_sense_file, df_biotypes_sorted_DESeq2_min_sense_file])
    combined_DOG_biotypes.to_csv(output_prefix + "/All/DOG/biotypes/" + "combined_DOG_with_biotypes.csv", sep="\t", index=None)

    dfunpackedplu = get_unpacked_csv(input_DOG_DESeq2_plu_sense_file, df_plu_DOG_csv, "DOG")
    dfunpackedmin = get_unpacked_csv(input_DOG_DESeq2_min_sense_file, df_min_DOG_csv, "DOG")
    dfunpacked = pd.concat([dfunpackedplu,dfunpackedmin])
    dfunpacked.to_csv(output_prefix + "/All/DOG/biotypes/" + "combined_DOG_with_biotypes_UNPACKED.csv", sep="\t", index=None)


    #All RUNIN
    df_plu_ADOG, df_plu_ADOG_biotypes  = get_merged_csv(input_ADOG_DESeq2_plu_antisense_file, df_plu_ADOG_csv,"ADOG")
    df_plu_ADOG_bed = df_plu_ADOG [ bed_list]
    df_plu_ADOG_bed.to_csv(output_prefix + "/All/ADOG/bed/" + "plu_antisense_ADOG.bed", sep="\t", index=None,   header=None)

    df_min_ADOG, df_min_ADOG_biotypes  = get_merged_csv(input_ADOG_DESeq2_min_antisense_file, df_min_ADOG_csv,"ADOG")
    df_min_ADOG_bed = df_min_ADOG [ bed_list]
    df_min_ADOG_bed.to_csv(output_prefix + "/All/ADOG/bed/" + "min_antisense_ADOG.bed", sep="\t", index=None,   header=None)

    combined_ADOG = pd.concat([df_plu_ADOG, df_min_ADOG])
    combined_ADOG = combined_ADOG [combined_ADOG_list]
    combined_ADOG.to_csv(output_prefix + "/All/ADOG/" + "combined_ADOG.csv", sep="\t", index=None)
    combined_ADOG_biotypes = pd.concat([df_plu_ADOG_biotypes, df_min_ADOG_biotypes])
    combined_ADOG.to_csv(output_prefix + "/All/ADOG/biotypes/" + "combined_ADOG_with_biotypes.csv", sep="\t", index=None)


    #############################################################
    def get_sig_runion_csv(df, padj, log2FoldChange):
        """This function will take in a DSeq2 normalized matrix and make a gtf of all non-sig genes"""
        df["significant"] = np.where(df["padj"] < padj, True, False)
        df_final = df[df.significant]  #Keep all True values. (significant)
        del df_final["significant"]
        df_final = df_final.copy()
        #Treatment vs Control

        df_final["T_vs_C_up"] = np.where(df_final["log2FoldChange"] > log2FoldChange, True, False)
        df_T_vs_C_up = df_final[df_final.T_vs_C_up]  #Keep all True values. (up)
        del df_T_vs_C_up["T_vs_C_up"]

        df_final["T_vs_C_down"] = np.where(df_final["log2FoldChange"] < log2FoldChange, True, False)
        df_T_vs_C_down = df_final[df_final.T_vs_C_down]  #Keep all True values. (down)
        del df_T_vs_C_down["T_vs_C_up"]
        del df_T_vs_C_down["T_vs_C_down"]
        return df_T_vs_C_up, df_T_vs_C_down


    #SIGNIFICANT and 2LFC RUNON
    # #Get sig genes from deseq2
    df_plu_sense_sig_T_vs_C_up, df_plu_sense_sig_T_vs_C_down = get_sig_runion_csv(df_no_biotypes_DESeq2_plu_sense_file, padj, log2FoldChange)
    df_min_sense_sig_T_vs_C_up, df_min_sense_sig_T_vs_C_down = get_sig_runion_csv(df_no_biotypes_DESeq2_min_sense_file, padj, log2FoldChange)

    combined_sense_up = pd.concat([df_plu_sense_sig_T_vs_C_up, df_min_sense_sig_T_vs_C_up])
    combined_sense_down = pd.concat([df_plu_sense_sig_T_vs_C_down, df_min_sense_sig_T_vs_C_down])

    combined_sense_up = combined_sense_up [combined_DOG_list]
    combined_sense_up.to_csv(output_prefix + "/Significant/DOG/" + greater_output + greater_output_name + "_combined_DOG.csv", sep="\t", index=None)
    combined_sense_down = combined_sense_down [combined_DOG_list]
    combined_sense_down.to_csv(output_prefix + "/Significant/DOG/" + less_output + less_output_name + "_combined_DOG.csv", sep="\t", index=None)

    #Make bed files

    df_plu_sense_sig_T_vs_C_up_bed   = df_plu_sense_sig_T_vs_C_up [  bed_list]
    df_plu_sense_sig_T_vs_C_down_bed = df_plu_sense_sig_T_vs_C_down [bed_list]
    df_min_sense_sig_T_vs_C_up_bed   = df_min_sense_sig_T_vs_C_up [  bed_list]
    df_min_sense_sig_T_vs_C_down_bed = df_min_sense_sig_T_vs_C_down [bed_list]
    
    df_plu_sense_sig_T_vs_C_up_bed.to_csv(output_prefix +   "/Significant/DOG/bed/" + greater_output + greater_output_name + "_plu_sense_DOG.bed", sep="\t", index=None, header=None)
    df_plu_sense_sig_T_vs_C_down_bed.to_csv(output_prefix + "/Significant/DOG/bed/" + less_output + less_output_name + "_plu_sense_DOG.bed", sep="\t", index=None, header=None)
    df_min_sense_sig_T_vs_C_up_bed.to_csv(output_prefix +   "/Significant/DOG/bed/" + greater_output + greater_output_name + "_min_sense_DOG.bed", sep="\t", index=None, header=None)
    df_min_sense_sig_T_vs_C_down_bed.to_csv(output_prefix + "/Significant/DOG/bed/" + less_output + less_output_name + "_min_sense_DOG.bed", sep="\t", index=None, header=None)

    #Biotypes
    df_plu_sense_sig_biotypes_T_vs_C_up, df_plu_sense_sig_biotypes_T_vs_C_down = get_sig_runion_csv(df_biotypes_sorted_DESeq2_plu_sense_file, padj, log2FoldChange)
    df_min_sense_sig_biotypes_T_vs_C_up, df_min_sense_sig_biotypes_T_vs_C_down = get_sig_runion_csv(df_biotypes_sorted_DESeq2_min_sense_file, padj, log2FoldChange)
    combined_sense_up_biotypes = pd.concat([df_plu_sense_sig_biotypes_T_vs_C_up, df_min_sense_sig_biotypes_T_vs_C_up])
    combined_sense_down_biotypes = pd.concat([df_plu_sense_sig_biotypes_T_vs_C_down, df_min_sense_sig_biotypes_T_vs_C_down])


    dfunpackedup = dfunpacked[dfunpacked["gene_id_name"].isin(combined_sense_up["gene_id_name"])]
    dfunpackedup.to_csv(output_prefix + "/Significant/DOG/biotypes/" + greater_output + greater_output_name + "combined_DOG_with_biotypes_UNPACKED.csv", sep="\t", index=None)
    dfunpackeddown = dfunpacked[dfunpacked["gene_id_name"].isin(combined_sense_down["gene_id_name"])]
    dfunpackeddown.to_csv(output_prefix + "/Significant/DOG/biotypes/" + less_output + less_output_name + "combined_DOG_with_biotypes_UNPACKED.csv", sep="\t", index=None)


    combined_sense_up_biotypes.to_csv(output_prefix + "/Significant/DOG/biotypes/" + greater_output + greater_output_name + "_combined_DOG_with_biotypes.csv", sep="\t", index=None)
    combined_sense_down_biotypes.to_csv(output_prefix + "/Significant/DOG/biotypes/" + less_output + less_output_name + "_combined_DOG_with_biotypes.csv", sep="\t", index=None)


    ############################################################
    ############################################################
    ############################################################
    ############################################################
        #SIGNIFICANT and 2LFC RUNIN
    # #Get sig genes from deseq2
    df_plu_antisense_sig_T_vs_C_up, df_plu_antisense_sig_T_vs_C_down = get_sig_runion_csv(df_plu_ADOG, padj, log2FoldChange)
    df_min_antisense_sig_T_vs_C_up, df_min_antisense_sig_T_vs_C_down = get_sig_runion_csv(df_min_ADOG, padj, log2FoldChange)
    combined_antisense_up = pd.concat([df_plu_antisense_sig_T_vs_C_up, df_min_antisense_sig_T_vs_C_up])
    combined_antisense_down = pd.concat([df_plu_antisense_sig_T_vs_C_down, df_min_antisense_sig_T_vs_C_down])
    #Write to csv
    df_plu_antisense_sig_T_vs_C_up.to_csv(output_prefix + "/Significant/ADOG/" + greater_output + greater_output_name + "_plu_antisense_ADOG.csv", sep="\t", index=None)
    df_plu_antisense_sig_T_vs_C_down.to_csv(output_prefix + "/Significant/ADOG/" + less_output + less_output_name + "_plu_antisense_ADOG.csv", sep="\t", index=None)
    df_min_antisense_sig_T_vs_C_up.to_csv(output_prefix + "/Significant/ADOG/" + greater_output + greater_output_name + "_min_antisense_ADOG.csv", sep="\t", index=None)
    df_min_antisense_sig_T_vs_C_down.to_csv(output_prefix + "/Significant/ADOG/" + less_output + less_output_name + "_min_antisense_ADOG.csv", sep="\t", index=None)


    combined_antisense_up = combined_antisense_up[combined_ADOG_list]
    combined_antisense_up.to_csv(output_prefix + "/Significant/ADOG/" + greater_output + greater_output_name + "_combined_ADOG.csv", sep="\t", index=None)
    combined_antisense_down = combined_antisense_down[combined_ADOG_list]
    combined_antisense_down.to_csv(output_prefix + "/Significant/ADOG/" + less_output + less_output_name + "_combined_ADOG.csv", sep="\t", index=None)

    #Make bed files

    df_plu_antisense_sig_T_vs_C_up_bed   = df_plu_antisense_sig_T_vs_C_up [   bed_list]
    df_plu_antisense_sig_T_vs_C_down_bed = df_plu_antisense_sig_T_vs_C_down [ bed_list]
    df_min_antisense_sig_T_vs_C_up_bed   = df_min_antisense_sig_T_vs_C_up [   bed_list]
    df_min_antisense_sig_T_vs_C_down_bed = df_min_antisense_sig_T_vs_C_down [ bed_list]
    
    df_plu_antisense_sig_T_vs_C_up_bed.to_csv(output_prefix +   "/Significant/ADOG/bed/" + greater_output + greater_output_name + "_plu_antisense_ADOG.bed", sep="\t", index=None, header=None)
    df_plu_antisense_sig_T_vs_C_down_bed.to_csv(output_prefix + "/Significant/ADOG/bed/" + less_output + less_output_name + "_plu_antisense_ADOG.bed", sep="\t", index=None, header=None)
    df_min_antisense_sig_T_vs_C_up_bed.to_csv(output_prefix +   "/Significant/ADOG/bed/" + greater_output + greater_output_name + "_min_antisense_ADOG.bed", sep="\t", index=None, header=None)
    df_min_antisense_sig_T_vs_C_down_bed.to_csv(output_prefix + "/Significant/ADOG/bed/" + less_output + less_output_name + "_min_antisense_ADOG.bed", sep="\t", index=None, header=None)

    #Biotypes
    df_plu_antisense_sig_biotypes_T_vs_C_up, df_plu_antisense_sig_biotypes_T_vs_C_down = get_sig_runion_csv(df_plu_ADOG_biotypes, padj, log2FoldChange)
    df_min_antisense_sig_biotypes_T_vs_C_up, df_min_antisense_sig_biotypes_T_vs_C_down = get_sig_runion_csv(df_min_ADOG_biotypes, padj, log2FoldChange)
    combined_antisense_up_biotypes = pd.concat([df_plu_antisense_sig_biotypes_T_vs_C_up, df_min_antisense_sig_biotypes_T_vs_C_up])
    combined_antisense_down_biotypes = pd.concat([df_plu_antisense_sig_biotypes_T_vs_C_down, df_min_antisense_sig_biotypes_T_vs_C_down])

    combined_antisense_up_biotypes.to_csv(output_prefix + "/Significant/ADOG/biotypes/" + greater_output + greater_output_name + "_combined_ADOG_with_biotypes.csv", sep="\t", index=None)
    combined_antisense_down_biotypes.to_csv(output_prefix + "/Significant/ADOG/biotypes/" + less_output + less_output_name + "_combined_ADOG_with_biotypes.csv", sep="\t", index=None)
