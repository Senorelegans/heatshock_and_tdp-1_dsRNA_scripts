import pandas as pd
import argparse
import csv
import os
import numpy as np
import string

print("Parsing gtf and putting output into gtf folder...")
parser = argparse.ArgumentParser(description='Part I: Parsing the gtf to take out genes inside of genes.')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('--annotation_file_with_path', action= 'store', metavar='annotation_file_with_path', help='ensembl or wormbase annotation file for your organism.  The file may be downloaded from http://www.ensembl.org/info/data/ftp/index.html for ensembl or ftp://ftp.wormbase.org/pub/wormbase/releases/ for wormbase')
args = parser.parse_args()

annotation_file_with_path = args.annotation_file_with_path

my_col = ["chr","source","type","start","end","dot","strand","dot2","gene_id"]
df_gtf = pd.read_csv(annotation_file_with_path, sep="\t",header=None,names=my_col, comment="#", dtype={"chr" : str, "source" : str,"type" : str, "start" : int,"end" : int,"dot" : str,"strand" : str, "dot2" : str,"gene_id" : str})
df_gtf = df_gtf[~df_gtf["chr"].str.contains("\.") ]    # Take out patches
df_gtf.to_csv(annotation_file_with_path[:-4] + "_NO_PATCHES.gtf", sep="\t", index=None,quoting=csv.QUOTE_NONE)
df_gtf = df_gtf[df_gtf.type.str.contains("gene")]      # Keep everything that says gene

df_gtf["length"] = df_gtf["end"] - df_gtf["start"]
df_gtf = df_gtf.dropna()
chr_names=df_gtf['chr'].unique().tolist()


def header_int_to_char(df):
    """this function will change all numbers to letters for pandas"""
    import string
    x = 0  #rename int headers to letters
    for num in (list(df)):
        df.rename(columns={int(num) : str(string.ascii_uppercase[x])}, inplace=True)
        x = x + 1

def get_biotypes_names(df_gtf):
    """This function will take a gtf and unpack the gene id and biotypes column"""
    #Expand the gene id column, then merge back to get gene_id name in new column
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

    df_merge["gene_biotype_name"] = (df_merge[biotype_header]
                        .str
                        .extract(p, expand=True))

    #Some protein coding truncated end, change nan's
    x = 0  #rename int headers to letters
    for col in (list(df_merge)):
        for x in range(len(string.ascii_uppercase)):
            if col == string.ascii_uppercase[x]:
                del df_merge[col]

    df_merge = df_merge.replace(np.nan, "protein_coding", regex=True)
    names_biotype = df_merge["gene_biotype_name"].unique().tolist()
    return df_merge, names_biotype

def get_chrm_dic(df_gtf, chr_names, plu_or_min, gene_biotype):

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

def get_biotypes_inside_gene(df_biotypes, gtf_start, gtf_end):
    try:
        rows_df = df_biotypes[ (gtf_start < df_biotypes["start"]) & (df_biotypes["end"] < gtf_end) ]
        if len(rows_df) == 0:
            return "NA"
        else:
            return rows_df["gene_id_name"].tolist()
    except:
        print("*** Biotype overlap NA ***")


df_gtf, names_biotype = get_biotypes_names(df_gtf)
df_gtf.to_csv(annotation_file_with_path[:-4] + "_ALL_GENES_NO_PATCHES.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)
df_gtf_plu = df_gtf[ df_gtf["strand"].str.contains("\+") ].copy()
df_gtf_min = df_gtf[ df_gtf["strand"].str.contains("-") ].copy()
print(names_biotype)

strand_list = ["+", "-"]
for strand in strand_list:
    print("strand is", strand)
    for biotype in names_biotype:
        print("*"*20 + "\n")
        print("working on biotype:" + biotype)

        if strand == "+":
            d_same_strand_gtf_chr =     get_chrm_dic(df_gtf, chr_names, "+", biotype)
            d_opposite_strand_gtf_chr = get_chrm_dic(df_gtf, chr_names, "-", biotype)
            if d_same_strand_gtf_chr != "NA":
                print("working on inside_gene_same_strand...")
                df_gtf_plu['inside_gene_same_strand:' + biotype] = df_gtf_plu.apply(lambda row: get_biotypes_inside_gene(d_same_strand_gtf_chr[row['chr']], row["start"], row["end"]), axis=1)
            if d_opposite_strand_gtf_chr != "NA":
                print("working on inside_gene_opposite_strand...")
                df_gtf_plu['inside_gene_opposite_strand:' + biotype] = df_gtf_plu.apply(lambda row: get_biotypes_inside_gene(d_opposite_strand_gtf_chr[row['chr']], row["start"], row["end"]), axis=1)

        if strand == "-":
            d_same_strand_gtf_chr = get_chrm_dic(df_gtf, chr_names, "-", biotype)
            d_opposite_strand_gtf_chr = get_chrm_dic(df_gtf, chr_names, "+", biotype)
            if d_same_strand_gtf_chr != "NA":
                print("working on inside_gene_same_strand...")
                df_gtf_min['inside_gene_same_strand:' + biotype] = df_gtf_min.apply(lambda row: get_biotypes_inside_gene(d_same_strand_gtf_chr[row['chr']], row["start"], row["end"]), axis=1)
            if d_opposite_strand_gtf_chr != "NA":
                print("working on inside_gene_opposite_strand...")
                df_gtf_min['inside_gene_opposite_strand:' + biotype] = df_gtf_min.apply(lambda row: get_biotypes_inside_gene(d_opposite_strand_gtf_chr[row['chr']], row["start"], row["end"]), axis=1)

#Write out total gtf
df_gtf_plu.to_csv(annotation_file_with_path[:-4] + "_plu_ALL_with_inside_genes.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)
df_gtf_min.to_csv(annotation_file_with_path[:-4] + "_min_ALL_with_inside_genes.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)

df_gtf_plu = pd.read_csv(annotation_file_with_path[:-4] + "_plu_ALL_with_inside_genes.txt", sep="\t", comment="#")
df_gtf_min = pd.read_csv(annotation_file_with_path[:-4] + "_min_ALL_with_inside_genes.txt", sep="\t", comment="#")

def get_unpacked_inside_genes_same_strand(df_gtf, names_biotype):
    """This function will loop through all biotypes,
    get rid of nan values, unpack the list of genes
    return as a single col df"""
    final_inside_genes_list = []
    for biotype in names_biotype:
        try:
            inside_genes = df_gtf['inside_gene_same_strand:' + biotype].astype(object).tolist() #Make list of inside genes
            inside_genes = [x for x in inside_genes if str(x) != 'nan']          #Remove nan
            inside_genes = [i.strip('[]') for i in inside_genes]
            inside_genes = [i.strip() for i in inside_genes]
            inside_genes = [w.replace(" ", "",) for w in inside_genes]
            inside_genes = [w.replace("\'", "",) for w in inside_genes]
            final_inside_genes_list.extend(inside_genes)
        except:
            print("No genes for inside_gene_same_strand: ", biotype)
    #create a dataframe of the inside genes, some columns will have more than one gene
    df = pd.DataFrame({'tags': final_inside_genes_list})
    #Make a new list out of the columns
    all_genes_list = []
    for x in range(len(df)):
        list = df["tags"].iloc[x].split(",")
        all_genes_list.extend(list)
    df = pd.DataFrame({'inside_gene_same_strand_unpacked': all_genes_list})
    return df

print("Finished finding all genes inside other genes...")
#-------------------------------------------------------#

def get_unpacked_inside_genes_opposite_strand(df_gtf, names_biotype):
    """This function will loop through all biotypes,
    get rid of nan values, unpack the list of genes
    return as a single col df"""
    final_inside_genes_list = []
    for biotype in names_biotype:
        try:
            inside_genes = df_gtf['inside_gene_opposite_strand:' + biotype].astype(object).tolist() #Make list of inside genes
            inside_genes = [x for x in inside_genes if str(x) != 'nan']          #Remove nan
            inside_genes = [i.strip('[]') for i in inside_genes]
            inside_genes = [i.strip() for i in inside_genes]
            inside_genes = [w.replace(" ", "",) for w in inside_genes]
            inside_genes = [w.replace("\'", "",) for w in inside_genes]
            final_inside_genes_list.extend(inside_genes)
        except:
            print("No genes for inside_gene_opposite_strand: ", biotype)
    #create a dataframe of the inside genes, some columns will have more than one gene
    df = pd.DataFrame({'tags': final_inside_genes_list})
    #Make a new list out of the columns
    all_genes_list = []
    for x in range(len(df)):
        list = df["tags"].iloc[x].split(",")
        all_genes_list.extend(list)
    df = pd.DataFrame({'inside_gene_opposite_strand_unpacked': all_genes_list})
    return df


print("Finished finding all genes inside other genes...")
#-------------------------------------------------------#

df_gtf_plu_inside_genes = get_unpacked_inside_genes_same_strand(df_gtf_plu, names_biotype)
df_gtf_plu_inside_genes.to_csv(annotation_file_with_path[:-4] + "_plu_inside_genes_same_strand.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)

df_gtf_plu_NO_inside_genes = df_gtf_plu[~df_gtf_plu.gene_id_name.isin(df_gtf_plu_inside_genes.inside_gene_same_strand_unpacked)]
df_gtf_plu_NO_inside_genes.to_csv(annotation_file_with_path[:-4] + "_plu_NO_inside_genes_same_strand_with_overlap.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)

df_gtf_min_inside_genes = get_unpacked_inside_genes_same_strand(df_gtf_min, names_biotype)
df_gtf_min_inside_genes.to_csv(annotation_file_with_path[:-4] + "_min_inside_genes_same_strand.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)

df_gtf_min_NO_inside_genes = df_gtf_min[~df_gtf_min.gene_id_name.isin(df_gtf_min_inside_genes.inside_gene_same_strand_unpacked)]
df_gtf_min_NO_inside_genes.to_csv(annotation_file_with_path[:-4] + "_min_NO_inside_genes_same_strand_with_overlap.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)

df_gtf_plu_NO_inside_genes_bed = df_gtf_plu_NO_inside_genes[['chr', "start", "end","gene_id_name"]]
df_gtf_plu_NO_inside_genes_bed.to_csv(annotation_file_with_path[:-4] + "_plu_NO_inside_genes_same_strand.BedGraph", sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)

df_gtf_min_NO_inside_genes_bed = df_gtf_min_NO_inside_genes[['chr', "start", "end","gene_id_name"]]
df_gtf_min_NO_inside_genes_bed.to_csv(annotation_file_with_path[:-4] + "_min_NO_inside_genes_same_strand.BedGraph", sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)


df_gtf_plu_NO_inside_genes_bed = df_gtf_plu_NO_inside_genes[["chr","source","type","start","end","dot","strand","dot2","gene_id"]]
df_gtf_plu_NO_inside_genes_bed.to_csv(annotation_file_with_path[:-4] + "_plu_NO_inside_genes_same_strand.gtf", sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)

df_gtf_min_NO_inside_genes_bed = df_gtf_min_NO_inside_genes[["chr","source","type","start","end","dot","strand","dot2","gene_id"]]
df_gtf_min_NO_inside_genes_bed.to_csv(annotation_file_with_path[:-4] + "_min_NO_inside_genes_same_strand.gtf", sep="\t", index=None, header=None,   quoting=csv.QUOTE_NONE)


df_gtf_plu_inside_genes = get_unpacked_inside_genes_opposite_strand(df_gtf_plu, names_biotype)
df_gtf_plu_inside_genes.to_csv(annotation_file_with_path[:-4] + "_plu_inside_genes_opposite_strand.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)

df_gtf_min_inside_genes = get_unpacked_inside_genes_opposite_strand(df_gtf_min, names_biotype)
df_gtf_min_inside_genes.to_csv(annotation_file_with_path[:-4] + "_min_inside_genes_opposite_strand.txt", sep="\t", index=None,quoting=csv.QUOTE_NONE)
print("Finished writing out gtf's...")


