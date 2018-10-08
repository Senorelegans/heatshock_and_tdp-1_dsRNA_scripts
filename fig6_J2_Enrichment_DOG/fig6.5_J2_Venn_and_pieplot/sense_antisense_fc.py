#!python3.5
# This script will get anti-sense/sense ratios. Sort for the highest.
#load pandas and numpy and regrex
import pandas as pd
import numpy as np
import os
import sys
import re
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles





f1 = "HS_combined_DOG_noOperon_unique_cleaned.csv"
f2 = "HS_combined_ADOG_noOperon_unique_cleaned.csv"

f3 = "OK_combined_DOG_noOperon_unique_cleaned.csv"
f4 = "OK_combined_ADOG_noOperon_unique_cleaned.csv"



fN = f1[:-4]
f1out = fN + "_less.05padj.txt"

def get_above_below_2lfc_p05(f1):
    """This function will take in a DSeq2 normalized matrix and filter for everything LT 0.05
    Give back two df of above 2lfc and below 2lfc"""
    df = pd.read_csv(f1, index_col=None, sep=',')
    print("reading : ", f1)
    #df.rename(columns={"Unnamed: 0" : "gene_id"}, inplace=True)
    # df = df[df["baseMean"] > 20]
    # df = df[df["gene_biotype_name"] != "rRNA"]
    df_padj_less_05 = df[df["padj"] < 0.05]
    df_above_2_lfc = df_padj_less_05[df_padj_less_05['log2FoldChange'] > 0].sort_values(by='log2FoldChange', ascending=False)
    df_below_2_lfc = df_padj_less_05[df_padj_less_05['log2FoldChange'] < -0].sort_values(by='log2FoldChange', ascending=True)
    return df_above_2_lfc , df_below_2_lfc



df_HS_sense_above_2_lfc, df_HS_sense_below_2_lfc = get_above_below_2lfc_p05(f1)
df_HS_anti_sense_above_2_lfc, df_HS_anti_sense_below_2_lfc = get_above_below_2lfc_p05(f2)
df_OK_sense_above_2_lfc, df_OK_sense_below_2_lfc = get_above_below_2lfc_p05(f3)
df_OK_anti_sense_above_2_lfc, df_OK_anti_sense_below_2_lfc = get_above_below_2lfc_p05(f4)

df_HS_sense_above_2_lfc.to_csv("df_HS_sense_above.csv",sep="\t",index=None)
df_HS_sense_below_2_lfc.to_csv("df_HS_sense_below.csv",sep="\t",index=None)
df_OK_sense_above_2_lfc.to_csv("df_OK_sense_above.csv",sep="\t",index=None)
df_OK_sense_below_2_lfc.to_csv("df_OK_sense_below.csv",sep="\t",index=None)
df_HS_anti_sense_above_2_lfc.to_csv("df_HS_anti_sense_above.csv",sep="\t",index=None)
df_HS_anti_sense_below_2_lfc.to_csv("df_HS_anti_sense_below.csv",sep="\t",index=None)
df_OK_anti_sense_above_2_lfc.to_csv("df_OK_anti_sense_above.csv",sep="\t",index=None)
df_OK_anti_sense_below_2_lfc.to_csv("df_OK_anti_sense_below.csv",sep="\t",index=None)


df_HS_sense_above_2_lfc_GO = df_HS_sense_above_2_lfc["GeneID"]
df_OK_sense_above_2_lfc_GO = df_OK_sense_above_2_lfc["GeneID"]
df_OK_sense_above_2_lfc_GO.to_csv("df_OK_above_gene_id_GO.csv",sep="\t",index=None,header=False)
df_HS_sense_above_2_lfc_GO.to_csv("df_HS_above_gene_id_GO.csv",sep="\t",index=None,header=False) 

#MAKE NEW DATAFRAMES
df_concat  = pd.merge(df_HS_sense_above_2_lfc, df_OK_sense_above_2_lfc, how="inner", on="GeneID")
df_HS_only = df_HS_sense_above_2_lfc[~df_HS_sense_above_2_lfc.GeneID.isin(df_concat.GeneID)]
df_OK_only = df_OK_sense_above_2_lfc[~df_OK_sense_above_2_lfc.GeneID.isin(df_concat.GeneID)]




df_concat_GO = df_concat["GeneID"]
df_concat_GO.to_csv("df_HS_OK____above_gene_id_GO.csv",sep="\t",index=None)



df_concat.to_csv("df_HS_OK____sense_above.csv",sep="\t",index=None)
df_HS_only.to_csv("df_HS_only_sense_above.csv",sep="\t",index=None)
df_OK_only.to_csv("df_OK_only_sense_above.csv",sep="\t",index=None)


# Subset sizes
s = (
    len(df_HS_only),  # Ab
    len(df_OK_only),  # aB
    len(df_concat),  # AB
)

v = venn2(subsets=s, set_labels=('', ''))

for text in v.set_labels:
    text.set_x(text.get_position()[0] + 0.15)    #Move along x
    text.set_y(text.get_position()[1] + 0.15)    #Move along y
    text.set_fontsize(18)
#label.set_family('serif')
#label.set_x(label.get_position()[0] + 0.1)

#v.get_label_by_id('A').set_text('$x^2$') # Those are set labels
#v.get_label_by_id('A').set_fontsize(22)


hs_label = "Heatshock\n" + str(len(df_HS_only))

# Subset labels
v.get_label_by_id('10').set_text(str(len(df_HS_only)))
v.get_label_by_id('10').set_fontsize(20)
v.get_label_by_id('01').set_text(str(len(df_OK_only)))
v.get_label_by_id('11').set_text(str(len(df_concat)))

#v.get_label_by_id('tdp-1(ok803)').set_fontsize(20)




# Subset colors
v.get_patch_by_id('10').set_color('red')
v.get_patch_by_id('01').set_color('gray')
v.get_patch_by_id('11').set_color('white')

# Subset alphas
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)

# Border styles
c = venn2_circles(subsets=s, linestyle='solid')
c[0].set_ls('dashed')  # Line style
c[0].set_lw(2.0)       # Line width
plt.title("RUNON Heatshock and tdp-1(ok803) \n Sense J2 Enrichment compared to WT \n  log2FoldChange > 1 \n (Padj) < 0.05, \n basemean >20", fontsize=20)
plt.savefig('RUNON Heatshock and tdp-1(ok803) Sense J2 Enrichment compared to WT.png',bbox_inches='tight')
#savefig('foo.png', bbox_inches='tight')
#plt.show()










df_HS_anti_sense_above_2_lfc_GO = df_HS_anti_sense_above_2_lfc["GeneID"]
df_OK_anti_sense_above_2_lfc_GO = df_OK_anti_sense_above_2_lfc["GeneID"]
df_OK_anti_sense_above_2_lfc_GO.to_csv("df_OK_above_antisense_gene_id_GO.csv",sep="\t",index=None,header=False)
df_HS_anti_sense_above_2_lfc_GO.to_csv("df_HS_above_antisense_gene_id_GO.csv",sep="\t",index=None,header=False)

#MAKE NEW DATAFRAMES
df_concat  = pd.merge(df_HS_anti_sense_above_2_lfc, df_OK_anti_sense_above_2_lfc, how="inner", on="GeneID")
df_HS_only = df_HS_anti_sense_above_2_lfc[~df_HS_anti_sense_above_2_lfc.GeneID.isin(df_concat.GeneID)]
df_OK_only = df_OK_anti_sense_above_2_lfc[~df_OK_anti_sense_above_2_lfc.GeneID.isin(df_concat.GeneID)]




df_concat_GO = df_concat["GeneID"]
df_concat_GO.to_csv("df_HS_OK____above_antisense_gene_id_GO.csv",sep="\t",index=None)


df_concat.to_csv("df_HS_OK____anti_sense_above.csv",sep="\t",index=None)
df_HS_only.to_csv("df_HS_only_anti_sense_above.csv",sep="\t",index=None)
df_OK_only.to_csv("df_OK_only_anti_sense_above.csv",sep="\t",index=None)


# Subset sizes
s = (
    len(df_HS_only),  # Ab
    len(df_OK_only),  # aB
    len(df_concat),  # AB
)

v = venn2(subsets=s, set_labels=('', ''))

for text in v.set_labels:
    text.set_x(text.get_position()[0] + 0.15)    #Move along x
    text.set_y(text.get_position()[1] + 0.15)    #Move along y
    text.set_fontsize(18)
#label.set_family('serif')
#label.set_x(label.get_position()[0] + 0.1)

#v.get_label_by_id('A').set_text('$x^2$') # Those are set labels
#v.get_label_by_id('A').set_fontsize(22)


hs_label = "Heatshock\n" + str(len(df_HS_only))

# Subset labels
v.get_label_by_id('10').set_text(str(len(df_HS_only)))
v.get_label_by_id('10').set_fontsize(20)
v.get_label_by_id('01').set_text(str(len(df_OK_only)))
v.get_label_by_id('11').set_text(str(len(df_concat)))

#v.get_label_by_id('tdp-1(ok803)').set_fontsize(20)




# Subset colors
v.get_patch_by_id('10').set_color('red')
v.get_patch_by_id('01').set_color('gray')
v.get_patch_by_id('11').set_color('white')

# Subset alphas
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)

# Border styles
c = venn2_circles(subsets=s, linestyle='solid')
c[0].set_ls('dashed')  # Line style
c[0].set_lw(2.0)       # Line width
plt.title("AntisenseRUNON Heatshock and tdp-1(ok803) \n Sense J2 Enrichment compared to WT \n  log2FoldChange > 0 \n (Padj) < 0.05, \n basemean >20", fontsize=20)
plt.savefig('Antisense RUNON Heatshock and tdp-1(ok803) Sense J2 Enrichment compared to WT.png',bbox_inches='tight')
#savefig('foo.png', bbox_inches='tight')
#plt.show()



