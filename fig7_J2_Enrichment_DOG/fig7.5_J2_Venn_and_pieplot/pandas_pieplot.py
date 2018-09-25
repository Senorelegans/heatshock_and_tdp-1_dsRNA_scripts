#!/usr/bin/env python
import argparse
import csv
import numpy as np
import os
import sys
import pandas as pd
import string
import re
import glob

import pandas as pd
import matplotlib.pyplot as plt
#matplotlib.style.use('ggplot')




# Running in parallel

##################
#### PART 2 TAKE out run_on_overlap_opposite_genes from longest contig


# for file in glob.glob("Runon/less_than_padj_0.05/*/*nooppositestrand.csv"):
#     print(file)
#     df = pd.read_csv(file, sep="\t")
#     df.to_csv(file[:-4]+"_nooppositestrand.csv",sep="\t",index=None)



file = "df_HS_sense_above_padj.05.csv"

df = pd.read_csv(file, sep="\t")
x = str(len(df))
count_df = df.groupby("gene_biotype_name").count()
count_df2 = count_df[["GeneID"]].copy()
#count_df2.plot(kind="pie", figsize=(15,15), autopct='%.2f', fontsize=20,legend = False, title=x + " enriched sense run on in Heatshock vs Wildtype"   ,use_index=False, subplots=True)#, colormap="Pastel1")

fig = count_df2.plot(kind="pie", figsize=(15,15), autopct='%.2f', fontsize=25,legend = False,use_index=False, subplots=True)#, colormap="Pastel1")
plt.title("Biotypes of runon sections \n sense J2 enriched in heatshock", fontsize=30)
plt.axis('off')
plt.savefig(file[:-4] + "_PIE.png")


file = "df_HS_sense_below_padj.05.csv"
df = pd.read_csv(file, sep="\t")
x = str(len(df))
count_df = df.groupby("gene_biotype_name").count()
count_df2 = count_df[["GeneID"]].copy()
count_df2.plot(kind="pie", figsize=(15,15), autopct='%.2f', fontsize=25,legend = False,use_index=False, subplots=True)
plt.title("Biotypes of runon sections \n antisense J2 enriched in heatshock", fontsize=30)
#ax = count_df2.plot(kind="pie", figsize=(15,15), autopct='%.2f', fontsize=10,legend = True, title=x + " enriched antisense run on in Heatshock vs Wildtype"   ,use_index=False, subplots=True)#,
# colormap="Pastel1")
plt.axis('off')
plt.savefig(file[:-4] + "_PIE.png")
