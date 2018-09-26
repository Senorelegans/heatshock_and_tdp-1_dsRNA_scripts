
# coding: utf-8

# In[107]:

import pandas as pd
import os

import os
cwd = os.getcwd()

# print(os.listdir("../data/Dogcatcher_Out/"))

print(os.listdir("../data/Dogcatcher_Out/OK/N2_vs_OK_2ndtime/FINAL_OUT/ALL"))
# data/Dogcatcher_Out/OK/N2_vs_OK_2ndtime/FINAL_OUT/All/DOG/

# In[117]:


# PATH="/Users/M/Google_Drive/Scripts/hs/worm/EMBO/writeup/Figures/fig7_J2_Enrichment_DOG/unique_DOG"
# # os.chdir(PATH)
multi_or_unique = "unique"



# In[118]:
def clean_df(f1, fout):
    """This function will take in a DSeq2 normalized matrix and filter for everything LT 0.05
    Give back two df of above 2lfc and below 2lfc"""
    df = pd.read_csv(f1, index_col=None, sep='\t')
    # df_c = pd.read_csv(count_file, index_col=None, sep='\t')
    # df = pd.merge(df,df_c,how="left",left_index = True,right_index = True,)
    # del df["GeneID"]
    df.rename(columns={"gene_id_name" : "GeneID"}, inplace=True)
    #df = df.dropna()
    print("length : ", len(df))
    # df= df[~df.GeneID.isin(rrna_df.GeneID)] #Get rid of ribosomal reads
    # print("length after rrna removal : ", len(df))
    # df= df[df["gene_biotype_name"] != "rRNA"] #Get rid of ribosomal reads
    # print("length after rrna removal : ", len(df))
    dfsig = df[df["padj"] < 0.05]
    dfsig.to_csv(fout[:-4]+ "_"+multi_or_unique+"_under_padj05.csv", sep="\t", index=None)
    dfup = dfsig[dfsig["log2FoldChange"] > 0 ]
    print("up", len(dfup))
    dfdo = dfsig[dfsig["log2FoldChange"] < 0 ]
    print("down", len(dfdo))
    dfup.to_csv(fout[:-4]+ "_"+multi_or_unique+"_under_padj05_up.csv")
    dfdo.to_csv(fout[:-4]+ "_"+multi_or_unique+"_under_padj05_do.csv")
    df.to_csv(fout[:-4]+ "_"+multi_or_unique+"_cleaned.csv",index=None)
    return df

HSPATH = "../data/Dogcatcher_Out/HS/N2_vs_HS_2ndtime/FINAL_OUT/ALL"
OKPATH = "../data/Dogcatcher_Out/OK/N2_vs_OK_2ndtime/FINAL_OUT/ALL"


print("HS")
f1 = clean_df(HSPATH + "/DOG/combined_DOG_noOperon.csv", "HS_combined_DOG_noOperon.csv")
print("OK")
f1 = clean_df(OKPATH + "/DOG/combined_DOG_noOperon.csv", "OK_combined_DOG_noOperon.csv")

print("ADOGS.....")

print("HS")
f1 = clean_df(HSPATH + "/ADOG/combined_ADOG_noOperon.csv", "HS_combined_ADOG_noOperon.csv")
print("OK")
f1 = clean_df(OKPATH + "/ADOG/combined_ADOG_noOperon.csv", "OK_combined_ADOG_noOperon.csv")

