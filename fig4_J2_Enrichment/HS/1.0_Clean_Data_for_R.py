
# coding: utf-8

# In[103]:


import pandas as pd
import os



# In[104]:


os.getcwd()
SCRIPTS = "/Users/M/Google_Drive/Scripts/"
# PATH="/Users/M/Google_Drive/Scripts/hs/worm/EMBO/writeup/Figures/fig4_J2_Enrichment/unique/HS/"
# os.chdir(PATH)



HSPATH = "../../data/Dogcatcher_Out/HS/initial_Rsubread"

df_biotype  = pd.read_csv("../../data/c_elegans.PRJNA13758.WS258.canonical_geneset_ALL_GENES_NO_PATCHES.txt", sep="\t")
df_name = pd.read_csv("../../data/WS258_WB_gene_names_unique.txt", sep="\t", names=["GeneID","name"])



p = os.listdir(HSPATH)
print(p)
f1 = HSPATH + "/DESeq2_sense.csv"
f2 = HSPATH + "/DESeq2_antisense.csv"

LFC = 0



def clean_df(f1):
    """This function will take in a DSeq2 normalized matrix and filter for everything LT 0.05
    Give back two df of above 2lfc and below 2lfc"""
    df = pd.read_csv(f1, index_col=None, sep=',')
    print(f1, "********")
    # df_c = pd.read_csv(count_file, index_col=None, sep='\t')
    # df = pd.merge(df,df_c,how="left",left_index = True,right_index = True,)
    # del df["GeneID"]
    df.rename(columns={"Unnamed: 0" : "GeneID"}, inplace=True)
    print("length : ", len(df))
    # df= df[~df.GeneID.isin(rrna_df.GeneID)] #Get rid of ribosomal reads
    # print("length after rrna removal : ", len(df))
    df = df.dropna()
    del df["pvalue"]
    del df["lfcSE"]
    del df["stat"]
    #if J2 enriched for group, does group have over mean 100
    df = df.sort_values(by='baseMean')
    print("Length of df is", len(df))
    # df = df[df['baseMean'] < 20]
    # dfbelow20 = df[df['baseMean'] < 20]
    # print(dfbelow20)
    # print("amount below 20 basemean is", len(dfbelow20))
    # print("Length of df over 20 basemean is", len(df))
    dfsig = df[df["padj"]<0.05]
    # dfsig.to_csv(f1[:-4]+"_sig.csv")
    print("length after dropNA : ", len(df))
    print("length sig is : ", len(dfsig))
    df_above_2_lfc = dfsig[dfsig['log2FoldChange'] > LFC].sort_values(by='log2FoldChange', ascending=False)
    df_below_2_lfc = dfsig[dfsig['log2FoldChange'] < -LFC].sort_values(by='log2FoldChange', ascending=True)
    print("Total length up padj 0.05 :",  len(df_above_2_lfc))
    print("Total length down padj 0.05 :",len(df_below_2_lfc))
    return df



df = clean_df(f1)
df2 = clean_df(f2)

df = df.merge(df2,on="GeneID")
#print(df.head())

# #Get significant sense and antisense and LFC
#Make significant
df["sense_significant_and_antisense_significant"] = (df["padj_x"] < 0.05) & (df["padj_y"] < 0.05)

df["log2FoldChange_up"] = (df["log2FoldChange_x"] > LFC) & (df["log2FoldChange_y"] > LFC)
df["log2FoldChange_down"] = (df["log2FoldChange_x"] < LFC) & (df["log2FoldChange_y"] < LFC)

df["log2FoldChange_sig_and_sense_significant_and_antisense_significant"] = (df["sense_significant_and_antisense_significant"] == True) & ( (df["log2FoldChange_up"] ==True) | (df["log2FoldChange_down"] ==True) )

# #Merge back with regular files
#Reload dataframes
df_sig = df[[ "GeneID","log2FoldChange_sig_and_sense_significant_and_antisense_significant"] ]


df = clean_df(f1)
df2 = clean_df(f2)

df = df.merge(df_sig,on="GeneID", how="left")
df2 = df2.merge(df_sig,on="GeneID", how="left")

df = df.fillna(False)
df2 = df2.fillna(False)

#put in biotype name


df_biotype["GeneID"] = df_biotype["gene_id_name"]
df_biotype=df_biotype[["GeneID","gene_biotype_name"] ]
# print(df_biotype.head())


df = df.merge(df_biotype, on="GeneID",how="left")
df2 = df2.merge(df_biotype, on="GeneID",how="left")

#Put in regular worm names

df = df.merge(df_name, on="GeneID", how="left")
df2 = df2.merge(df_name, on="GeneID", how="left")
df = df.fillna("NA")
df2 = df2.fillna("NA")

#Write out files for deseq2

df.to_csv("DESeq2_senseunique.csv"[:-4] + "_cleaned.txt", sep="\t", index=None)
df2.to_csv("DESeq2_antisenseunique.csv"[:-4] + "_cleaned.txt", sep="\t", index=None)

df = df[ df['log2FoldChange_sig_and_sense_significant_and_antisense_significant'] == True]
df2 = df2[ df2['log2FoldChange_sig_and_sense_significant_and_antisense_significant'] == True]



df.to_csv("DESeq2_senseunique.csv"[:-4] + "_cleaned_dsRNAonly.txt", sep="\t", index=None)
df2.to_csv("DESeq2_antisenseunique.csv"[:-4] + "_cleaned_dsRNAonly.txt", sep="\t", index=None)




print("DSRNA counts")
df_above_2_lfc = df[df['log2FoldChange'] > LFC].sort_values(by='log2FoldChange', ascending=False)
df_below_2_lfc = df[df['log2FoldChange'] < -LFC].sort_values(by='log2FoldChange', ascending=True)
print("Total length up padj 0.05 :",  len(df_above_2_lfc))
print("Total length down padj 0.05 :",len(df_below_2_lfc))
print("antisense")
# print(df.head())


