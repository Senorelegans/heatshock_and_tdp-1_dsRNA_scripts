import pandas as pd


f1 = "HS_combined_ADOG_Annotated_sequences.csv"
df = pd.read_csv(f1, sep="\t")






def writeFasta(df,name_out,length,condition):
    df = df[df["Annotated_Condition"]==condition]
    length = str(length)
    with open(name_out,"w") as outfile:
        for row in range(len(df)):
            outfile.write(">"+df["gene_id_name"].iloc[row] + "\n" + df['fasta_'+length+'_upstream'].iloc[row] + "\n")

length = 500
writeFasta(df,f1[:-4]+"_HS_"+str(length)+".fa",length,"HS")
writeFasta(df,f1[:-4]+"_WT_"+str(length)+".fa",length,"WT")
writeFasta(df,f1[:-4]+"_BO_"+str(length)+".fa",length,"BOTH")




