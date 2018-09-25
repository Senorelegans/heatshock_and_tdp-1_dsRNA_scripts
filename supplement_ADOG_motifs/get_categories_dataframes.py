import pandas as pd


f1 = "HS_combined_ADOG_Annotated_sequences.csv"
df = pd.read_csv(f1, sep="\t")


dfHS = df[df["Annotated_Condition"]=="HS"]
dfWT = df[df["Annotated_Condition"]=="WT"]
dfBO = df[df["Annotated_Condition"]=="BOTH"]

print(len(dfHS))
print(len(dfWT))
print(len(dfBO))


