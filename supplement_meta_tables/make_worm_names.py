import pandas as pd

df = pd.read_excel("S5_file.xlsx",sep="\t", comment='#')

names_list = []
for col in list(df):
    if "names" in col:
        names_list.append(col)


print(names_list)
# print(df.head())

for x in range(len(df)):
    for col in names_list:
        df_list = pd.DataFrame()
        l = df[col].iloc[x]


