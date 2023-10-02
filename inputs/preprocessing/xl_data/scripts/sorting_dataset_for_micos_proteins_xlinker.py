### to save crosslinks from each datsets based on the crosslinker type ###
import os, sys
import glob
import pandas as pd

PIR_df = pd.DataFrame()
DSSO_df = pd.DataFrame()
BDP_df = pd.DataFrame()

count = 0
for data in sorted(glob.glob('../xlinkdb-datasets/*.txt')):
    count += 1
    df1 = pd.read_table(data)
    df2 = df1[(df1['uniprotA'].str.startswith('MIC')) & (df1['uniprotB'].str.startswith('MIC')) & (~df1['uniprotA'].str.contains('MICA'))].copy() #it has all the "MIC" containing rows
    cols_ = df2.filter(['uniprotA', 'modposA','uniprotB','modposB'])
    print(data, len(df2)) # printing the number of crosslinks from each datasets

    if data.startswith("../xlinkdb-datasets/Caudal") or data.startswith("../xlinkdb-datasets/MCF7_") or data.startswith("../xlinkdb-datasets/iqPIR"):
        PIR_df = pd.concat([PIR_df,cols_])

    elif data.startswith("../xlinkdb-datasets/Chavez") or data.startswith("../xlinkdb-datasets/HeLa") or data.startswith("../xlinkdb-datasets/Keller") or data.startswith("../xlinkdb-datasets/Schweppe"):
        BDP_df = pd.concat([BDP_df,cols_])

    elif data.startswith("../xlinkdb-datasets/Liu") or data.startswith("../xlinkdb-datasets/Mango") or data.startswith("../xlinkdb-datasets/Yeast"):
        DSSO_df = pd.concat([DSSO_df,cols_])

PIR_df = PIR_df.to_csv('../inputs/PIR_xl.csv', index = False)
BDP_df = BDP_df.to_csv('../inputs/BDP_xl.csv', index = False)
DSSO_df = DSSO_df.to_csv('../inputs/DSSO_xl.csv', index = False)
