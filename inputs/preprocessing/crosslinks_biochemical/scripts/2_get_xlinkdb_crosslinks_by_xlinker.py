### to save ../crosslinks from each datsets based on the crosslinker type ###

# Datasets ChemBiol17AAG_Bruce and Liu2015NatureMethods_Heck are not considered further.
import os, sys
import glob
import pandas as pd

PIR_df = pd.DataFrame()
DSSO_df = pd.DataFrame()
BDP_df = pd.DataFrame()

count = 0
for data in sorted(glob.glob('../crosslinks/xlinkdb/xlinkdb_datasets/*.txt')):
    count += 1
    df1 = pd.read_table(data)
    df2 = df1[(df1['uniprotA'].str.startswith('MIC')) & (df1['uniprotB'].str.startswith('MIC')) & (~df1['uniprotA'].str.contains('MICA'))].copy() #it has all the "MIC" containing rows
    cols_ = df2.filter(['uniprotA', 'modposA','uniprotB','modposB'])
    print(data, len(df2)) # printing the number of rosslinks from each datasets

    if data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/Caudal") or data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/MCF7_") or data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/iqPIR"):
        PIR_df = pd.concat([PIR_df,cols_])

    elif data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/Chavez") or data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/HeLa") or data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/Keller") or data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/Schweppe"):
        BDP_df = pd.concat([BDP_df,cols_])

    elif data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/Liu") or data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/Mango") or data.startswith("../crosslinks/xlinkdb/xlinkdb_datasets/Yeast"):
        DSSO_df = pd.concat([DSSO_df,cols_])

PIR_df = PIR_df.to_csv('../mapping_to_homologs/inputs/PIR.csv', index = False)
BDP_df = BDP_df.to_csv('../mapping_to_homologs/inputs/BDP.csv', index = False)
DSSO_df = DSSO_df.to_csv('../mapping_to_homologs/inputs/DSSO.csv', index = False)
