import pandas as pd
import os, sys

media =  sys.argv[1]
df1 = pd.read_csv(media, index_col = False)

micos_only = df1[(df1['Protein 1'].str.startswith('MIC')) & (df1['Protein 2'].str.startswith('MIC'))].copy() # taking only the micos proteins

out_df = pd.DataFrame(columns = ['Protein1','Residue1', 'Protein2','Residue2']) #output DataFrame

for row in range(len(micos_only)):

    out_list = [] #new list for each row
    out_list.append(micos_only['Protein 1'].str.split('_').iloc[row][0])
    out_list.append(micos_only['Residue1'].iloc[row])
    out_list.append(micos_only['Protein 2'].str.split('_').iloc[row][0])
    out_list.append(micos_only['Residue2'].iloc[row])

    out_df.loc[len(out_df)] = out_list # appending list as row

out_df.to_csv(sys.argv[2], index = False)
