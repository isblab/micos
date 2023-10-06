import pandas as pd
import os, sys

media =  sys.argv[1]
df1 = pd.read_csv(media, index_col = False)

micos_only = df1[(df1['Protein.Descriptions.A'].str.startswith('MIC')) & (df1['Protein.Descriptions.B'].str.startswith('MIC'))].copy() # taking only the micos proteins
out_df = pd.DataFrame(columns = ['Protein1','Residue1', 'Protein2','Residue2']) #output DataFrame

for row in range(len(micos_only)):

    out_list = [] #new list for each row
    out_list.append(micos_only['Protein.Descriptions.A'].str.split('_').iloc[row][0])
    out_list.append(micos_only['Leading.Protein.Position.A'].iloc[row])
    out_list.append(micos_only['Protein.Descriptions.B'].str.split('_').iloc[row][0])
    out_list.append(micos_only['Leading.Protein.Position.B'].iloc[row])

    out_df.loc[len(out_df)] = out_list # appending list as row

out_df.to_csv(sys.argv[2], index = False)
