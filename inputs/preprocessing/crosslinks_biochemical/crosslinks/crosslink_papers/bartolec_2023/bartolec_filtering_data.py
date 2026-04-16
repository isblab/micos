import pandas as pd

medias =  ['DHSO', 'DSSO']

for media in medias:
    df1 = pd.read_csv(f'bartolec_xlinks_{media}.csv', index_col = False)
    print(df1)
    micos_only = df1[(df1['Protein.Descriptions.A'].str.startswith('MIC')) & (df1['Protein.Descriptions.B'].str.startswith('MIC'))].copy() # taking only the micos proteins
    out_df = pd.DataFrame(columns = ['Protein1','Residue1', 'Protein2','Residue2']) #output DataFrame

    for row in range(len(micos_only)):

        out_list = [] #new list for each row
        out_list.append(micos_only['Protein.Descriptions.A'].str.split('_').iloc[row][0])
        out_list.append(micos_only['Leading.Protein.Position.A'].iloc[row])
        out_list.append(micos_only['Protein.Descriptions.B'].str.split('_').iloc[row][0])
        out_list.append(micos_only['Leading.Protein.Position.B'].iloc[row])

        out_df.loc[len(out_df)] = out_list # appending list as row

    out_df.to_csv(f'../../../mapping_to_homologs/inputs/bartolec_{media}.csv', index = False)
