import pandas as pd
import sys

supp_file = sys.argv[1]
df = pd.read_excel(supp_file, sheet_name= 'Supp Data 1_Combined (DSSO)')

micos_prot = ['IMMT', 'MINOS1', 'CHCHD3', 'QIL1']
xl = []

for prot in micos_prot:
    series_a = df.loc[df['gene_a'].astype(str).str.contains(prot, case=False, na=False), 'crosslinks_ab']
    xl.extend(series_a.astype(str).tolist())  # convert to list of strings and append

    series_b = df.loc[df['gene_b'].astype(str).str.contains(prot, case=False, na=False), 'crosslinks_ba']
    xl.extend(series_b.astype(str).tolist())


micos_list = []
nonmicos_list = []

for s in xl:
    data = s.split('#')
    df = pd.DataFrame([x.split('-') for x in data], columns = ['Protein1', 'Residue1', 'Protein2', 'Residue2'])
    df_micos = df[(df['Protein1'].isin(micos_prot)) & (df['Protein2'].isin(micos_prot))]
    df_nonmicos = df[~(df['Protein1'].isin(micos_prot)) | ~(df['Protein2'].isin(micos_prot))]

    micos_list.append(df_micos)
    nonmicos_list.append(df_nonmicos)


micos_xlinks = pd.concat(micos_list, ignore_index=True)
nonmicos_xlinks = pd.concat(nonmicos_list, ignore_index=True)

micos_xlinks = micos_xlinks.drop_duplicates()
nonmicos_xlinks = nonmicos_xlinks.drop_duplicates()
micos_xlinks.reset_index(inplace=True)
nonmicos_xlinks.reset_index(inplace=True)

## Remove all mic60 crosslinks in N and TM domains (residues less than 410)
to_drop = []
for i in micos_xlinks.index:
    if micos_xlinks.iloc[i, 1] == 'IMMT' and micos_xlinks.iloc[i, 3] == 'IMMT':
        if (int(micos_xlinks.iloc[i, 2]) < 410) or (int(micos_xlinks.iloc[i, 4]) < 410):
            to_drop.append(i)

    ## Rename proteins
    if micos_xlinks.iloc[i, 1] == 'IMMT' or micos_xlinks.iloc[i, 3] == 'IMMT':
        micos_xlinks.iloc[i, 1] = 'MIC60'
        micos_xlinks.iloc[i, 3] = 'MIC60'
    elif micos_xlinks.iloc[i, 1] == 'CHCHD3' or micos_xlinks.iloc[i, 3] == 'CHCHD3':
        micos_xlinks.iloc[i, 1] = 'MIC19'
        micos_xlinks.iloc[i, 3] = 'MIC19'
    elif micos_xlinks.iloc[i, 1] == 'MINOS1' or micos_xlinks.iloc[i, 3] == 'MINOS1':
        micos_xlinks.iloc[i, 1] = 'MIC10'
        micos_xlinks.iloc[i, 3] = 'MIC10'

micos_xlinks = micos_xlinks.drop(to_drop).reset_index(drop=True)


# Save the files as excel sheet
micos_xlinks.to_excel('micos_xlinks.xlsx', index=False)
nonmicos_xlinks.to_excel('nonmicos_xlinks.xlsx', index=False)