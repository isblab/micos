import pandas as pd
import glob
# file name format is Crosslinker_reference study_source(human).csv
out_df = pd.DataFrame(columns = ['Protein1','Residue1', 'Protein2','Residue2','Crosslinker type','Reference']) #output DataFrame

for files in glob.glob('*.csv'):
    df = pd.read_csv(files)
    components = files.split("_")
    Crosslinker = components[0]
    Reference = components[1]

    for index, row in df.iterrows():
        out_list = [] #new list for each row
        out_list.append(row[0])
        out_list.append(row[1])
        out_list.append(row[2])
        out_list.append(row[3])
        out_list.append(Crosslinker)
        out_list.append(Reference)

        out_df.loc[len(out_df)] = out_list # appending list as row

out_df.to_csv('combined.csv', index = False)
