### to print list of datasets which have micos proteins crosslinks ###
import os
import glob
import pandas as pd

curr_datasets = set()

for file in glob.glob('../xlinkdb-Mic/xlinkdb-Mic*.txt'):
    df = pd.read_table(file)
    unique_datasets = pd.unique(df['Datasets(s)'])
    for entry in unique_datasets:
        curr_datasets.update(entry.split(','))

curr_datasets_list = list(curr_datasets)
print(curr_datasets_list)
