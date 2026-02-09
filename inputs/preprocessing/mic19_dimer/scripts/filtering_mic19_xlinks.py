import pandas as pd
import glob
import os, sys
self_links = pd.DataFrame()
monomeric_links = pd.DataFrame()
input_file = sys.argv[1]
output_file = sys.argv[2]

data = pd.read_csv(input_file)
df_filtered = data[(data['Protein1'].str.startswith('MIC19')) & (data['Protein2'].str.startswith('MIC19'))].copy()

low, high = 60, 175

df_filtered['Residue1'] = pd.to_numeric(df_filtered['Residue1'], errors='coerce')
df_filtered['Residue2'] = pd.to_numeric(df_filtered['Residue2'], errors='coerce')

same_residue = df_filtered['Residue1'] == df_filtered['Residue2']

within_range = (
    df_filtered['Residue1'].between(low, high) &
    df_filtered['Residue2'].between(low, high)
)

self_links = df_filtered[same_residue & within_range]
monomeric_links = df_filtered[~same_residue & within_range]

# same_residue = df_filtered['Residue1'] == df_filtered['Residue2']

# self_links = pd.concat([self_links, df_filtered[same_residue]])
# monomeric_links = pd.concat([monomeric_links, df_filtered[~same_residue]])

print(input_file, 
      "Total mic19 cc:", (len(df_filtered[same_residue & within_range]) + len(df_filtered[~same_residue & within_range])),
      "Total:", len(df_filtered), 
      "Self-links:", len(df_filtered[same_residue & within_range]), 
      "monomeric_links:", len(df_filtered[~same_residue & within_range]))

self_links.to_csv(f"../self_links/self_{output_file}.csv", index=False)
monomeric_links.to_csv(f"../monomeric_links/monomeric_{output_file}.csv", index=False)
