import os, sys
import pandas as pd

intracrosslink = pd.DataFrame()
intercrosslink = pd.DataFrame()

input_file = sys.argv[1]
output_file = sys.argv[2]

data = pd.read_csv(input_file)

intracrosslink = data[data['Protein1'] == data['Protein2']]
intercrosslink = data[data['Protein1'] != data['Protein2']]

intracrosslink.to_csv(f"intracrosslink_{output_file}.csv", index=False)
intercrosslink.to_csv(f"intercrosslink_{output_file}.csv", index=False)
