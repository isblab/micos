import pandas as pd
import os, sys

input_file = sys.argv[1]
protein1_start = int(sys.argv[2]) # according to name of csv file, protein1 and protein2
protein2_start = int(sys.argv[3])

df = pd.read_csv(input_file, sep=' ', header=None)
header_values = [str(i) for i in range(protein2_start, df.shape[1] + 1)]
df.columns = header_values
df.index = range(protein1_start, protein1_start + len(df))

# Save values which are equal to 1 along with headers and index to a csv file
output_csv_values_file = f'{os.path.splitext(os.path.basename(input_file))[0]}_contacts_20A.csv'
df_values_equal_to_1 = df[df == 1].dropna(how='all').dropna(axis=1, how='all')
df_values_equal_to_1.to_csv(output_csv_values_file, index=True, sep=',')
