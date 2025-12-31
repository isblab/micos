import pandas as pd
import sys
import matplotlib.colors as mcolors

input_csv = sys.argv[1]
output_csv = sys.argv[2]
threshold = float(sys.argv[3])

df = pd.read_csv(input_csv)

if "Minimum distance" not in df.columns:
    raise ValueError("Input CSV must contain a 'Minimum distance' column.")

grey = mcolors.to_hex("grey").upper()
red = mcolors.to_hex("red").upper()

# Define bins: <= threshold → grey, > threshold → red
bins = [-float('inf'), threshold, float('inf')]
labels = [grey, red]

df['Fixed Colour'] = pd.cut(
    df['Minimum distance'],
    bins=bins,
    labels=labels,
    include_lowest=True
)
df.to_csv(output_csv, index=False)
