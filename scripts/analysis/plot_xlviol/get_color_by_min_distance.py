import pandas as pd
import sys
import matplotlib.colors as mcolors

input_csv = sys.argv[1]
output_csv = sys.argv[2]

df = pd.read_csv(input_csv)

if "Minimum distance" not in df.columns:
    raise ValueError("Input CSV must contain a 'Distance' column.")

bins = [0, 10, 20, 30, 40, 50, 60, float('inf')]
colors = ["red", "blue","orange", "grey", "yellow", "purple", "green"]
hex_colors = [mcolors.to_hex(c).upper() for c in colors]

# Create a new column 'Fixed Colour' based on distance
df['Fixed Colour'] = pd.cut(df['Minimum distance'], bins=bins, labels=hex_colors, include_lowest=True)

# Save to CSV
df.to_csv(output_csv, index=False)
