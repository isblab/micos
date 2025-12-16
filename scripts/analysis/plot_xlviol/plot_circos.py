import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt
import sys

input_csv = sys.argv[1]
df = pd.read_csv(input_csv)

protein_lengths = {
    'MIC10': 78,
    'MIC13': 118,
    'MIC60': 758,
    'MIC19': 227
}

circos = Circos(sectors=protein_lengths, space=5)

protein_colors = {
    'MIC10': '#d47a00',
    'MIC13': '#2ed766',
    'MIC60': '#1e90ff',
    'MIC19': '#fa8072',
}

for sector in circos.sectors:
    track = sector.add_track((90, 100))
    track.axis(
        fc=protein_colors.get(sector.name, "#CCCCCC"),
        ec="black",
        lw=0.5
    )
    track.text(sector.name, size=20, r=105)

# Links
for _, row in df.iterrows():
    r1 = int(row['Residue1'])
    r2 = int(row['Residue2'])

    circos.link(
        (row['Protein1'], r1, r1 + 1),
        (row['Protein2'], r2, r2 + 1),
        color=row['Fixed Colour'],
        lw=1
    )

circos.plotfig(figsize=(8, 8))
plt.savefig(f'plot_{input_csv.split('.csv')[0]}.png')
# plt.show()
