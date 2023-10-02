import csv
import pandas as pd

lines = []
with open('Ryl_supp_S2B.csv', 'r') as readFile: #save th einput file as csv otherwise it will give index error

    reader = csv.reader(readFile)
    for row in reader:
        countpsm = list(row[i] for i in [13])
        content = list(row[i] for i in [4,7,6,8])
        #changing column values from protein description to just protein names
        for i in countpsm:
            try:
                i = int(i)
                if i > 1:
                    if 'MICOS' in content[0] and 'MICOS'in content[2]:
                        for j in [0,2]:
                            if "MIC60" in content[j]:
                                content[j] = "MIC60"
                            if "MIC19" in content[j]:
                                content[j] = "MIC19"
                            if "MIC10" in content[j]:
                                content[j] = "MIC10"
                            if "MIC26" in content[j]:
                                content[j] = "MIC26"
                        lines.append(content)

            except:
                ValueError

with open('../../xl_data/scripts/ryl_micos_BS3.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerow(["Protein1","Residue1","Protein2", "Residue2"])
    writer.writerows(lines)
