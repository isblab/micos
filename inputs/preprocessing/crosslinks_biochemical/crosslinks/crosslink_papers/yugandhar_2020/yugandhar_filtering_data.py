import csv

lines = []

with open('yugandhar_supp_table4.csv', 'r') as readFile: #save th einput file as csv otherwise it will give index error

    reader = csv.reader(readFile)
    for row in reader:
        content = list(row[i] for i in [10,11,14,15])
        #changing column values from protein description to just protein names

        if 'MICOS' in content[0] and 'MICOS'in content[2]:
            for j in [0,2]:
                if "MIC60" in content[j]:
                    content[j] = "MIC60"
                if "MIC19" in content[j]:
                    content[j] = "MIC19"

            i = content[1].split('_')[-1] #TODO

            j = content[3].split('_')[-1] #TODO


            content[1]=i

            content[3]=j

            lines.append(content)

with open('../../../mapping_to_homologs/inputs/yu_DSSO.csv', 'w') as writeFile:

    writer = csv.writer(writeFile)
    writer.writerow(["Protein1","Residue1","Protein2", "Residue2"])
    writer.writerows(lines)
