import csv
import os, sys


def changing_names(content):
    content[0]=content[0].split("_")[0]

    content[2]=content[2].split("_")[0]

    return content

with open(sys.argv[1], 'r') as data: #input file

    reader = csv.reader(data)
    lines_h, lines_m, lines_y = [], [], []
    for row in reader:
        content = list(row[i] for i in [0,1,2,3,])

        if 'HUMAN' in content[0]:
            #TODO fix return value
            new_content = changing_names(content)
            lines_h.append(new_content)

        if 'MOUSE' in content[0]:
            new_content = changing_names(content)
            lines_m.append(new_content)

        if 'YEAST' in content[0]:
            new_content = changing_names(content)
            lines_y.append(new_content)

with open(sys.argv[2], 'w') as x_file, open(sys.argv[3], 'w') as m_file, open(sys.argv[4], 'w') as y_file: #output files with human, mouse, yeast in file_name
    writer = csv.writer(x_file)
    writer.writerow(["Protein1","Residue1","Protein2", "Residue2"])
    writer.writerows(lines_h)

    writer = csv.writer(m_file)
    writer.writerow(["Protein1","Residue1","Protein2", "Residue2"])
    writer.writerows(lines_m)

    writer = csv.writer(y_file)
    writer.writerow(["Protein1","Residue1","Protein2", "Residue2"])
    writer.writerows(lines_y)
