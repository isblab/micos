'''
Given an CSV file for XLs in Protein1,Residue1,Protein2,Residue2 format,
this script will split it in the given ratio for use in Nested Sampling
'''
import os,sys
import random
import pandas as pd

xl_file = sys.argv[1]
perc_to_evi = 0.1

xls = []
# header = None
with open(xl_file,'r') as xlf:
    for ln in xlf.readlines():
        if (not ln.startswith('Protein1')) and (not ln.startswith('Linker')):
            xls.append(ln)
        else:
            header = ln

sampling, evi_calc = [], []
for link in xls:
    rng = random.random()
    if rng<perc_to_evi:
        evi_calc.append(link)
    else:
        sampling.append(link)

fname = xl_file.split('/')[-1]
dir_path = xl_file.split('/')
if len(dir_path)>1:
    dir_path = '/'.join(dir_path[0:-1])
else:
    dir_path = '.'
with open(f'{dir_path}/sampling_{fname}','w') as sf:
    if not header is None:
        sf.write(header)
    for lnk in sampling:
        sf.write(lnk)

with open(f'{dir_path}/evicalc_{fname}','w') as evif:
    if not header is None:
        evif.write(header)
    for lnk in evi_calc:
        evif.write(lnk)

####
# intracrosslink = pd.DataFrame()
# intercrosslink = pd.DataFrame()
#
# input_file = sys.argv[1]
# output_file = sys.argv[2]
#
# data = pd.read_csv(input_file)
#
# intracrosslink = data[data['Protein1'] == data['Protein2']]
# intercrosslink = data[data['Protein1'] != data['Protein2']]
#
# intracrosslink.to_csv(f"intracrosslink_{output_file}.csv", index=False)
# intercrosslink.to_csv(f"intercrosslink_{output_file}.csv", index=False)
