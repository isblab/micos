## to map both xlinks or pairwise biochemical restraints
## write input file, output file, query species name and target species name as sys.argv

import Bio
from Bio.SeqIO import parse
from Bio import AlignIO
import os,sys
import csv

# MSA input file from MAFFT
MIC10 = "Mic10_mafft_fasta.txt"
MIC13 = "Mic13_mafft_fasta.txt"
MIC60 = "Mic60_mafft_fasta.txt"
MIC19 = "Mic19_mafft_fasta.txt"

def getting_residue_range(query,target,start,end,p1_name,p2_name): # for mapping residue ranges (pariwise biochemical restraint)
    start_res, end_res = 0,0

    for res in range(len(query)):
        if query[res] == start: # checking that residue number in the query sequence list
            start_res = res #remember this is index, actual position of start_res in query list is res+1

        if query[res] == end:
            end_res = res

    output = []

    for id_ in range(start_res,end_res+1):
        if target[id_] != 0:
            output.append(target[id_])

    if len(output) != 0: # to fix value error; if the length of target sequence is less than that of query, don't print it in the output file
        return (sys.argv[4], p1_name,min(output),max(output),sys.argv[4],p2_name)

def mapping_to_human(col1,p1,p2):
    p1_name = p1
    p2_name = p2
    if p1 == 'MIC10':
        p1 = MIC10
    elif p1 == 'MIC13':
        p1 = MIC13
    elif p1 == 'MIC60':
        p1 = MIC60
    elif p1 == 'MIC19' or p1 == 'MIC25':
        p1 = MIC19

    if p2 == 'MIC10':
        p2 = MIC10
    elif p2 == 'MIC13':
        p2 = MIC13
    elif p2 == 'MIC60':
        p2 = MIC60
    elif p2 == 'MIC19' or p2 == 'MIC25':
        p2 = MIC19

    if len(data_list[0]) == 8: #for pairwise restraint, we are checking only one protein
        query = converting_to_list(col1,p1) # the species name should be present in the fasta id
        target = converting_to_list(sys.argv[4],p1)
        return getting_residue_range(query,target,start,end,p1_name,p2_name)

    else:
        if p1 == p2: #here we are checking two proteins and their residues (xlink restraints) (intra-crosslinks)
            query = converting_to_list(col1,p1) # the species name should be present in the fasta id, col1 is the query species name
            target = converting_to_list(sys.argv[4],p1)

        else:
            for i in [p1,p2]: # inter-crosslinks
                query = converting_to_list(col1,i) # the species name should be present in the fasta id
                target = converting_to_list(sys.argv[4],i)
        return getting_residues(query, target,residue1,residue2,p1_name,p2_name)

def getting_residues(query,target,start_range,end_range,p1_name,p2_name):
    start_res, end_res = 0,0
    for res in range(len(query)):

        if query[res] == start_range: # checking that residue number in the query sequence list
            start_res = res #remember this is index, actual position of start_res in query list is res+1

        if query[res] == end_range:
            end_res = res

    output = []

    if target[start_res] != 0:
        if target[end_res] != 0:
            output.append(p1_name)
            output.append(target[start_res])
            output.append(p2_name)
            output.append(target[end_res])

    return (output) #printing the output in p1,r1,p2,r2 form


def converting_to_list(species_name, protein_name): #converting the query and target sequences to list
    _list = []
    count = 0
    i = 0

    for record in parse(protein_name, "fasta"):
        if species_name in record.id:
            while i < len(record.seq):
                if record.seq[i] == '-':
                    _list.append(0) # adding 0 for gaps
                elif record.seq[i] != '-':
                    count += 1
                    _list.append(count)
                i += 1
    return _list

final = []

with open (sys.argv[1], "r") as data_file:
      data_file = csv.reader(data_file)
      data_list = list(data_file)

      if "MIC" in data_list[0][0]:
          print("this is xlink file") #to print which type of file is there in input
          col1 = sys.argv[3]
          for i in range(len(data_list)):
              if str(data_list[i][0]) != 'MIC27' and str(data_list[i][0]) != 'MIC26' and str(data_list[i][2]) != 'MIC27'and str(data_list[i][2]) != 'MIC26':
                  protein1 = str(data_list[i][0])
                  protein2 = str(data_list[i][2])
                  residue1 = int(data_list[i][1])
                  residue2 = int(data_list[i][3])

                  corresponding_res_in_human = mapping_to_human(col1,protein1,protein2)
                  final.append(corresponding_res_in_human)



      elif len(data_list[0]) == 8: # because its a
          print("this is Pairwise_table")
          for i in range(len(data_list)):
              if str(data_list[i][0]) == sys.argv[4]:
                  final.append(data_list[i])

              else:
                  col1 = sys.argv[3]
                  if str(data_list[i][1]) != 'MIC26' and str(data_list[i][5]) != 'MIC26':
                      protein1 = str(data_list[i][1])
                      protein2 = str(data_list[i][5])
                      start = int(data_list[i][2])
                      if data_list[i][3] == '': # if only one residue is known to interact in biochemical restraints (not a range)
                          end = 0
                      else:
                          end = int(data_list[i][3])

                      corresponding_res_in_human = mapping_to_human(col1,protein1,protein2)
                      if corresponding_res_in_human is not None: 
                          final.append(corresponding_res_in_human)

with open(sys.argv[2], 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerows(final)
