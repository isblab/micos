## to map both xlinks or pairwise biochemical restraints
## write input file, output file, query species name and target species name as sys.argv

import Bio
from Bio.SeqIO import parse
from Bio import AlignIO
import os,sys
import pandas as pd

# function to get query and target species sequences ( MSA aligned ) as lists
def map_to_MSA(protein):

    protein_query, protein_target = [],[]

    if protein == 'MIC10':
        protein_query = mapping_MSA["MIC10_query"]
        protein_target = mapping_MSA["MIC10_target"]

    elif protein == 'MIC19':
        protein_query = mapping_MSA["MIC19_query"]
        protein_target = mapping_MSA["MIC19_target"]

    elif protein == 'MIC13':
        protein_query = mapping_MSA["MIC13_query"]
        protein_target = mapping_MSA["MIC13_target"]

    elif protein == 'MIC60':
        protein_query = mapping_MSA["MIC60_query"]
        protein_target = mapping_MSA["MIC60_target"]

    elif protein == 'MIC25':
        protein_query = mapping_MSA["MIC25_query"]
        protein_target = mapping_MSA["MIC25_target"]

    elif protein == 'MIC26':
        protein_query = mapping_MSA["MIC26_query"]
        protein_target = mapping_MSA["MIC26_target"]

    elif protein == 'MIC27':
        protein_query = mapping_MSA["MIC27_query"]
        protein_target = mapping_MSA["MIC27_target"]

    return [protein_query, protein_target]


def getting_residue_range(query,target,start,end,p1_name): # for mapping residue ranges; for crosslinks start and end are same
    start_res, end_res = 0,0

    if start == 0 and end == 0: # otherwise it will try to align start_res = 0 to a position in MSA
        return (p1_name,start,end) # so we will print it as it is


    for res in range(len(query)):
        if query[res] == start: # checking that residue number in the query sequence list
            start_res = res #remember this is index, actual position of start_res in query list is res+1

        if query[res] == end:
            end_res = res

    output = []

    for id_ in range(start_res,end_res+1):
        if target[id_] != 0:
            output.append(target[id_])

    if len(output) == 0: # to fix value error; if the length of target sequence is less than that of query, don't print it in the output file
        return (p1_name,0)
    else:
        if min(output) == max(output): # this is for crosslinks
            return (p1_name,output[0])
        else:
            return (p1_name,min(output),max(output)) #this is for mapping range

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



# here is the first function which would be read by python
def file_parsing(data):
    total_restraints = pd.DataFrame()

    if len(data.columns) == 4: #this shows it's a crosslink file, with query species name written in the name of the file and we are giving it as a sys.argv
        print("this is xlink file")
        final = []

        for i in range(len(data)):
            new = ()
            for j in [0,2]: # reading one protein and its residue at a time
              protein1 = data.iloc[i, j]
              residue1 = data.iloc[i, j+1]
              # if protein1 != 'MIC27' and protein1 != 'MIC26': # we are not interested in these proteins
              list_seq = map_to_MSA(protein1) # going to this function to get two aligned seq as lists: query and target
              new += getting_residue_range(list_seq[0],list_seq[1],residue1,residue1,protein1) # they are needed here to map the residues; note: residue1 is parsed twice
            # new_tuple = new + (data.iloc[i,4],) #adding the xlinktype entry


            if len(new) == 4: # to make sure both proteins are "our protein of interest"
                # print(new)
                if new[1] != 0 and new[3] != 0: # to exclude xlinks mapped to a gap (in which case it will print 0)
                    final.append(new)
                    total_restraints = pd.DataFrame(final, columns = ['Protein1', 'Residue1', 'Protein2', 'Residue2'])

    elif len(data.columns) == 7: #this shows pairwise file or single protein but its written twice in the file
        print("this is pariwise file")
        final = []

        for i in range(len(data)) :
            new1 = (data.iloc[i,0],) # checking the species name

            if new1[0] == target_species: # for the restraint in target_species only, write that as it is to the output file
                final.append(tuple(data.iloc[i,:],))

            else:
                new = (target_species,)
                for j in [1,4]:
                    protein1 = data.iloc[i, j]
                    residue_start = data.iloc[i, j+1]
                    residue_end = data.iloc[i, j+2]
                    # if protein1 != 'MIC27' and protein1 != 'MIC26':
                    list_seq = map_to_MSA(protein1)
                    new += getting_residue_range(list_seq[0],list_seq[1],residue_start,residue_end,protein1)

                final.append(new)
        total_restraints = pd.DataFrame(final, columns = ['species','Protein1', 'start','end','Protein2','start','end'])


    return total_restraints


input_file = sys.argv[1]
output_file = sys.argv[2]
query_species = sys.argv[3]
target_species = sys.argv[4]

# MSA input file from MAFFT
MIC10_MSA = "/home/muskaan/Documents/modeling_micos_complex/supp_data/MSA/MAFFT/Mic10_mafft_fasta.txt"
MIC13_MSA = "/home/muskaan/Documents/modeling_micos_complex/supp_data/MSA/MAFFT/Mic13_mafft_fasta.txt"
MIC60_MSA = "/home/muskaan/Documents/modeling_micos_complex/supp_data/MSA/MAFFT/Mic60_mafft_fasta.txt"
MIC19_MSA = "/home/muskaan/Documents/modeling_micos_complex/supp_data/MSA/MAFFT/Mic19_mafft_fasta.txt"
MIC25_MSA = "/home/muskaan/Documents/modeling_micos_complex/supp_data/MSA/MAFFT/Mic25_mafft_fasta.txt"
MIC26_MSA = "/home/muskaan/Documents/modeling_micos_complex/supp_data/MSA/MAFFT/Mic26_mafft_fasta.txt"
MIC27_MSA = "/home/muskaan/Documents/modeling_micos_complex/supp_data/MSA/MAFFT/Mic27_mafft_fasta.txt"


#dictionary for getting sequence as list for all the proteins in one go
mapping_MSA = {"MIC10_query" :  converting_to_list(query_species, MIC10_MSA),"MIC10_target" : converting_to_list(target_species, MIC10_MSA),
               "MIC13_query" : converting_to_list(query_species, MIC13_MSA),"MIC13_target" : converting_to_list(target_species, MIC13_MSA),
               "MIC19_query" : converting_to_list(query_species, MIC19_MSA),"MIC19_target" : converting_to_list(target_species, MIC19_MSA),
               "MIC25_query" : converting_to_list(query_species, MIC25_MSA),"MIC25_target" : converting_to_list(target_species, MIC25_MSA),
               "MIC60_query" : converting_to_list(query_species, MIC60_MSA),"MIC60_target" : converting_to_list(target_species, MIC60_MSA),
               "MIC26_query" : converting_to_list(query_species, MIC26_MSA),"MIC26_target" : converting_to_list(target_species, MIC26_MSA),
               "MIC27_query" : converting_to_list(query_species, MIC27_MSA),"MIC27_target" : converting_to_list(target_species, MIC27_MSA)}

# print(mapping_MSA["MIC10_query"])
data = pd.read_csv(input_file, index_col = None)

df = file_parsing(data).to_csv(output_file, index = False)
