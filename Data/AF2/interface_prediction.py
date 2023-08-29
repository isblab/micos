import numpy as np
import sys
import pickle
from Bio.PDB import PDBParser
import json
import pandas as pd

#TODO 5 A, 10 A
#  
path = sys.argv[1] # path to output of AF2 multimer
outf = sys.argv[2] # path to where output should be stored
name = sys.argv[3] # name of complex e.g. dp-pg


pdb = f'{path}/ranked_0.pdb'
with open(f'{path}/ranking_debug.json', 'r') as f:
    o = json.loads(f.read())
best_pkl = o['order'][0]
pkl = f'{path}/result_{best_pkl}.pkl'
with open(pkl, 'rb') as f:
    data = pickle.load(f)

models = PDBParser().get_structure('pdb', pdb)

f2 = open(f'{outf}/{name}_distance_cutoff.txt', 'w')
# f3 = open(f'{outf}/{name}_interface_pae.txt', 'w')

all_pairs = []
interface_pairs = []

## defining functions##

def custom_sort(item):
    return (item[0], item[2], -item[4])


def calculate_distance (res1, res2, res3, res4):

    for model in models: #TODO do not read models again and again. 

        chains = [c for c in model.get_chains()]

        if res1 in chains[0] and res2 in chains[0]:
            dist1 = abs(chains[0][res1]['CA'] - chains[0][res2]['CA'])

        if res3 in chains[1] and res4 in chains[1]:
            dist2 = abs(chains[1][res3]['CA'] - chains[1][res4]['CA'])


    return create_patch(res1, res2, res3, res4, dist1, dist2)


def create_patch(res1, res2, res3, res4, dist1, dist2):

     # '''take list of res_pairs
     #    check if the distance between two pairs are less than 5
     #    if yes, save in one patch
     #    if not, create and save in new patch
     #    do it for all res_pairs going through all patches, if none of them matches, then create new patch
     #    return list of all created patches'''

    lists = []
    new_list = []

    if dist1 < 10.0 and dist2 < 10.0:
        lists.append((res1, res3, res2, res4))

    else:
        new_list.append((res1, res3))
        new_list.append((res2, res4))

    return lists, new_list



## to calculate PAE in a given patch
# def calculate_pae(patch):
#     # '''for all res pairs in input (which is a patch), calculate pae all vs all with two for loops
#     # increment the values in a variable
#     # take the average and that will be the average pae for the patch
#     # return that and use it for cutoff '''
#
#     var = 0
#     count = 0
#     confident_patches = []
#
#     for i in range(len(patch)):
#         for j in range(1, len(patch)):
#             var += data['predicted_aligned_error'][patch[i][0]][patch[j][1]]
#             count += 1
#
#     avg_pae = var / count
#
#     for pair in patch:
#         if data['predicted_aligned_error'][pair[0]][pair[1]] < avg_pae:
#             confident_patches.append((pair[0], pair[1], data['predicted_aligned_error'][pair[0]][pair[1]]))
#
#     return confident_patches


## reading the pdb structure

for model in models:

    chains = [c for c in model.get_chains()]

    len_chain_1 = len([r for r in chains[0]]) #TODO remove if not used 

    for i,resa in enumerate(chains[0]):

        sorted_pairs = 0
        all_pairs = [] {(i+1,resa,j+1,resb):atoma-atomb} 2. [i+1,resa]
        for atoma in resa:
            if atoma.get_bfactor() > 70:
                for j, resb in enumerate(chains[1]):
                    for atomb in resb:
                        if atomb.get_bfactor() > 70:
                            if atoma-atomb < 10.0:
                                if not (i+1, atoma, j+1, atomb) in all_pairs:
                                all_pairs.append((i+1, atoma, j+1, atomb, atoma-atomb))
                                #TODO chains[0][i]; resa.get_resnum(), resa.get_com() 
                                #TODO store chains A, B. resa.get_resnum(), A, resa.get_com(),
                                #TODO need a break/continue 
                                #TODO linalg.norm()

        # to sort the interface residue pairs and selecting the residue pair with the shortest distance between atoms
        sorted_pairs = sorted(all_pairs, key = custom_sort)
        if len(sorted_pairs) != 0:
            for pair in range(len(sorted_pairs)):
                try:
                    if sorted_pairs[pair][0] == sorted_pairs[pair+1][0] and sorted_pairs[pair][2] != sorted_pairs[pair+1][2]:
                        interface_pairs.append(sorted_pairs[pair])
                    if sorted_pairs[pair][0] != sorted_pairs[pair+1][0]:
                        interface_pairs.append(sorted_pairs[pair])

                except:
                    IndexError
                    interface_pairs.append(sorted_pairs[pair])

# calculating the distance between adjacent residues
#TODO look at all previous patches. 
for i in range(len(interface_pairs)-1):

    if (interface_pairs[i][0]-interface_pairs[i+1][0]) < 10 A and ([2]):
        

    dist = calculate_distance(interface_pairs[i][0], interface_pairs[i+1][0], interface_pairs[i][2], interface_pairs[i+1][2])
    # print(interface_pairs[i][0], interface_pairs[i+1][0], interface_pairs[i][2], interface_pairs[i+1][2], dist)


    f2.write(f'{interface_pairs[i][0], interface_pairs[i+1][0], interface_pairs[i][2], interface_pairs[i+1][2], dist}\n') #TODO do not write to file till the end 

print(dist)

#TODO implement all vs all PAE and interface pairs? 
#TODO implement 5 A or 10 A 

# file format is
