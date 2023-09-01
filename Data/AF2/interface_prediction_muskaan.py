import numpy as np
import sys
import pickle
from Bio.PDB import PDBParser
import json
from itertools import product
from sklearn.cluster import AffinityPropagation
from scipy.spatial import KDTree
import argparse

parser = argparse.ArgumentParser(description="predciting interface residues")
parser.add_argument("--interface_cutoff","-i", help="defining interface distance cutoff", default=10.0, type=float)
parser.add_argument("--intra_chain_residue_distance_cutoff","-d", help="defining distance between two residue pairs cutoff", default=10.0, type=float)
parser.add_argument("--path","-p", help="path to output of AF2 multimer", default=None)
parser.add_argument("--output_path","-o", help="output path", default=None)
parser.add_argument("--out_name","-n", help="output file name", default=None)
# parser.add_argument("--method","-m", help="which method", default='algorithm')

args = parser.parse_args()

path = args.path # path to output of AF2 multimer
outf = args.output_path # path to where output should be stored
name = args.out_name # name of complex e.g. dp-pg
interface_cutoff = args.interface_cutoff
intra_chain_residue_distance_cutoff = args.intra_chain_residue_distance_cutoff
# method = args.method

# Get the AF2 metrics for the predicted complex.
pdb = f'{path}/ranked_0.pdb'
with open(f'{path}/ranking_debug.json', 'r') as f:
    o = json.loads(f.read())
best_pkl = o['order'][0]
pkl = f'{path}/result_{best_pkl}.pkl'
with open(pkl, 'rb') as f:
    data = pickle.load(f)

models = PDBParser().get_structure('pdb', pdb)

# f2 = open(f'{outf}/{name}_interface_residues.txt', 'w')

interface_pairs = []

## Step 1. Get confidently predicted interface residue pairs

for model in models:
    chains = [c for c in model.get_chains()]
    for i, resa in enumerate(chains[0]):
        for j, resb in enumerate(chains[1]):
            for atoma, atomb in product(resa, resb):
                if atoma.get_bfactor() > 70 and atomb.get_bfactor() > 70:
                    if atoma-atomb < interface_cutoff:
                        interface_pairs.append((i, resa.get_id()[1], 'A', resa.center_of_mass(), j, resb.get_id()[1], 'B', resb.center_of_mass()))  #TODO add indx as well as residue number
                        break

# print(interface_pairs)
# Step 2. Get distance matrix to cluster interface residue pairs
# N = len(interface_pairs)
# distance_matrix = np.zeros((N, N))
# max_num = 100000
#
# for i in range(N):
#     for j in range(N):
#         if i != j:
#             dist1 = np.linalg.norm(interface_pairs[i][3] - interface_pairs[j][3])
#             dist2 = np.linalg.norm(interface_pairs[i][7] - interface_pairs[j][7])
#
#             if dist1 < intra_chain_residue_distance_cutoff and dist2 < intra_chain_residue_distance_cutoff:
#                 max_dist = max(dist1, dist2)
#                 distance_matrix[i][j] = max_dist
#             else:
#                 distance_matrix[i][j] = max_num
import numpy as np
import matplotlib.pyplot as plt

# Your code for calculating the distance matrix
N = len(interface_pairs)
distance_matrix = np.zeros((N, N))
max_num = 100000

for i in range(N):
    for j in range(N):
        if i != j:
            dist1 = np.linalg.norm(interface_pairs[i][3] - interface_pairs[j][3])
            dist2 = np.linalg.norm(interface_pairs[i][7] - interface_pairs[j][7])

            if dist1 < intra_chain_residue_distance_cutoff and dist2 < intra_chain_residue_distance_cutoff:
                max_dist = max(dist1, dist2)
                distance_matrix[i][j] = max_dist
            else:
                distance_matrix[i][j] = max_num

# Plotting the distance matrix
plt.imshow(distance_matrix, cmap='viridis', origin='upper', interpolation='nearest')
plt.colorbar(label='Distance')
plt.title('Distance Matrix')
plt.xlabel('Interface Residue Pairs')
plt.ylabel('Interface Residue Pairs')
plt.show()

# print(distance_matrix)

# Step 3. make clusters of the interface residues based on the distance matrix
affinity_propagation = AffinityPropagation(affinity='precomputed', random_state = 0 )
affinity_propagation.fit(-distance_matrix)  # Use negative distances as input

cluster_labels = affinity_propagation.labels_
n_clusters = len(set(cluster_labels))

# print(n_clusters)

#Step 4. calculate average PAE of all the interface pairs in each cluster and output is mean PAE
cluster_avg_pae =[[] for _ in range(n_clusters)]

for i in range(N): #each interface residue pair
    pae = data['predicted_aligned_error'][interface_pairs[i][0]][interface_pairs[i][4]] #TODO test this

    cluster_label = affinity_propagation.labels_[i]
    cluster_avg_pae[cluster_label].append(pae)

# print(len(cluster_avg_pae))
for i, cluster in enumerate(cluster_avg_pae,1):
    average = sum(cluster)/len(cluster)
    print(f'cluster_{i}', average)
