import numpy as np
import sys
import pickle
from Bio.PDB import PDBParser
import json
from itertools import product
from sklearn.cluster import AffinityPropagation
import argparse

parser = argparse.ArgumentParser(description="predciting interface residues")
parser.add_argument("--interface_cutoff","-i", help="defining interface distance", default=10)
parser.add_argument("--residue_distance_cutoff","-d", help="defining distance between two residue pairs", default=10)
parser.add_argument("--path","-p", help="path to output of AF2 multimer", default=None)
parser.add_argument("--output_path","-o", help="output path", default=None)
parser.add_argument("--out_name","-n", help="output file name", default=None)

args = parser.parse_args()

path = args.path # path to output of AF2 multimer
outf = args.output_path # path to where output should be stored
name = args.out_name # name of complex e.g. dp-pg
interface_cutoff = int(args.interface_cutoff)
residue_distance_cutoff = int(args.residue_distance_cutoff)


pdb = f'{path}/ranked_0.pdb'
with open(f'{path}/ranking_debug.json', 'r') as f:
    o = json.loads(f.read())
best_pkl = o['order'][0]
pkl = f'{path}/result_{best_pkl}.pkl'
with open(pkl, 'rb') as f:
    data = pickle.load(f)

models = PDBParser().get_structure('pdb', pdb)

f2 = open(f'{outf}/{name}_interface_residues.txt', 'w')

interface_pairs = []

## reading the pdb structure

for model in models:
    chains = [c for c in model.get_chains()]
    for resa in (chains[0]):
        for resb in (chains[1]):
            for atoma, atomb in product(resa, resb):
                if atoma.get_bfactor() > 70 and atomb.get_bfactor() > 70:
                    if atoma-atomb < interface_cutoff:
                        interface_pairs.append((resa.get_id()[1], 'A', resa.center_of_mass(), resb.get_id()[1], 'B', resb.center_of_mass()))
                        break

# print(interface_pairs)
# calculating half the traingle ?? not done yet

N = len(interface_pairs)
distance_matrix = np.zeros((N, N))
max_num = 100000

for i in range(N):
    for j in range(N):
        if i != j:
            dist1 = np.linalg.norm(interface_pairs[i][2] - interface_pairs[j][2])
            dist2 = np.linalg.norm(interface_pairs[i][5] - interface_pairs[j][5])

            if dist1 < residue_distance_cutoff and dist2 < residue_distance_cutoff:
                max_dist = max(dist1, dist2)
                distance_matrix[i][j] = max_dist
            else:
                distance_matrix[i][j] = max_num

# print(distance_matrix)

### affinity propogation clustering ####
affinity_propagation = AffinityPropagation(affinity='precomputed', random_state = 0 )
affinity_propagation.fit(-distance_matrix)  # Use negative distances as input

cluster_labels = affinity_propagation.labels_
n_clusters = len(set(cluster_labels))

clustered_interface_pairs = [[] for _ in range(n_clusters)]

for i in range(N):
    cluster_label = cluster_labels[i]
    j = i + 1
    while j < N:
        if cluster_labels[j] == cluster_label:
            clustered_interface_pairs[cluster_label].append((i, j))
        j += 1

### Kd Tree clustering



## function to calculate pae in both ways
def calculate_pae(patch):
# all vs all PAE
    var = 0
    count = 0
    pae = 0
    confident_patches_all = []
    confident_patches_pairwise = []

    for i in range(len(patch)):
        for j in range(len(patch)):
            var += data['predicted_aligned_error'][patch[i][0]][patch[j][2]]
            count += 1

    avg_pae_all = var / count

    for i in range(len(patch)):
        pae += data['predicted_aligned_error'][patch[i][0]][patch[i][2]]

# pairwise
    avg_pae_pairwise = pae/len(patch)
    # print(avg_pae_pairwise)

    for pair in patch:
        if data['predicted_aligned_error'][pair[0]][pair[2]] < avg_pae_all:
            confident_patches_all.append((pair[0], pair[2], data['predicted_aligned_error'][pair[0]][pair[2]]))

        if data['predicted_aligned_error'][pair[0]][pair[2]] < avg_pae_pairwise:
            confident_patches_pairwise.append((pair[0], pair[2], data['predicted_aligned_error'][pair[0]][pair[2]]))


    return len(confident_patches_all), len(confident_patches_pairwise)


## to make one list without duplicate entries in each cluster and calculate pae of that cluster
for cluster_label, pairs_list in enumerate(clustered_interface_pairs):
    print(f"Cluster {cluster_label}:")
    patches = []
    for pairs in pairs_list:
        i, j = pairs
        pair1 = (interface_pairs[i][0], interface_pairs[i][1], interface_pairs[i][3], interface_pairs[i][4])
        pair2 = (interface_pairs[j][0], interface_pairs[j][1], interface_pairs[j][3], interface_pairs[j][4])

        if pair1 not in patches:
            patches.append(pair1)

        if pair2 not in patches:
            patches.append(pair2)

    print(calculate_pae(patches))


#TODO implement all vs all PAE and interface pairs?
#TODO implement 5 A or 10 A
