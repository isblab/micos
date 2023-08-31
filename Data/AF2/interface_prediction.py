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
parser.add_argument("--interface_cutoff","-i", help="defining interface distance", default=10)
parser.add_argument("--residue_distance_cutoff","-d", help="defining distance between two residue pairs", default=10)
parser.add_argument("--path","-p", help="path to output of AF2 multimer", default=None)
parser.add_argument("--output_path","-o", help="output path", default=None)
parser.add_argument("--out_name","-n", help="output file name", default=None)
parser.add_argument("--method","-m", help="which method", default='algorithm')

args = parser.parse_args()

path = args.path # path to output of AF2 multimer
outf = args.output_path # path to where output should be stored
name = args.out_name # name of complex e.g. dp-pg
interface_cutoff = int(args.interface_cutoff)
residue_distance_cutoff = int(args.residue_distance_cutoff)
method = args.method


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
                    if method == 'PICKluster':
                        interface_pairs.append((resa.get_id()[1], 'A', resa.center_of_mass()))
                        interface_pairs.append((resb.get_id()[1], 'B', resb.center_of_mass()))
                        break

                    if method == 'algorithm':
                        if atoma-atomb < interface_cutoff:
                            interface_pairs.append((resa.get_id()[1], 'A', resa.center_of_mass(), resb.get_id()[1], 'B', resb.center_of_mass()))
                            break

# print(interface_pairs)

# function to calculate pae in both ways
def calculate_pae(patch, method):

    var = 0
    count = 0
    pae = 0
    confident_patches_all = []
    confident_patches_pairwise = []


# all vs all PAE
    for i in range(len(patch)):
        for j in range(len(patch)):
            if method == 'algorithm':
                var += data['predicted_aligned_error'][patch[i][0]][patch[j][2]]
            elif method == 'PICKluster':
                var += data['predicted_aligned_error'][patch[i][0]][patch[j][0]]
            count += 1

    avg_pae_all = var / count

# pairwise
    #
    # for i in range(len(patch)):
    #     pae += data['predicted_aligned_error'][patch[i][0]][patch[i][2]]
    # avg_pae_pairwise = pae/len(patch)

# using the average pae as cutoff to get confident patches
    if len(patch) >1 :

        for pair in patch:
            # print(pair)
            if method == 'algorithm':
                if data['predicted_aligned_error'][pair[0]][pair[2]] < avg_pae_all:
                    confident_patches_all.append((pair[0], pair[2], data['predicted_aligned_error'][pair[0]][pair[2]]))

                if data['predicted_aligned_error'][pair[0]][pair[2]] < avg_pae_pairwise:
                    confident_patches_pairwise.append((pair[0], pair[2], data['predicted_aligned_error'][pair[0]][pair[2]]))

                    return len(confident_patches_all), len(confident_patches_pairwise)


            elif method == 'PICKluster':
                if data['predicted_aligned_error'][patch[0][0]][patch[1][0]] < avg_pae_all:
                    confident_patches_all.append((patch, data['predicted_aligned_error'][patch[0][0]][patch[1][0]]))

                    return confident_patches_all


# calculating half the traingle ?? not done yet
if method == 'algorithm':
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

    ## to make one list without duplicate entries in each cluster and calculate pae of that cluster

    for cluster_label, pairs_list in enumerate(clustered_interface_pairs):
        print(f"Cluster {cluster_label}:")
        patches = []
        for pairs in pairs_list:
            i, j= pairs
            pair1 = (interface_pairs[i][0], interface_pairs[i][1], interface_pairs[i][3], interface_pairs[i][4])
            pair2 = (interface_pairs[j][0], interface_pairs[j][1], interface_pairs[j][3], interface_pairs[j][4])

            if pair1 not in patches:
                patches.append(pair1)

            if pair2 not in patches:
                patches.append(pair2)

        print(calculate_pae(patches, 'algorithm'))


if method == 'PICKluster':
    interface_coords = [pair[2] for pair in interface_pairs]
    interface_coords_array = np.array(interface_coords)
    kdtree = KDTree(interface_coords_array)
    clusters = kdtree.query_ball_tree(kdtree, interface_cutoff)

    # Initialize a list to store clusters
    clustered_interface_pairs = []

    for i, cluster in enumerate(clusters):
        print(f"Cluster {i + 1}:")
        cluster_entries = []
        for idx in cluster:
            if (interface_pairs[idx][0],interface_pairs[idx][1]) not in cluster_entries:
                cluster_entries.append((interface_pairs[idx][0], interface_pairs[idx][1]))

        print(calculate_pae(cluster_entries, 'PICKluster'))
