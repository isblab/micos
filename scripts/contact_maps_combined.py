import os
import glob
import time
import shutil
import argparse

import RMF
import IMP
import IMP.rmf
import IMP.core
import IMP.atom
import numpy as np
from tqdm import tqdm
import dmaps_functions
import concurrent.futures
from matplotlib import pyplot as plt
import pandas as pd
import fit_to_binding_data_functions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Get contact maps for a given protein pair"
    )
    parser.add_argument(
        "--rmf_file",
        "-rf",
        dest="rfile",
        required=True,
        help="Sampcon cluster extracted .rmf3 file",
    )
    parser.add_argument(
        "--nprocs",
        "-p",
        dest="nprocs",
        help="Number of processes to be launched",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--resolution",
        "-r",
        dest="res",
        help="Resolution for computing distance maps",
        default=1,
        required=False,
    )
    parser.add_argument(
        "--binarization_threshold",
        "-t",
        dest="threshold",
        type=float,
        help="Distance threshold for identifying a contact",
        default=1,
        required=False,
    )

    return parser.parse_args()

def get_patches(arrays):
    patches = []
    single = []
    temp = []
    l2, master_cont = [], []

    i = 0
    while i < len(arrays[0]):
        value = arrays[0][i]
        temp.append((arrays[0][i], arrays[1][i]))
        while i + 1 < len(arrays[0]) and arrays[0][i + 1] == value:
            i += 1
            temp.append((arrays[0][i], arrays[1][i]))
        temp.sort(key=lambda x: x[1])
        l2.append(temp.copy())
        temp.clear()
        i += 1

    for l in range(len(l2)):
        for i in range(len(l2[l]) - 1):
            v1, v2 = l2[l][i]
            continous = []
            # print(v1, v2)
            for k in range(1, len(l2)):
                for m in range(len(l2[k])):
                    if l2[k][m][0] == v1+1 and l2[k][m][1] == v2+1:
                        # print(l2[k][m], 'nn')
                        continous.append((v1,v2)) # TODO remove redundant entries of v1, v2
                        continous.append(l2[k][m])
                        v1 = l2[k][m][0]
                        v2 = l2[k][m][1]
                        break

            else:
                if (v1,v2) not in single and (v1,v2) not in continous:
                    single.append((v1,v2))

            if len(continous) > 0:
                master_cont.append(continous)
            # elif len(continous) !=0:
            #     single.append(tuple(continous))

    for lis in master_cont:
        min_1 = min(lis, key = lambda x:x[0])
        max_1 = max(lis, key = lambda x:x[0])
        patches.append((min_1, max_1))

    return patches, single

def compute_dmaps():
    args = parse_args()

    for adirectory in [
        "distance_matrices",
        "distance_maps",
        "binarized_distance_matrices",
        "binarized_distance_maps",
    ]:
        if not os.path.isdir(os.path.join(os.getcwd(), adirectory)):
            os.mkdir(os.path.join(os.getcwd(), adirectory))

    tic = time.time()
    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(args.rfile)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)

    all_proteins = dmaps_functions.get_protein_names(hier)
    sizes_dict = dmaps_functions.get_protein_sizes(hier, all_proteins)
    print(all_proteins)
    # print(sizes_dict)
    # TODO will not work if the last particle in selection is not the last numbered bead.
    # TODO this can be the case for missing residues in PDB coming at the end of selection.
    # TODO updating this to take the min and max instead. Also handles the first res not being 1 issue.

    nmodels = rmf_fh.get_number_of_frames()
    mdl_ids = [i for i in range(nmodels)]

    del rmf_fh

    print(f"Total number of models: {nmodels}")
    dmaps_functions.split_models_into_subsets(args.rfile, mdl_ids, args.nprocs)

    concatenated_rmfs = glob.glob("concatenated_models/*rmf3")
    print("\nStarting distance calculations")

    i = 1
    done_prot_pairs = []
    protein_domain_dict = {}
    for p1 in all_proteins:
        for p2 in tqdm(all_proteins, desc=f"Processing interactions of {p1}"):
            if ((p1, p2) in done_prot_pairs) or ((p2, p1) in done_prot_pairs):
                # TODO  Shorter way to implement uses itertools.combinations or some other iterator over list
                # TODO to get pairs with yx
                """No xy -> yx repetitions"""
                continue

            print(
                f"\nProcessing the distance calculations for {p1}/{p2} pair. Hold on..."
            )
            s1, s2 = sizes_dict[p1], sizes_dict[p2]

            all_distances = []

            with concurrent.futures.ProcessPoolExecutor(args.nprocs) as executor:
                for distances in executor.map(
                    dmaps_functions.measure_beadwise_distances,
                    [p1 for _ in range(args.nprocs)],
                    [p2 for _ in range(args.nprocs)],
                    [s1 for _ in range(args.nprocs)],
                    [s2 for _ in range(args.nprocs)],
                    concatenated_rmfs,
                    [1 for _ in range(args.nprocs)],
                ):
                    all_distances.append(distances)

            all_distances = np.concatenate(all_distances, axis=2)
            protein_domain_dict[p1+'-'+p2] = all_distances
            # print(protein_domain_dict)
            l = fit_to_binding_data_functions.get_matrix(protein_domain_dict)


            exit()
            # mean_distances = all_distances.mean(axis=2)
            # mean_distances = np.delete(np.delete(mean_distances, 0, axis=0), 0, axis=1)
            # np.savetxt(
            #     os.path.join("distance_matrices", f"{p1}-{p2}_all_distances.csv"),
            #     mean_distances,
            # )
            # if p1 != p2:
            #     interacting_beads = dmaps_functions.get_interacting_beads(mean_distances)
            #     mapped_res_ = []
            #
            #     for k in range(2):
            #         mapped_array = [x + s1[0] if k == 0 else x + s2[0] for k in interacting_beads[k]]
            #         mapped_res_.append(mapped_array)
            #
            #     patches = get_patches(mapped_res_)
            #     df = pd.DataFrame(patches[0], columns=['start', 'end'])
            #     outfile = os.path.join("binarized_distance_matrices", f"{p1}-{p2}_merged_patches.csv")
            #     df.to_csv(outfile, index = False)
            #
                # np.savetxt(
                #     os.path.join("distance_matrices", f"{p1}-{p2}_mean_distances.csv"),
                #     mean_distances,
                # )
            #
            #
            # plt.figure(i, dpi=1200)
            # plt.imshow(mean_distances, cmap="hot")
            # plt.xlabel(p2)
            # plt.ylabel(p1)
            # plt.colorbar()
            # plt.savefig(
            #     os.path.join("distance_maps", f"{p1}-{p2}_dmap_w_new_script.png"),
            #     dpi=600,
            # )
            # plt.close(fig=i)
            #
            # i += 1
            # plt.figure(i, dpi=1200)
            # plt.imshow(binarized_distance_matrix, cmap="Greys")
            # plt.xlabel(p2)
            # plt.ylabel(p1)
            # # plt.colorbar()
            # plt.savefig(
            #     os.path.join(
            #         "binarized_distance_maps",
            #         f"{p1}-{p2}_binarized_dmap_w_new_script.png",
            #     ),
            #     dpi=600,
            # )
            # plt.close(fig=i)
            #
            # i += 1
            # done_prot_pairs.append((p1, p2))
            #     #TODO Remove these breaks
        #     break
        # break


    toc = time.time()
    print(
        f"Processed {len(done_prot_pairs)} pairs of proteins from {nmodels} models in {toc-tic} seconds"
    )
    shutil.rmtree(os.path.join(os.getcwd(), "concatenated_models"))
    if "__pycache__" in os.listdir("./"):
        shutil.rmtree(os.path.join(os.getcwd(), "__pycache__"))


compute_dmaps()

# dmat = np.loadtxt(
#     os.path.join("distance_matrices", "WDR76.0-WDR76.0_mean_distances.csv")
# )
# interfaces = np.where((dmat < 20.0) & (dmat != 0))
# print(interfaces[0].shape)
# print(dmat[interfaces])
