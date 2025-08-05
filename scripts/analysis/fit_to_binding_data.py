import matplotlib.pyplot as plt
import numpy as np
import os, sys
import glob
import RMF
import IMP
import IMP.rmf
import IMP.core
import IMP.atom
import pandas as pd


def calculate_distances(domain1, domain2):
    dist = []
    for frame_id in range(num_frames):
        IMP.rmf.load_frame(rmf_fh, frame_id)
        mdl.update()
        x=[]
        for d1 in domain1:
            for d2 in domain2:
                p1 = IMP.core.get_distance(IMP.core.XYZR(mdl, d1), IMP.core.XYZR(mdl, d2))
                if p1 < 0:
                    p1 =0
                x.append(p1)
        dist.append(min(x))
    return dist

def plot_fit_to_data(df, output, ):
    selections = []
    labels = []

    for idx, row in df.iterrows():
        p1 = row['Protein1']
        p2 = row['Protein2']
        r1 = row['Region1']
        r2 = row['Region2']
        
        if '-' in str(r1):
            start1, end1 = map(int, r1.split('-'))
        else:
            start1 = end1 = int(r1)

        if '-' in str(r2):
            start2, end2 = map(int, r2.split('-'))
        else:
            start2 = end2 = int(r2)

        d1 = IMP.atom.Selection(hierarchy=hier, molecule=str(p1), residue_indexes=range(start1, end1 + 1)).get_selected_particles()
        d2 = IMP.atom.Selection(hierarchy=hier, molecule=str(p2), residue_indexes=range(start2, end2 + 1)).get_selected_particles()
        selections.append((d1, d2))
        labels.append(f'{p1}_{r1}_{p2}_{r2}')

    distance_dict = {}

    for label, (d1, d2) in zip(labels, selections):
        dist = calculate_distances(d1, d2)
        distance_dict[label] = dist


    all_data = list(distance_dict.values())
    label_names = list(distance_dict.keys())

    fig, ax = plt.subplots(figsize=(12, 5))
    violinplot = ax.violinplot(all_data, showmeans=True, showextrema=True, showmedians=True)

    ax.set_xticks(np.arange(1, len(label_names) + 1))
    ax.set_xticklabels(label_names, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Minimum Distance (Å)')

    plt.tight_layout()
    plt.savefig(f'{output}.png', dpi=600)
    # plt.show()


rmf_file = sys.argv[1]
data_used = pd.read_csv(sys.argv[2]) # csv file for biochemical data and AF3 predictions used in modeling
data_not_used = pd.read_csv(sys.argv[3]) # csv file for biochemical data not used in modeling

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_file) # input sampcon_0_extracted.rm3 file
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
num_frames = rmf_fh.get_number_of_frames()

mdl_ids = [i for i in range(num_frames)]

### Data used in modeling -----------------------------------------------
plot_fit_to_data(data_used, 'fit_to_data_used')

### Data not used in modeling -----------------------------------------------
plot_fit_to_data(data_not_used, 'fit_to_data_not_used')