import matplotlib.pyplot as plt
import numpy as np
import os, sys
import glob
import RMF
import IMP
import IMP.rmf
import IMP.core
import IMP.atom


rmf_file = sys.argv[1]
mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_file) # input sampcon_0_extracted.rm3 file
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
num_frames = rmf_fh.get_number_of_frames()

mdl_ids = [i for i in range(num_frames)]

##specify domains
### Data not used in modeling -----------------------------------------------
# d1 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(15,20)).get_selected_particles()
# d2 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(1,79)).get_selected_particles()

# d3 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(84,104)).get_selected_particles()
# d4 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',copy_indexes = [0,1,2,3],residue_indexes = range(410,759)).get_selected_particles()

# d5 = IMP.atom.Selection(hierarchy=hier, molecule='MIC19',residue_indexes = range(172,222)).get_selected_particles()
# d6 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(410,759)).get_selected_particles()

# d7 = IMP.atom.Selection(hierarchy=hier, molecule='MIC19',residue_indexes = range(1,228)).get_selected_particles()
# d8 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',copy_indexes = [0,1,2,3],residue_indexes = range(410,759)).get_selected_particles()

# d9 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(46,53)).get_selected_particles()
# d10 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(1,79)).get_selected_particles()

# d11 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(24,29)).get_selected_particles()
# d12 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(1,79)).get_selected_particles()


# def calculate_distances(domain1, domain2, mdl_ids, rmf_file):
#     dist = []
#     for frame_id in range(num_frames):
#         IMP.rmf.load_frame(rmf_fh, frame_id)
#         mdl.update()
#         x=[]
#         for d1 in domain1:
#             for d2 in domain2:
#                 p1 = IMP.core.get_distance(IMP.core.XYZR(mdl, d1), IMP.core.XYZR(mdl, d2))
#                 if p1 < 0:
#                     p1 =0
#                 x.append(p1)
#         dist.append(min(x))

# #     # print(dist)
# #     # print('\n',len(dist))

#     return dist

# d1_d2 = calculate_distances(d1,d2,mdl_ids,rmf_file)
# d3_d4 = calculate_distances(d3,d4,mdl_ids,rmf_file)
# d5_d6 = calculate_distances(d5,d6,mdl_ids,rmf_file)
# d7_d8 = calculate_distances(d7,d8,mdl_ids,rmf_file)
# d9_d10 = calculate_distances(d9,d10,mdl_ids,rmf_file)
# d11_d12 = calculate_distances(d11,d12,mdl_ids,rmf_file)

# print('\nnow plotting')

# all_data = [d1_d2, d3_d4, d5_d6, d7_d8, d9_d10, d11_d12] 
# fig, ax = plt.subplots()

# tick_positions = np.arange(1, len(all_data) + 1)
# violinplot = ax.violinplot(all_data)

# ax.set_xticks(tick_positions)
# # ax.set_xticklabels(['MIC10-MIC60', 'MIC10-MIC60', 'MIC13-MIC10','MIC10-MIC10', 'MIC10-MIC10','MIC13-MIC60','MIC13-MIC60','MIC13-MIC10'], rotation=30, ha='right',fontsize=6) #'MIC10-MIC13', 'MIC13-MIC60',
# # ax.set_xticklabels(['MIC13-MIC60'], rotation=30, ha='right',fontsize=6)
# ax.set_xlabel('Domain Pairs')
# ax.set_ylabel('Minimum Distance')
# ax.set_title('Plot for MDBPR')
# # plt.show()
# plt.savefig('fit_to_binding_data_not_used', dpi=600)



### Data used in modeling -----------------------------------------------
def calculate_distances(domain1, domain2, num_frames, rmf_fh, mdl):
    dist = []
    for frame_id in range(num_frames):
        IMP.rmf.load_frame(rmf_fh, frame_id)
        mdl.update()
        distances = []
        for d1 in domain1:
            for d2 in domain2:
                p = IMP.core.get_distance(IMP.core.XYZR(mdl, d1), IMP.core.XYZR(mdl, d2))
                distances.append(p)
        dist.append(min(distances))
    return dist

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read CSV
df = pd.read_csv(sys.argv[2])

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

    print(p1, p2, start1, end1)
    # d1 = IMP.atom.Selection(hierarchy=hier, molecule=str(p1), residue_indexes=range(start1, end1 + 1)).get_selected_particles()
    # d2 = IMP.atom.Selection(hierarchy=hier, molecule=str(p2), residue_indexes=range(start2, end2 + 1)).get_selected_particles()
    d1 = IMP.atom.Selection(hierarchy=hier, molecule=str(p1),residue_indexes = range(15,20)).get_selected_particles()
    d2 = IMP.atom.Selection(hierarchy=hier, molecule=str(p2),residue_indexes = range(1,79)).get_selected_particles()
    print(d1, d2)
    exit()
    selections.append((d1, d2))
    labels.append(f'{p1}_{r1}_{p2}_{r2}')

distance_dict = {}

for label, (d1, d2) in zip(labels, selections):
    dist = calculate_distances(d1, d2, num_frames, rmf_fh, mdl)
    distance_dict[label] = dist


all_data = list(distance_dict.values())
label_names = list(distance_dict.keys())

fig, ax = plt.subplots(figsize=(12, 5))
violinplot = ax.violinplot(all_data, showmeans=True, showextrema=True, showmedians=True)

ax.set_xticks(np.arange(1, len(label_names) + 1))
ax.set_xticklabels(label_names, rotation=45, ha='right', fontsize=7)

ax.set_xlabel('Domain Pairs')
ax.set_ylabel('Minimum Distance (Å)')
ax.set_title('Minimum Distances from CSV-based Predictions')

plt.tight_layout()
plt.savefig('af3_csv_distances_violinplot.png', dpi=600)
# plt.show()
