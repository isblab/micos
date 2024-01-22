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
# d1 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(23,28)).get_selected_particles()
# d2 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(18,23)).get_selected_particles()

# d3 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(88,92)).get_selected_particles()
# d4 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',residue_indexes = range(657,668)).get_selected_particles()

d5 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(3,4)).get_selected_particles()
d6 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',residue_indexes = range(639,643)).get_selected_particles()

d7 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(2,4)).get_selected_particles()
d8 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(720,725)).get_selected_particles()

d9 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(81,86)).get_selected_particles()
d10 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(1,79)).get_selected_particles()

d11 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(24,29)).get_selected_particles()
d12 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(1,79)).get_selected_particles()

d13 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(46,53)).get_selected_particles()
d14 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(1,79)).get_selected_particles()

d15 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(2,27)).get_selected_particles()
d16 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',copy_indexes = [0,1,2,3],residue_indexes = range(410,759)).get_selected_particles()

d17 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(84,104)).get_selected_particles()
d18 = IMP.atom.Selection(hierarchy=hier, molecule='MIC60',copy_indexes = [0,1,2,3],residue_indexes = range(410,759)).get_selected_particles()

d19 = IMP.atom.Selection(hierarchy=hier, molecule='MIC13',residue_indexes = range(15,20)).get_selected_particles()
d20 = IMP.atom.Selection(hierarchy=hier, molecule='MIC10',residue_indexes = range(1,79)).get_selected_particles()



def calculate_distances(domain1, domain2, mdl_ids, rmf_file):
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

    # print(dist)
    # print('\n',len(dist))

    return dist


# d1_d2 = [calculate_distances(d1,d2, mdl_ids,rmf_file)]
# # d3_d4 = calculate_distances(d3,d4)
d5_d6 = calculate_distances(d5,d6,mdl_ids,rmf_file)
d7_d8 = calculate_distances(d7,d8,mdl_ids,rmf_file)
d9_d10 = calculate_distances(d9,d10,mdl_ids,rmf_file)
d11_d12 = calculate_distances(d11,d12,mdl_ids,rmf_file)
d13_d14 = calculate_distances(d13,d14,mdl_ids,rmf_file)
d15_d16 = calculate_distances(d15,d16,mdl_ids,rmf_file)
d17_18 = calculate_distances(d17,d18,mdl_ids,rmf_file)
d19_20 = calculate_distances(d19,d20,mdl_ids,rmf_file)

print('\nnow plotting')
# # print(d15_d16)
all_data = [d5_d6, d7_d8, d9_d10, d11_d12, d13_d14, d15_d16, d17_18, d19_20] #d1_d2, d3_d4,
fig, ax = plt.subplots()
#
tick_positions = np.arange(1, len(all_data) + 1)
violinplot = ax.violinplot(all_data)
#
ax.set_xticks(tick_positions)
# ax.set_xticklabels(['MIC10-MIC60', 'MIC10-MIC60', 'MIC13-MIC10','MIC10-MIC10', 'MIC10-MIC10','MIC13-MIC60','MIC13-MIC60','MIC13-MIC10'], rotation=30, ha='right',fontsize=6) #'MIC10-MIC13', 'MIC13-MIC60',
# ax.set_xticklabels(['MIC13-MIC60'], rotation=30, ha='right',fontsize=6)
ax.set_xlabel('Domain Pairs')
ax.set_ylabel('Minimum Distance')
ax.set_title('Plot for MDBPR')
# plt.show()
plt.savefig('fit_to_binding_data', dpi=600)
