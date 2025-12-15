import IMP
import RMF
import IMP.rmf
import IMP.atom
import IMP.core
import IMP.algebra
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

import sys

rmf_file = sys.argv[1]

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_file)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
mdl.update()

nmodels = rmf_fh.get_number_of_frames()
mdl_ids = [i for i in range(nmodels)]


all_angles = []

for frame_id in tqdm(range(rmf_fh.get_number_of_frames())):
    IMP.rmf.load_frame(rmf_fh, frame_id)
    
    p1 = IMP.atom.Selection(
        hierarchy=hier, molecule='MIC13', residue_index = 80
    ).get_selected_particles()

    p2 = IMP.atom.Selection(
        hierarchy=hier, molecule='MIC13', residue_index = 95
    ).get_selected_particles()

    p3 = IMP.atom.Selection(
        hierarchy=hier, molecule='MIC13', residue_index = 9
    ).get_selected_particles()

    p4 = IMP.atom.Selection(
        hierarchy=hier, molecule='MIC13', residue_index = 23
    ).get_selected_particles()


    axis1 = IMP.core.XYZ(mdl, p1[0]).get_coordinates() - IMP.core.XYZ(mdl, p2[0]).get_coordinates()
    axis2 = IMP.core.XYZ(mdl, p3[0]).get_coordinates() - IMP.core.XYZ(mdl, p4[0]).get_coordinates()

    angle = np.degrees(np.arccos(np.dot(axis1, axis2)/(np.linalg.norm(axis1)*np.linalg.norm(axis2))))
    all_angles.append(angle)


plt.hist(all_angles, bins=20)
plt.xlabel('Angle (degrees)')
plt.ylabel('Number of models')
plt.savefig('Mic13_TM-C_helix_angle')
plt.show()
exit()