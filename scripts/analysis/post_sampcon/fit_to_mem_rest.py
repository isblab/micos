import matplotlib.pyplot as plt
import math
import sys
import RMF
import IMP
import IMP.rmf
import IMP.core
import IMP.atom
import pandas as pd
import numpy as np
import tqdm

def calculate_distances(domain1, rest_type):

    dist = []
    for frame_id in tqdm.tqdm(range(num_frames)):
        IMP.rmf.load_frame(rmf_fh, frame_id)
        mdl.update()
        all_dist=[]
        for d1 in domain1:
            for d in d1:
                x,y = IMP.core.XYZR(mdl, d).get_x(), IMP.core.XYZR(mdl, d).get_y()
                p1 = math.sqrt(x**2 + y**2)

                if rest_type == 'ilr':
                    if p1 < r:
                        p1 =0

                elif rest_type == 'tlr':
                    if R>p1>r:
                        p1 =0

                elif rest_type == 'mlr':
                    if p1 > R:
                        p1 =0
                elif rest_type == 'zar_above':
                    p1 = IMP.core.XYZR(mdl, d).get_z()
                    if max_above_cj<p1< on_cj:
                        p1 =0
                elif rest_type == 'zar_oncj':
                    p1 = IMP.core.XYZR(mdl, d).get_z()
                    if on_cj<p1<max_in_cj:
                        p1 =0
                else:
                    exit('define thresholds')

                all_dist.append(p1)
        dist.append(min(all_dist))

    return dist

def get_particles(domains, region):
    out = []
    for prot, regions in domains.items():
        try:
            for start, end in regions[region]:
                residue_range = range(start, end + 1)
                if prot == "MIC19":
                    copy_indexes = [0, 1] if end <= 185 else [0, 2]
                elif prot == "MIC60":
                    copy_indexes = [0, 1, 2, 3]
                elif prot == "MIC10":
                    copy_indexes = [0,1]
                else:
                    copy_indexes = [0]

                sel = IMP.atom.Selection(
                    hierarchy=hier,
                    molecule=prot,
                    residue_indexes=residue_range,
                    copy_indexes=copy_indexes
                ).get_selected_particles()

                out.append(sel)
        
        except KeyError:
            print('Domain not present')

    return out


rmf_file = sys.argv[1]

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_file) # input sampcon_0_extracted.rm3 file
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
num_frames = rmf_fh.get_number_of_frames()

mdl_ids = [i for i in range(num_frames)]

R = 165  # radius of the outer cylinder
r = 125  # radius of the inner cylinder
on_cj = 0
max_in_cj = 70
max_above_cj = -50

# Membrane regions of each modeled protein
domains = {'MIC10': {'IMS': [[1,12], [61, 78]], 'TM': [[13,36], [40, 60]], 'Matrix': [[37, 39]], 'zar_below': [[1,12], [61,78]]}, 
            'MIC13': {'IMS': [[24, 118]], 'TM': [[8, 23]], 'Matrix': [[1,7]], 'zar_below': [[1, 118]]} , 
            'MIC19': {'IMS': [[1,227]], 'zar_sam50': [[1,14]]},
            'MIC60': {'IMS': [[410,758]], 'zar_mic60_cc': [[410,582]], 'zar_below': [[410,758]]}}


ims_particles = get_particles(domains, 'IMS')
tm_particles = get_particles(domains, 'TM')
matrix_particles = get_particles(domains, 'Matrix')
zar_sam50_particles = get_particles(domains, 'zar_sam50')
zar_mic60_cc_particles = get_particles(domains, 'zar_mic60_cc')
zar_below_particles = get_particles(domains, 'zar_below')

# print(zar_sam50_particles, zar_mic60_cc_particles, zar_below_particles)

score_ilr = calculate_distances(ims_particles, 'ilr')
score_mlr = calculate_distances(matrix_particles, 'mlr')
score_tlr = calculate_distances(tm_particles, 'tlr')
score_zar_above = calculate_distances(zar_sam50_particles, 'zar_above')
score_zar_oncj = calculate_distances(zar_mic60_cc_particles, 'zar_oncj')
score_zar_below = calculate_distances(zar_below_particles, 'zar_oncj')

def report_3sd(name, scores):
    scores = np.asarray(scores)
    mean = scores.mean()
    std = scores.std()
    within = (scores >= mean - 3*std) & (scores <= mean + 3*std)

    n_total = len(scores)
    n_within = within.sum()
    percent = 100 * n_within / n_total

    return f"{name}: {n_within}/{n_total} ({percent:.2f}%) within ±3 SD"


print(
    report_3sd("ILR", score_ilr), "\n",
    report_3sd("MLR", score_mlr), "\n",
    report_3sd("TLR", score_tlr), "\n",
    report_3sd("zar above", score_zar_above), "\n",
    report_3sd("zar on CJ", score_zar_oncj), "\n",
    report_3sd("zar below", score_zar_below),
)