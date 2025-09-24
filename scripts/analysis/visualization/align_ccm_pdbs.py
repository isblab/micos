import os, sys
import IMP
import RMF
import IMP.core
import IMP.rmf
import IMP.pmi.analysis
import glob
from Bio.PDB import *


###################################################################################################
############################################# Inputs ##############################################
###################################################################################################

input_file = sys.argv[1]

pdb_files = ['offset_corrected/offset_corrected_MIC60_tet.pdb',\
            'offset_corrected/offset_corrected_MIC60_19.pdb',\
            'offset_corrected/offset_corrected_MIC60_0.pdb',\
            'offset_corrected/offset_corrected_MIC60_1.pdb',\
            'offset_corrected/offset_corrected_MIC19_0.pdb',\
            'offset_corrected/offset_corrected_MIC19_1.pdb',\
            'offset_corrected/offset_corrected_MIC10_0.pdb',\
            'offset_corrected/offset_corrected_MIC10_1.pdb',\
            'offset_corrected/offset_corrected_MIC10_dimer_TM_MIC13_TM.pdb',\
            'offset_corrected/offset_corrected_MIC13_0.pdb',\
            'offset_corrected/offset_corrected_MIC13_1.pdb',\
            'offset_corrected/offset_corrected_MIC13_2.pdb']


# The all proteins list has the following architecture:
# [{protein:{chain_id,residue range}}, {protein:{chain_id,residue range}]
# The order of entries in the offset list must be the same as that in the pdb_files list

all_proteins = [{'MIC60':{'A':[0,range(410,583)]},'MIC60':{'B':[1,range(410,583)]},'MIC60':{'C':[2,range(410,583)]},'MIC60':{'D':[3,range(410,583)]}},\
            {'MIC60':{'A':[0,range(627,649)]},'MIC60':{'A':[0,range(683,759)]},'MIC60':{'C':[1,range(627,649)]},'MIC60':{'C':[1,range(683,759)]},'MIC19':{'B':[0,range(186,227)]},'MIC19':{'D':[2,range(186,227)]}},\
            {'MIC60':{'A':[0,range(588,627)]}},\
            {'MIC60':{'A':[1,range(588,627)]}},\
            {'MIC19':{'A':[0,range(59,175)]}},\
            {'MIC19':{'A':[1,range(59,175)]}},\
            {'MIC10':{'A':[0,range(2,13)]}},\
            {'MIC10':{'A':[1,range(2,13)]}},\
            {'MIC10':{'A':[0,range(13,37)]},'MIC10':{'A':[0,range(40,61)]}, 'MIC10':{'B':[1,range(13,37)]},'MIC10':{'B':[0,range(40,61)]}, 'MIC13':{'C':[0,range(8,24)]}},\
            {'MIC13':{'A':[0,range(2,8)]}},\
            {'MIC13':{'A':[0,range(24,69)]}},\
            {'MIC13':{'A':[0,range(79,118)]}}]


###################################################################################################
##################################### Get transformations #########################################
###################################################################################################

for file_index in range(len(pdb_files)):
    print(f"Aligning: {pdb_files[file_index]}")
    pdb_file = pdb_files[file_index]
    proteins = all_proteins[file_index]
    # print(proteins)
    ccm_mdl = IMP.Model()
    ccm = RMF.open_rmf_file_read_only(input_file)
    hier = IMP.rmf.create_hierarchies(ccm, ccm_mdl)[0]
    IMP.rmf.load_frame(ccm, 0)
    ccm_mdl.update()
    pdb_ca_mdl = IMP.Model()
    pdb_ca = IMP.atom.read_pdb(pdb_file,pdb_ca_mdl,IMP.atom.CAlphaPDBSelector())
    pdb_ca_mdl.update()

    new_mdl = IMP.Model()
    reload = IMP.atom.read_pdb(pdb_file, new_mdl)

    coords_pdb_ca = {}
    coords_ccm = {}

    for prot in proteins.keys():
        for chain_id in proteins[prot]:
            protein_name = prot

            sel_ca_pdb = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_id,residue_indexes=[i for i in proteins[prot][chain_id][1]]).get_selected_particles()
            sel_ccm = IMP.atom.Selection(hier,resolution=1,molecule=protein_name,copy_index=proteins[prot][chain_id][0],residue_indexes=[i for i in proteins[prot][chain_id][1]]).get_selected_particles()
            print(len(sel_ca_pdb),len(sel_ccm))
            print(protein_name, proteins[prot][chain_id][0], proteins[prot][chain_id][1])

            # Remove coarse grained beads
            new_ccm_sel = []
            for selection in sel_ccm:
                if not IMP.atom.Fragment.get_is_setup(selection):
                    new_ccm_sel.append(selection)
                    # print(selection)

            # print(len(sel_ca_pdb),'\n\n', len(new_ccm_sel))
            # for i in range(len(new_ccm_sel)):
                # print(new_ccm_sel[i],sel_ca_pdb[i])


            coords_pdb_ca[protein_name] = [IMP.core.XYZ(i).get_coordinates() for i in sel_ca_pdb]
            coords_ccm[protein_name] = [IMP.core.XYZ(i).get_coordinates() for i in new_ccm_sel]
            print(len(coords_pdb_ca[protein_name]),len(coords_ccm[protein_name]))
    _, transformation = IMP.pmi.analysis.Alignment(query=coords_pdb_ca, template=coords_ccm).align()
    print(transformation)


    ###################################################################################################
    #################################### Transform and write PDB ######################################
    ###################################################################################################

    IMP.atom.transform(reload, transformation)
    IMP.atom.write_pdb(reload, f"./aligned_{pdb_file.split('/')[-1]}.pdb")
