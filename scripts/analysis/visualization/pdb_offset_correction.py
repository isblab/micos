from Bio import PDB
import glob


########################################################################################
######################################## Inputs ########################################
########################################################################################

pdb_files = ['pdb/MIC60_tet.pdb',\
            'pdb/MIC60_19.pdb',\
            'pdb/MIC60.pdb',\
            'pdb/MIC10.pdb',\
            'pdb/MIC13.pdb']

offsets = [{'A':450,'B':303,'C':156,'D':9},\
            {'A':626,'B':52,'C':457,'D':-117 },\
            {'A':0},\
            {'A':0},\
            {'A':0}]


########################################################################################
###################################### Actual work #####################################
########################################################################################

for file_index in range(len(pdb_files)):
    pdbfile = pdb_files[file_index]
    offset = offsets[file_index]
    print(file_index,pdbfile,offsets[file_index])
    pdb_io = PDB.PDBIO()
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(" ", pdbfile)

    for model in structure:
        for chain in model:
            print(chain.id)
            for res in chain.get_residues():
                res_id_list = list(res.id)
                res_id_list[1] = res_id_list[1] + offset[chain.id]
                res.id = tuple(res_id_list)

    pdb_io.set_structure(structure)
    pdb_io.save('pdb/offset_corrected/offset_corrected_'+pdbfile.split('/')[-1])