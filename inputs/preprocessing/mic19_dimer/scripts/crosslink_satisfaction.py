from Bio import PDB
import pandas as pd
import sys, os

input_file = sys.argv[1]
output_file = sys.argv[2]
violation_ = float(sys.argv[3])
pdb_file = sys.argv[4]
flag = sys.argv[5]

def calculate_distance(pdb_file, residue1, residue2, violation):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    try:
        ca_atom1 = structure[0]['A'][residue1]['CA']
        ca_atom2 = structure[0]['B'][residue2]['CA']

        distance_r1_A_r2_B = ca_atom1 - ca_atom2  # this single distance is enough for self links

        min_distance = distance_r1_A_r2_B
        chainA, chainB = "A", "B"

        if flag == 'monomeric':
            ca_atom1 = structure[0]['A'][residue2]['CA']
            ca_atom2 = structure[0]['B'][residue1]['CA']
            distance_r1_B_r2_A = ca_atom1 - ca_atom2

            if distance_r1_B_r2_A < min_distance:
                min_distance = distance_r1_B_r2_A
                chainA = "B"
                chainB = "A"

            ca_atom1 = structure[0]['A'][residue1]['CA']
            ca_atom2 = structure[0]['A'][residue2]['CA']
            distance_r1_A_r2_A = ca_atom1 - ca_atom2

            if distance_r1_A_r2_A < min_distance:
                min_distance = distance_r1_A_r2_A
                chainA = "A"
                chainB = "A"

            ca_atom1 = structure[0]['B'][residue1]['CA']
            ca_atom2 = structure[0]['B'][residue2]['CA']
            distance_r1_B_r2_B = ca_atom1 - ca_atom2

            if distance_r1_B_r2_B < min_distance:
                min_distance = distance_r1_B_r2_B
                chainA = "B"
                chainB = "B"

        if min_distance < violation:
            return min_distance, chainA, chainB, True
        else:
            return min_distance, chainA, chainB, False

    except KeyError:
        return None

df = pd.read_csv(input_file)
df2 = pd.DataFrame(columns=["Residue1", "Residue2", "Distance"])
df3 = pd.DataFrame(columns=["Residue1", "Residue2", "Distance"])

if __name__ == "__main__":
    pdb_file = f"../servers/{pdb_file}"
    count_not_present = 0
    rows2 = []
    rows3 = []

    for index, row in df.iterrows():
        res1 = int(row['Residue1'])
        res2 = int(row['Residue2'])

        if 59 < res1 < 175 and 59 < res2 < 175:
            result = calculate_distance(pdb_file, res1, res2, violation_)

            if result is not None:
                distance, chainA, chainB, not_violated = result
                if not_violated:
                    rows2.append((res1, res2, distance))
                else:
                    rows3.append((res1, res2, distance))
            else:
                count_not_present += 1

    df2 = pd.DataFrame(rows2, columns=["Residue1", "Residue2", "Distance"])
    df3 = pd.DataFrame(rows3, columns=["Residue1", "Residue2", "Distance"])
    
print('not present in the pdb', count_not_present)
try:
    print(output_file, (len(df2) / (len(df) - count_not_present)) * 100)
except ZeroDivisionError:
    print('0%')
