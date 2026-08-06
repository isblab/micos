#!/usr/bin/env python3

from collections import defaultdict
from Bio import AlignIO
from pymsaviz import MsaViz
import matplotlib
matplotlib.use("pdf")
import sys
from Bio.Align import MultipleSeqAlignment

############################################################
# SETTINGS
############################################################

MSA_FILE = sys.argv[1]

CONSURF_FILE = sys.argv[2]

REFERENCE_ID = "Human"

KEEP_SPECIES = [
    "Human",
    "Mouse",
    "Xenopus",
    "Danio",
    "Drosophila",
    "Caenorhabditis",
    "Chaetomium",
    "Lachancea",
    "Saccharomyces"
]

OUTPUT = "msa_consurf.pdf"

############################################################
# ConSurf 9-color palette
############################################################

CONSURF_COLORS = {
    1: "#10C8D1",
    2: "#8CFFFF",
    3: "#D7FFFF",
    4: "#EAFFFF",
    5: "#FFFFFF",
    6: "#FCEDF4",
    7: "#FAC9DE",
    8: "#F07DAB",
    9: "#A02560",
}

############################################################
# Read ConSurf grades
############################################################

grades = {}

with open(CONSURF_FILE) as f:
    for line in f:

        if not line.strip():
            continue

        if line.startswith("#"):
            continue

        cols = line.split()

        try:
            residue_number = int(cols[0])
            grade = int(cols[3].replace("*", ""))
        except:
            continue

        grades[residue_number] = grade

############################################################
# Read alignment
############################################################

msa = AlignIO.read(MSA_FILE, "fasta")

############################################################
# Keep selected sequences
############################################################

############################################################
# Keep selected species (skip missing ones)
############################################################

selected = []
found_species = []

for species in KEEP_SPECIES:
    for rec in msa:
        if species.lower() in rec.id.lower() or species.lower() in rec.description.lower():
            selected.append(rec)
            found_species.append(species)
            break   # only keep the first matching sequence

missing = [s for s in KEEP_SPECIES if s not in found_species]

print("Found species:")
for s in found_species:
    print(f"  ✓ {s}")

if missing:
    print("\nMissing species (skipped):")
    for s in missing:
        print(f"  - {s}")

if len(selected) == 0:
    raise ValueError("None of the requested species were found in the alignment.")

msa = MultipleSeqAlignment(selected)

############################################################
# Find reference sequence
############################################################

reference = None

for rec in msa:
    if REFERENCE_ID.lower() in rec.id.lower():
        reference = rec
        break

if reference is None:
    raise ValueError("Reference sequence not found.")

for rec in msa:
    # UniProt FASTA headers: sp|ACC|ENTRY_NAME or tr|ACC|ENTRY_NAME
    if "|" in rec.id:
        rec.id = rec.id.split("|")[-1]
    # Remove description after the first space (if any)
    rec.description = rec.id

############################################################
# residue number -> alignment column
############################################################

column_grade = {}

residue_number = 0

for col, aa in enumerate(str(reference.seq)):

    if aa != "-":
        residue_number += 1

        if residue_number in grades:
            column_grade[col] = grades[residue_number]

############################################################
# custom color function
############################################################

def consurf_color(row, col, aa, msa):

    if aa == "-":
        return "#FFFFFF"

    grade = column_grade.get(col)

    if grade is None:
        return "#FFFFFF"

    return CONSURF_COLORS[grade]

############################################################
# Plot
############################################################

mv = MsaViz(
    msa,
    color_scheme="None",
    show_grid=False,
    show_seq_char=True,
    wrap_length=80,
)

mv.set_custom_color_func(consurf_color)

fig = mv.plotfig()

fig.savefig(OUTPUT, dpi=300, bbox_inches="tight")