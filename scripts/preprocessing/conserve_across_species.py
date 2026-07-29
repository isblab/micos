import argparse
import glob
import os

import pandas as pd
from Bio import AlignIO


##############################################################################
# Similar amino acid groups
##############################################################################

SIMILARITY_GROUPS = [
    set("STA"),
    set("NEQK"),
    set("NHQK"),
    set("NDEQ"),
    set("QHRK"),
    set("MILV"),
    set("MILF"),
    set("FYW"),
    set("HY"),
]


##############################################################################
# Utilities
##############################################################################

def find_alignment_file(protein, alignment_dir):
    """
    Search for an alignment file containing the protein name.
    """

    exts = ("*.fa", "*.faa", "*.fasta", "*.fas", "*.aln", "*.txt")

    for ext in exts:
        for f in glob.glob(os.path.join(alignment_dir, ext)):
            if protein.upper() in os.path.basename(f).upper():
                return f

    return None


def find_human_record(alignment):
    """
    Locate Homo sapiens sequence.
    """

    for rec in alignment:
        desc = rec.description

        if "OS=Homo sapiens" in desc:
            return rec

        if "_HUMAN" in rec.id.upper():
            return rec

    raise ValueError("No human sequence found.")


def residue_to_alignment_column(sequence, residue_number):
    """
    Convert human residue numbering to alignment column.
    """

    count = 0

    for i, aa in enumerate(sequence):

        if aa != "-":
            count += 1

        if count == residue_number:
            return i

    return None


def classify_site(residues):
    """
    Classify alignment column.
    """

    residues = [x for x in residues if x != "-"]

    if len(residues) == 0:
        return "Gap", 0.0

    unique = set(residues)

    reference = residues[0]

    identity = 100 * sum(r == reference for r in residues) / len(residues)

    if len(unique) == 1:
        return "Identical", identity

    for group in SIMILARITY_GROUPS:
        if unique.issubset(group):
            return "Conserved", identity

    return "Variable", identity


##############################################################################
# Analysis
##############################################################################

def analyse_position(protein, residue_number, alignment_dir):

    aln_file = find_alignment_file(protein, alignment_dir)
    fmt = "fasta"
    alignment = AlignIO.read(aln_file, fmt)

    human = find_human_record(alignment)

    col = residue_to_alignment_column(
        str(human.seq),
        residue_number
    )

    if col is None:
        raise ValueError(
            f"{protein}: residue {residue_number} not present "
            "in human sequence."
        )
    
    species = []

    # ------------------------------------------------
    # Compare Human vs Chaetomium thermophilum (CT)
    # ------------------------------------------------

    human_aa = human.seq[col]
    ct_record = None

    for rec in alignment:

        if "Chaetomium thermophilum" in rec.description:
            ct_record = rec
            break

    if ct_record is None:
        raise ValueError(
            "Chaetomium thermophilum sequence not found"
        )

    ct_aa = ct_record.seq[col]


    # Pairwise classification
    if human_aa == "-" or ct_aa == "-":
        status = "Gap"
        identity = 0
        print('gap')
    elif human_aa == ct_aa:
        status = "Identical"
        identity = 100
        print('identical')
    else:
        pair = {human_aa, ct_aa}
        conserved = False

        for group in SIMILARITY_GROUPS:
            if pair.issubset(group):
                conserved = True
                break

        if conserved:
            status = "Conserved"
        else:
            status = "Variable"
        identity = 0


    species = [
        f"Homo sapiens:{human_aa}",
        f"Chaetomium thermophilum:{ct_aa}"
    ]

    return {
        "AA": human.seq[col],
        "Status": status,
        "Identity": round(identity, 1),
        "SpeciesResidues": "; ".join(species),
    }


##############################################################################
# Main
##############################################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a",
        "--alignments",
        required=True,
        help="Folder containing FASTA/MSA files"
    )

    parser.add_argument(
        "-c",
        "--crosslinks",
        required=True,
        help="Crosslink CSV"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="crosslink_conservation.csv"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.crosslinks)
    
    # ---------------------------------------------------
    # Collect unique protein-residue sites from crosslinks
    # ---------------------------------------------------

    unique_sites = set()

    for _, row in df.iterrows():
        
        unique_sites.add((row["Protein1"], int(row["Residue1"])))
        unique_sites.add((row["Protein2"], int(row["Residue2"])))
    
    print("Number of unique sites:")
    print(len(unique_sites))

    for protein in sorted(set(x[0] for x in unique_sites)):
        print(
            protein,
            len([x for x in unique_sites if x[0] == protein])
        )

    summary = []

    for protein, residue in sorted(unique_sites):
        print(protein, residue)
        
        result = analyse_position(protein, residue, args.alignments)

        summary.append({
            "Protein": protein,
            "Residue": residue,
            "Human_AA": result["AA"],
            "Status": result["Status"],
            "Identity (%)": result["Identity"],
            "Species_Residues": result["SpeciesResidues"]
        })

    summary_df = pd.DataFrame(summary)
    summary_df = summary_df.sort_values(["Protein", "Residue"])

    summary_file = args.output.replace(".csv", "_protein_summary.csv")
    summary_df.to_csv(summary_file, index=False)

    print(f"Protein summary written to {summary_file}")

    # Protein-wise statistics
    protein_stats = (
        summary_df
        .groupby("Protein")["Status"]
        .value_counts()
        .unstack(fill_value=0)
    )

    # Ensure all columns exist
    for col in ["Identical", "Conserved", "Variable"]:
        if col not in protein_stats.columns:
            protein_stats[col] = 0

    protein_stats["Total Crosslinked Residues"] = (
        protein_stats["Identical"] +
        protein_stats["Conserved"] +
        protein_stats["Variable"]
    )

    protein_stats = protein_stats[
        ["Total Crosslinked Residues",
        "Identical",
        "Conserved",
        "Variable"]
    ].reset_index()

    stats_file = args.output.replace(
        ".csv",
        "_protein_statistics.csv"
    )

    protein_stats.to_csv(stats_file, index=False)

    print("\nProtein-wise statistics")
    print(protein_stats.to_string(index=False))
    print(f"\nStatistics written to {stats_file}")


if __name__ == "__main__":
    main()