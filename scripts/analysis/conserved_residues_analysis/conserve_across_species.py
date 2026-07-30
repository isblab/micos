import argparse
import glob
import os

import pandas as pd
from Bio import AlignIO


# Similar amino acid groups from BLOSSUM
SIMILARITY_GROUPS = [
    set("C"),
    set("GSTAP"),
    set("DENQ"),
    set("RHK"),
    set("LVMI"),
    set("FWY")]


def find_alignment_file(protein, alignment_dir):
    """
    Search for an alignment file containing the protein name.
    """

    for f in glob.glob(os.path.join(alignment_dir,  "*.txt")):
        filename = os.path.basename(f).upper()
        if filename.startswith(protein.upper()):
            return f

    return None


def find_record(alignment, species):
    """
    Locate a sequence record for a given species.

    Parameters
    ----------
    alignment : Bio.Align.MultipleSeqAlignment
    species : str
        Species name, e.g.
        "Homo sapiens"
        "Mus musculus"
        "Chaetomium thermophilum"
        "Lachancea thermotolerans"
    """

    for rec in alignment:

        if f"OS={species}" in rec.description:
            return rec

    raise ValueError(f"No sequence found for {species}.")


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
    alignment = AlignIO.read(aln_file, "fasta")

    human = find_record(alignment, "Homo sapiens")

    col = residue_to_alignment_column(
        str(human.seq),
        residue_number
    )

    human_aa = human.seq[col]

    mouse_record = find_record(alignment, "Mus musculus")
    lt_record = find_record(alignment, "Lachancea thermotolerans")
    ct_record = find_record(alignment, "Chaetomium thermophilum")

    mouse_aa = mouse_record.seq[col] if mouse_record else None
    lt_aa = lt_record.seq[col] if lt_record else None
    ct_aa = ct_record.seq[col] if ct_record else None

    # ------------------------------------------------
    # Pairwise classification helper
    # ------------------------------------------------

    def classify_pair(ref_aa, other_aa):

        if other_aa is None:
            return "Missing"

        if ref_aa == "-" or other_aa == "-":
            return "Gap"

        if ref_aa == other_aa:
            return "Identical"

        pair = {ref_aa, other_aa}

        for group in SIMILARITY_GROUPS:
            if pair.issubset(group):
                return "Conserved"

        return "Variable"


    mouse_status = classify_pair(human_aa, mouse_aa)
    lt_status = classify_pair(human_aa, lt_aa)
    ct_status = classify_pair(human_aa, ct_aa)

    # ------------------------------------------------
    # Classification across ALL species
    # ------------------------------------------------

    all_residues = []

    species = []

    for rec in alignment:

        aa = rec.seq[col]
        all_residues.append(aa)

        if "OS=" in rec.description:
            species_name = rec.description.split("OS=")[1].split(" OX=")[0]
        else:
            species_name = rec.id

        species.append(f"{species_name}:{aa}")

    nongap = [aa for aa in all_residues if aa != "-"]

    if len(nongap) == 0:
        all_status = "Gap"
        identity = 0.0

    elif len(set(nongap)) == 1:
        all_status = "Identical"
        identity = 100.0

    else:

        unique = set(nongap)
        all_status = "Variable"

        for group in SIMILARITY_GROUPS:
            if unique.issubset(group):
                all_status = "Conserved"
                break

        identity = (
            100.0 * nongap.count(human_aa) / len(nongap)
        )

    return {
        "AA": human_aa,

        "Mouse_Status": mouse_status,
        "LT_Status": lt_status,
        "CT_Status": ct_status,

        "All_Status": all_status,
        "Identity": round(identity, 1),

        "SpeciesResidues": "; ".join(species),
    }

def collect_unique_sites(crosslink_folder):

    import glob
    import os

    unique_sites = set()

    csv_files = glob.glob(
        os.path.join(crosslink_folder, "*.csv")
    )

    if len(csv_files) == 0:
        raise FileNotFoundError(
            "No crosslink CSV files found"
        )

    print("\nCrosslink files found:")

    for csv in csv_files:

        print(os.path.basename(csv))

        df = pd.read_csv(csv)

        for _, row in df.iterrows():

            unique_sites.add(
                (
                    row["Protein1"],
                    int(row["Residue1"])
                )
            )

            unique_sites.add(
                (
                    row["Protein2"],
                    int(row["Residue2"])
                )
            )

    return unique_sites

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
        help="Folder containing crosslink CSV files"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="crosslink_conservation.csv"
    )

    args = parser.parse_args()
    
    # ---------------------------------------------------
    # Collect unique protein-residue sites from crosslinks
    # ---------------------------------------------------

    # Collect unique residues from ALL crosslink files

    unique_sites = collect_unique_sites(
        args.crosslinks)    

    summary = []
    for protein, residue in sorted(unique_sites):        
        result = analyse_position(protein, residue, args.alignments)

        summary.append({
        "Protein": protein,
        "Residue": residue,
        "Human_AA": result["AA"],

        "Mouse": result["Mouse_Status"],
        "LT": result["LT_Status"],
        "CT": result["CT_Status"],

        "All_Species": result["All_Status"],
        "Identity (%)": result["Identity"],

        "Species_Residues": result["SpeciesResidues"]
    })

    summary_df = pd.DataFrame(summary)
    summary_df = summary_df.sort_values(["Protein", "Residue"])

    summary_file = args.output.replace(".csv", "_protein_summary.csv")
    summary_df.to_csv(summary_file, index=False)

    print(f"Protein summary written to {summary_file}")

    for status_col in ["Mouse", "LT", "CT", "All_Species"]:

        protein_stats = (
            summary_df
            .groupby("Protein")[status_col]
            .value_counts()
            .unstack(fill_value=0)
        )

        # add missing columns
        for col in ["Identical", "Conserved", "Variable", "Gap", "Missing"]:
            if col not in protein_stats.columns:
                protein_stats[col] = 0

        protein_stats["Total Crosslinked Residues"] = protein_stats.sum(axis=1)

        protein_stats = protein_stats.reset_index()

        protein_stats.to_csv(
            args.output.replace(
                ".csv",
                f"_{status_col.lower()}_statistics.csv"
            ),
            index=False
        )

if __name__ == "__main__":
    main()