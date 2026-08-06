#!/usr/bin/env python3

import argparse
from Bio import AlignIO
from Bio.Align import substitution_matrices


# BLOSUM62
BLOSUM62 = substitution_matrices.load("BLOSUM62")


SPECIES = [
    "Homo sapiens",
    "Mus musculus",
    "Xenopus",
    "Danio",
    "Drosophila",
    "Caenorhabditis",
    "Chaetomium thermophilum",
    "Lachancea thermotolerans",
    "Saccharomyces",
]


def find_species(alignment, keyword):
    for rec in alignment:
        desc = rec.description
        if keyword.lower() in desc.lower():
            return rec
    return None


def similarity_score(ref_seq, target_seq):
    similar = 0
    compared = 0

    for a,b in zip(ref_seq,target_seq):

        if a == "-" or b == "-":
            continue

        compared += 1

        if a == b:
            similar += 1

        else:
            try:
                score = BLOSUM62[a,b]
                # Positive BLOSUM score = conservative substitution
                if score > 0:
                    similar += 1

            except KeyError:
                pass


    if compared == 0:
        return 0,0,0


    return (
        compared,
        100*similar/compared
    )



def main():

    parser = argparse.ArgumentParser(
        description="Calculate MSA similarity between species"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="FASTA alignment"
    )

    parser.add_argument(
        "-r",
        "--reference",
        default="Homo sapiens",
        help="Reference species"
    )

    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="Reference residue start"
    )

    parser.add_argument(
        "--end",
        type=int,
        required=True,
        help="Reference residue end"
    )


    args = parser.parse_args()


    alignment = AlignIO.read(
        args.input,
        "fasta"
    )


    ref_record = find_species(alignment, args.reference)

    if ref_record is None:
        raise ValueError(f"{args.reference} not found")

    print("\nReference:")
    print(ref_record.description)

    #
    # Convert residue numbering to alignment columns
    #

    start_col = None
    end_col = None

    count = 0

    for i,aa in enumerate(ref_record.seq):

        if aa != "-":
            count += 1

        if count == args.start:
            start_col = i

        if count == args.end:
            end_col = i+1
            break


    if start_col is None or end_col is None:
        raise ValueError(
            "Residue range not found"
        )


    ref_seq = ref_record.seq[start_col:end_col]


    print(
        f"\nRange: {args.start}-{args.end}"
    )


    print(
        "\n%-20s %12s %15s"
        %
        (
            "Species",
            "Compared",
            "Similarity %"
        )
    )

    print("-"*55)



    for species in SPECIES:
        if species == args.reference:
            continue

        record = find_species(alignment, species)

        if record is None:
            continue
        target_seq = record.seq[start_col:end_col]

        compared, similarity = similarity_score(
            ref_seq,
            target_seq
        )

        print(
            "%-30s %12d %14.1f"
            % (
                species,
                compared,
                similarity,
            )
        )



if __name__ == "__main__":
    main()