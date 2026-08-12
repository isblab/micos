#!/usr/bin/env python3

from Bio.Align import substitution_matrices

# Load BLOSUM45

BLOSUM45 = substitution_matrices.load("BLOSUM45")

# Paste your sequences here

# SEQ1 = """
# KALEHHRSEIQAEQDRKIEEVRDAMENEMRTQLRRQAAAHTDHLRDVLRVQEQELKSEFEQNLSEKLSEQELQFRRLSQEQVDNFTLDINTAYARLRGIEQAVQSHAVAEEEARKAHQLWLSVEALKYS-MKTSSAETPTIPLGSAVE
# """

# SEQ2 = """
# IRAKSREVELTQQFLNEFNAFKAQLEKHSSEELASALKANEQALLAKQSNEVALLSMKQVEEFTKILSEKLDQERQGRLSKLEALNGSVQELAEAVDQVDTLVMKSEVLSQLSLLTTLLKNKLHA--ESSVKIDSELARLKTLCDILP
# """

# SEQ1 = """
# TRAELIDRFRRVANEVRKASLL
# """

# SEQ2 = """
# SEETLRARFYAVQKLARRVAMI
# """

# SEQ1 = """
# GDDVESILARTQAFLEEGDLDNAAREMNALTGWSKTLSRDWLAEVRKVLEVRQALEVIQAEARLQSLRE
# """

# SEQ2 = """
# DINTFKLLSYASYCIEHGDLELAAKFVNQLKGESRRVAQDWLKEARMTLETKQIVEILTAYASAVGIGTT
# """

SEQ1 = """
LP----------E-S--VEKARSEVVRCLREHDRRPLNCWQEVEAFKEEVRKLEKGW
"""

SEQ2 = """
EVEAKFK---RYESHPVCADLQAKILQCYRENTHQTLKCSALATQYMHCVNHAKQSM
"""


def calculate_scores(seq1, seq2):
    similar = 0
    identical = 0
    compared = 0

    for a, b in zip(seq1, seq2):

        # Ignore positions containing gaps
        if a == "-" or b == "-":
            continue

        compared += 1

        # Exact identity
        if a == b:
            # print(a,b, compared)
            identical += 1
            similar += 1

        else:
            try:
                score = BLOSUM45[a, b]

                # Positive BLOSUM45 score = conservative substitution
                if score > 0:
                    similar += 1

            except KeyError:
                pass

    if compared == 0:
        return 0, 0, 0

    identity = 100 * identical / compared
    similarity = 100 * similar / compared

    return compared, identity, similarity

def main():

    seq1 = SEQ1.strip()
    seq2 = SEQ2.strip()

    compared, identity, similarity = calculate_scores(
        seq1,
        seq2
    )

    print(f"Compared:   {compared}")
    print(f"Identity:   {identity:.1f}%")
    print(f"Similarity: {similarity:.1f}%")

if __name__ == "__main__":
    main()
