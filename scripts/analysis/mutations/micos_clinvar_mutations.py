import os
from string import Template
from typing import Any
import xmltodict
import yaml
import pprint
import warnings
import argparse
import pandas as pd
from tqdm import tqdm
from datetime import datetime
import xml.etree.ElementTree as ET
from IMP_Toolbox.sequence.sequence_alignment import (
    PairwiseSequenceAlignment,
)

from IMP_Toolbox.utils.file_helpers import (
    read_fasta,
)

from IMP_Toolbox.mutations.af_missense import (
    af_missense_df_to_dict,
    fetch_fasta_dict_for_af_missense,
    fetch_af_missense_data,
    export_af_missense_data,
    get_af_missense_attribute,
)
from IMP_Toolbox.constants.mutation_constants import (
    AF_MISSENSE_COLUMNS,
    AF_MISSENSE_CSV_SUFFIX,
    AF_MISSENSE_PAIR_ALN_SUFFIX,
    AF_MISSENSE_AA_SUBSTITUTIONS_TSV,
    CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE,
    CLINVAR_DF_COLUMNS,
    CLINVAR_TEMPLATE_QUERY_DETAIL,
    CLINVAR_TEMPLATE_QUERY_ID,
    DATE_FORMAT,
    API_URLS,
)

from IMP_Toolbox.mutations.clinvar_mutations import (
    VariantInfo,
    process_clinvar_variant_data,



)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Fetch missense variants from ClinVar."
    )
    parser.add_argument(
        "--config_file",
        type=str,
        default="/home/muskaan/projects/micos/github/inputs/mutations/configs.yml",
        help="Path to the configuration YAML file.",
    )
    parser.add_argument(
        "--include_VUS",
        action="store_true",
        default=False,
        help="Include variants of uncertain significance (VUS).",
    )
    parser.add_argument(
        "--clinvar_output_dir",
        type=str,
        default="/home/muskaan/projects/micos/github/results/mutations/clinvar",
        help="Directory to save ClinVar variant information.",
    )
    parser.add_argument(
        "--alpha_missense_dir",
        type=str,
        default="/home/muskaan/projects/micos/github/results/mutations/alpha_missense",
        help="Directory to save AlphaMissense variant information.",
    )
    parser.add_argument(
        "--pairwise_alignments_dir",
        type=str,
        default="/home/muskaan/projects/micos/github/results/sequence_alignments/pairwise_alignments",
        help="Directory to save pairwise alignments.",
    )
    parser.add_argument(
        "--include_AF_missense",
        action="store_true",
        default=True,
        help="Include AlphaMissense pathogenicity scores.",
    )
    parser.add_argument(
        "--af_missense_mode",
        type=str,
        choices=["online", "offline"],
        default="online", ## TSV file is not downloaded automatically. 
        help="Mode to fetch AlphaMissense data.",
    )
    parser.add_argument(
        "--af_missense_tsv",
        type=str,
        required=False,
        default=AF_MISSENSE_AA_SUBSTITUTIONS_TSV,
        help="Path to AlphaMissense aa substitutions TSV file for offline mode.",
    )
    parser.add_argument(
        "--protein_sequences_fasta",
        type=str,
        default="/home/muskaan/projects/micos/github/inputs/mutations/proteins.fasta",
        # default=protein_sequences_fasta,
        help="Fasta file containing modeled protein sequences.",
    )
    args = parser.parse_args()

    config_yaml = yaml.load(open(args.config_file, "r"), Loader=yaml.FullLoader)
    protein_uniprot_map = config_yaml["micos_protein_uniprot_map"]
    protein_gene_map = config_yaml["micos_protein_gene_map"]
    protein_sequences = read_fasta(args.protein_sequences_fasta)

    if args.include_VUS:
        CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE.append("Uncertain significance")

    os.makedirs(args.clinvar_output_dir, exist_ok=True)

    df_rows = []
    fasta_dict = fetch_fasta_dict_for_af_missense(
        os.path.join(args.alpha_missense_dir, "af_missense_sequences.fasta"),
        protein_uniprot_map,
    )
    uniprot_bases = [uid.split("-")[0] for uid in protein_uniprot_map.values()]

    #######################################################################
    # Fetch AlphaMissense data
    #######################################################################
     #######################################################################
    # Fetch AlphaMissense data
    #######################################################################
    if args.include_AF_missense:

        af_missense_df_gen = fetch_af_missense_data(
            args.alpha_missense_dir,
            uniprot_bases,
            mode=args.af_missense_mode,
            overwrite=False,
            af_missense_tsv=args.af_missense_tsv,
        )

        export_af_missense_data(
            args.alpha_missense_dir,
            af_missense_df_gen,
            overwrite=False,
        )

    df = process_clinvar_variant_data(
        protein_uniprot_map=protein_uniprot_map,
        protein_gene_map=protein_gene_map,
        protein_sequences=protein_sequences,
        clinvar_output_dir=args.clinvar_output_dir,
        pairwise_alignments_dir=args.pairwise_alignments_dir,
        include_AF_missense=args.include_AF_missense,
        include_VUS=args.include_VUS,
        alpha_missense_dir=args.alpha_missense_dir,
    )

    # print(df.head())

    out_name = "clinvar_missense_variants"
    if args.include_VUS:
        out_name += "_with_VUS"

    df_file = os.path.join(args.clinvar_output_dir, f"{out_name}.xlsx")
    df.to_excel(df_file, index=False)

    print(f"ClinVar missense variants saved to {df_file}")