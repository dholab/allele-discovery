#!/usr/bin/env python3

import argparse
import csv
import os

import pandas as pd
import utils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a FASTA file to use for genotyping an input set of reads.",
    )
    parser.add_argument(
        "-n",
        "--novel_alleles",
        type=str,
        required=True,
        help="Input FASTA file of putative novel allele sequences.",
    )
    parser.add_argument(
        "-r",
        "--known_alleles",
        type=str,
        required=True,
        help="Reference FASTA of previously known alleles.",
    )
    parser.add_argument(
        "-a",
        "--alignment",
        type=str,
        required=True,
        help="Alignment of reads to previously known alleles.",
    )
    parser.add_argument(
        "-c",
        "--closest_matches",
        type=str,
        required=True,
        help="Matrix of closest matches between reads and previously known alleles.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_command_line_args()

    # create dictionary of classifications
    classified = {}

    # add cdna matches
    with open(args.alignment) as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
    for row in reader:
        classified[row[0]] = ["extend-cdna", utils.removeSpecialCharacters(row[1])]

    # add novel matches
    if os.stat(args.closest_matches).st_size > 0:  # noqa: PTH116
        novel_df = pd.read_excel(args.closest_matches, index_col=0)

    for _, row in novel_df.iterrows():
        classified[str(row[0])] = [
            "novel",
            utils.removeSpecialCharacters(
                row[1] + "|" + row[2] + "|" + row[3] + "|" + row[4] + "|" + row[5],
            ),
        ]

    # create renamed FASTA file with updated names for genotyping
    with open(args.novel_alleles, "a") as handle:
        # add IPD gDNA sequences
        with open(args.known_alleles) as input_handle:
            sequences = SeqIO.parse(input_handle, "fasta")
            SeqIO.write(sequences, handle, "fasta")

        # concatenate cDNA extensions and novel sequences with known IPD gDNA sequences
        # this enables genotyping against an expanded gDNA library even when there aren't a huge number of gDNA matches in this specific set of samples
        for record in SeqIO.parse(args.novel_alleles, "fasta"):
            # get information for matching sequence
            allele_data = classified.get(record.name)

            if allele_data is not None:
                renamed_name = str(
                    allele_data[1] + "_" + allele_data[0] + "_" + record.name,
                )
                renamed_seq = record.seq
                new_record = SeqRecord(
                    renamed_seq,
                    id=renamed_name,
                    description="",
                )

                SeqIO.write(new_record, handle, "fasta")


if __name__ == "__main__":
    main()
