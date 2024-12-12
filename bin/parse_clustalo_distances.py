#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython",
#     "pandas",
# ]
# ///
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="",
    )
    parser.add_argument(
        "-d",
        "--distances",
        type=str,
        required=True,
        help="Text file of nucleotide distances.",
    )
    parser.add_argument(
        "-n",
        "--novel_alleles",
        type=str,
        required=True,
        help="Input FASTA file of putative novel allele sequences.",
    )
    parser.add_argument(
        "-c",
        "--cdna_matches",
        type=str,
        required=True,
        help="Input FASTA file of any identified cDNA matches.",
    )
    parser.add_argument(
        "-m",
        "--match_excel",
        type=str,
        required=False,
        default="novel_closest_matches.xlsx",
        help="Excel-formatted table of closest matches for each novel allele.",
    )
    parser.add_argument(
        "-t",
        "--distances_tmp",
        type=str,
        required=False,
        default="distances_tmp.txt",
        help="Temporary text file for intermediate distance computations.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_command_line_args()

    # test if there is data in novel.fasta
    # create output files if empty
    if Path.stat(args.novel_alleles).st_size == 0:
        Path(args.match_excel).touch()
        Path(args.distances_tmp).touch()
        sys.exit()

    # import clustal omega distance file
    # remove spaces in sequence identifiers that makes distances.txt file human readable
    with open(args.distances) as fin, open(args.distances_tmp, "w") as fout:
        for line in fin:
            fout.write(re.sub(" +", " ", line))

    # import into pandas dataframe
    clustal_df = pd.read_csv(
        args.distances_tmp,
        skiprows=1,
        header=None,
        sep=" ",
        index_col=0,
    )

    # use row_index names as column names
    ordered_values = clustal_df.index.tolist()
    clustal_df.columns = ordered_values
    clustal_df.columns = clustal_df.columns.astype(str)

    # convert column names that are cdna_matches to cdna names
    # this makes it possible to glean lineages for novel alleles that are closest matches to cdna extensions
    cdna_dict = {
        seq_record.name: seq_record.description.split(" ")[1]
        for seq_record in SeqIO.parse(args.cdna_matches, "fasta")
    }
    clustal_df = clustal_df.rename(cdna_dict, axis=1)

    # self-by-self comparisons always have a distance of 0.000000
    # for novel alleles, there will never be a perfect match to an existing sequence
    # so convert 0.000000	to 1 so self-matches don't report in nsmallest calculation
    clustal_df = clustal_df.replace(0.000000, 1)

    # get names of sequences in novel.fasta
    identifiers = [
        seq_record.id for seq_record in SeqIO.parse(args.novel_alleles, "fasta")
    ]

    # create subset df with only novel fasta records as rows
    novel_df = clustal_df.loc[identifiers, :]

    # remove columns that are also novel alleles since these are not informative for naming
    novel_df = novel_df.drop(columns=identifiers)

    # filter on best three matches
    # sort on rows with lowest values (closest matches)
    # from https://stackoverflow.com/questions/54923349/top-3-values-per-row-in-pandas
    c = ["Self", "Closest", "2", "3", "4", "5"]
    ranked_df = novel_df.apply(
        lambda x: pd.Series(x.nsmallest(6).index, index=c),
        axis=1,
    ).reset_index()

    # drop 'self' match which will always be the best match
    top5_df = ranked_df.drop(columns=["Self"])
    top5_df.columns = top5_df.columns.astype(str)

    # Rename the original index column (previously 'index' -> 'Query')
    # Assuming the first column in top5_df is the reset index from ranked_df
    old_index_col = top5_df.columns[0]
    top5_df = top5_df.rename(columns={old_index_col: "Query"})

    match_cols = ["Closest", "2", "3", "4", "5"]

    def append_distance(row):
        query = row["Query"]
        for col in match_cols:
            match_name = row[col]
            dist_val = novel_df.loc[query, match_name].round(3)
            row[col] = f"{match_name} ({dist_val})"
        return row

    top5_df = top5_df.apply(append_distance, axis=1)

    top5_df.to_excel(args.match_excel)


if __name__ == "__main__":
    main()
