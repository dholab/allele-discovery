#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython",
# ]
# ///

import csv
import sys

from Bio import SeqIO


def main() -> None:
    assert (
        len(sys.argv) > 1
    ), "Usage: python3 rename_cdna_matched.py <cdna matches fasta> <matches tsv>"

    # pull the paths to the required input path from the command line
    cdna_match_fasta = sys.argv[1]
    matches_tsv_path = sys.argv[2]

    # create an empty output file to append to
    with open("cdna_matches.fasta", "w") as _:
        pass

    # read input FASTA line-by-line
    for record in SeqIO.parse(cdna_match_fasta, "fasta"):
        # parse file with gdna sequences that match cdna sequences
        with open(matches_tsv_path) as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                # test if name of sequence in cdna match file matches gdna sequence name
                if row[0] == record.name:
                    # update name of sequence in output file
                    record.description = row[1] + "|" + row[0]

                    # write to file
                    with open("cdna_matches.fasta", "w") as handle:
                        SeqIO.write(record, handle, "fasta")


if __name__ == "__main__":
    main()
