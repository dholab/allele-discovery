#!/usr/bin/env python3

from __future__ import annotations

import argparse
from typing import TYPE_CHECKING

from Bio import SeqIO

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


def calculate_avg_quality(record: SeqRecord) -> int | float:
    """Calculate the average quality score for a single FASTQ record."""
    qualities = record.letter_annotations["phred_quality"]
    return sum(qualities) / len(qualities)


def get_top_reads(
    input_fastq: str,
    output_fastq: str,
    log_file: str,
    top_n: int,
) -> None:
    """Extract the top N reads with the highest average quality from a FASTQ file."""

    # Parse all the FASTQ records and calculate their average quality
    records = list(SeqIO.parse(input_fastq, "fastq"))

    if len(records) == 0:
        msg = f"No records found in {input_fastq}"
        raise Exception(msg)  # noqa: TRY002

    records_with_avg_quality = [
        (record, calculate_avg_quality(record)) for record in records
    ]

    # Sort the records by average quality score in descending order
    sorted_records = sorted(records_with_avg_quality, key=lambda x: x[1], reverse=True)

    # Get the top N records with the highest quality
    top_records = [record for record, _ in sorted_records[:top_n]]

    # Write the top N records to the output FASTQ file
    with open(output_fastq, "w") as output_handle:
        SeqIO.write(top_records, output_handle, "fastq")

    # Log the results
    with open(log_file, "w") as log_handle:
        log_handle.write(
            f"Extracted top {top_n} reads with the highest average quality from {input_fastq}.\n",
        )
        log_handle.write(f"Total records processed: {len(records)}\n")
        log_handle.write(f"Top {top_n} records written to {output_fastq}.\n")


def main() -> None:
    # Argument parser for command-line inputs
    parser = argparse.ArgumentParser(description="Extract top quality FASTQ reads.")
    parser.add_argument("--input", type=str, required=True, help="Input FASTQ file")
    parser.add_argument("--output", type=str, required=True, help="Output FASTQ file")
    parser.add_argument(
        "--log",
        type=str,
        required=True,
        help="Log file to write details",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        required=True,
        help="Number of top reads to extract",
    )

    args = parser.parse_args()

    # Run the extraction
    get_top_reads(args.input, args.output, args.log, args.top_n)


if __name__ == "__main__":
    main()
