#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pysam",
# ]
# ///

from __future__ import annotations

import argparse
import sys

import pysam


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        Filter BAM alignments to retain only those reads that:
        - Span the entire reference sequence (start at position 1 and cover full reference length).
        - Have no substitutions (no 'X' operations in CIGAR).
        - Allow indels (insertions and deletions) relative to the reference.
        - Contain soft clipping only at the ends of the CIGAR string.
        """,
    )
    parser.add_argument(
        "-i",
        "--input_sam",
        required=True,
        help="Input SAM file (unsorted)",
    )
    parser.add_argument(
        "-r",
        "--reference_fasta",
        required=True,
        help="Reference FASTA file (indexed with samtools faidx)",
    )
    parser.add_argument(
        "-o",
        "--output_bam",
        required=True,
        help="Output filtered BAM file",
    )
    return parser.parse_args()


def load_reference_lengths(reference_fasta: str) -> dict[str, str]:
    """
    Load reference sequence lengths from the .fai index file.

    Args:
        reference_fasta (str): Path to the reference FASTA file.

    Returns:
        dict: Mapping from reference name to its length.
    """
    fai_path = reference_fasta + ".fai"
    ref_lengths = {}
    try:
        with open(fai_path) as fai:
            for line in fai:
                parts = line.strip().split("\t")
                ref_name = parts[0]
                ref_length = int(parts[1])
                ref_lengths[ref_name] = ref_length
    except FileNotFoundError:
        sys.stderr.write(
            f"Error: FAI index file not found for {reference_fasta}. Please run 'samtools faidx {reference_fasta}' first.\n",
        )
        sys.exit(1)
    return ref_lengths


def is_soft_clipping_only_at_ends(cigar_tuples: list[tuple[int, int]]) -> bool:
    """
    Check if soft clipping occurs only at the ends of the CIGAR string.

    Args:
        cigar_tuples (list of tuples): CIGAR operations as (operation, length).

    Returns:
        bool: True if soft clipping is only at the ends, False otherwise.
    """
    if not cigar_tuples:
        return False

    num_ops = len(cigar_tuples)
    soft_clip_op = 4
    for idx, (op, _) in enumerate(cigar_tuples):
        if op == soft_clip_op and idx != 0 and idx != num_ops - 1:  # Soft clipping
            return False
    return True


def has_no_substitutions(cigar_tuples: list[tuple[int, int]]) -> bool:
    """
    Check that there are no substitution operations ('X') in the CIGAR string.

    Args:
        cigar_tuples (list of tuples): CIGAR operations as (operation, length).

    Returns:
        bool: True if no 'X' operations are present, False otherwise.
    """
    substitution_code = 8
    return all(op != substitution_code for op, _ in cigar_tuples)


def filter_alignments(input_sam: str, reference_fasta: str, output_bam: str) -> None:  # noqa: C901
    """
    Filter alignments based on the specified criteria and write to output BAM.

    Args:
        input_bam (str): Path to the input BAM file.
        reference_fasta (str): Path to the reference FASTA file.
        output_bam (str): Path to the output filtered BAM file.
    """
    ref_lengths = load_reference_lengths(reference_fasta)

    # Open input BAM
    try:
        samfile = pysam.AlignmentFile(input_sam, "r")
    except FileNotFoundError:
        sys.stderr.write(f"Error: Input BAM file '{input_sam}' not found.\n")
        sys.exit(1)

    # Open output BAM (copy header from input)
    outfile = pysam.AlignmentFile(output_bam, "wb", header=samfile.header)

    for read in samfile.fetch(until_eof=True):
        # Skip unmapped reads
        if read.is_unmapped:
            continue

        # Check for soft clipping only at the ends in CIGAR
        if read.cigartuples is None:
            continue

        if not is_soft_clipping_only_at_ends(read.cigartuples):
            continue

        # Ensure no substitutions ('X' operations) in CIGAR
        if not has_no_substitutions(read.cigartuples):
            continue

        # Check if alignment starts at position 1 (0-based)
        if read.reference_start != 0:
            continue

        # Calculate alignment length (sum of M, =, X operations)
        aln_length = 0
        for op, length in read.cigartuples:
            if op in [0, 7, 8]:  # M, =, X
                aln_length += length

        # Get reference length
        ref_name = read.reference_name
        ref_length = ref_lengths.get(str(ref_name))
        if ref_length is None:
            # Reference name not found in reference lengths
            continue

        # Check if alignment spans the entire reference
        if aln_length == ref_length:
            outfile.write(read)

    # Close files
    samfile.close()
    outfile.close()


def main() -> None:
    args = parse_arguments()
    filter_alignments(args.input_sam, args.reference_fasta, args.output_bam)


if __name__ == "__main__":
    main()
