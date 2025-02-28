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
from dataclasses import dataclass

import pysam
from loguru import logger


@dataclass
class FilterCauses:
    """ """

    unmapped: int = 0
    subs: int = 0
    internal_soft_clips: int = 0
    shorter_than_ref: int = 0
    ref_not_in_fasta: int = 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        Filter SAM alignments to retain only those reads that:
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
        "--output_sam",
        required=True,
        help="Output filtered SAM file",
    )
    return parser.parse_args()


def load_reference_lengths(reference_fasta: str) -> dict[str, int]:
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

    # Count how many soft clipping operations we find
    soft_clip_op = 4
    soft_clip_positions = [idx for idx, (op, _) in enumerate(cigar_tuples) if op == soft_clip_op]

    # If no soft clipping at all, return True
    if not soft_clip_positions:
        return True

    # Check if all soft clipping positions are at the ends
    return all(pos == 0 or pos == len(cigar_tuples) - 1 for pos in soft_clip_positions)


def has_no_substitutions(cigar_tuples: list[tuple[int, int]]) -> bool:
    """
    Check that there are no substitution operations ('X') in the CIGAR string.

    Args:
        cigar_tuples (list of tuples): CIGAR operations as (operation, length).

    Returns:
        bool: True if no 'X' operations are present, False otherwise.
    """
    # return whether all CIGAR operators do *not* equal the substitution operator, 8.
    substitution_code = 8
    return all(op != substitution_code for op, _ in cigar_tuples)


def log_filtering_stats(
    read_count: int,
    retained_tally: int,
    filter_causes: FilterCauses,
    reference_name: str,
) -> None:
    """
    Log statistics about the filtering process.

    Args:
        read_count (int): Total number of reads processed.
        retained_tally (int): Number of reads retained after filtering.
        filter_causes (FilterCauses): Object containing counts of different filter causes.
        reference_name (str): Path to reference FASTA file.
    """
    # create a tally of the number of reads filtered
    filtered_tally = read_count - retained_tally

    # log out the results
    logger.info(
        f"Filtered {filtered_tally} reads from the input reads, leaving {retained_tally} behind. In all, {(filtered_tally / (filtered_tally + retained_tally)) * 100}% of reads were filtered and will thus not be included in the final genotyping report.",
    )
    logger.info("Proceeding to filtering statistics:")
    logger.info(
        f"{filter_causes.unmapped} reads were filtered because they were unmapped.",
    )
    logger.info(
        f"{filter_causes.subs} were filtered because they had substitutions relative to a reference.",
    )
    logger.info(
        f"{filter_causes.internal_soft_clips} were filtered because they had soft clips internal to the alignment.",
    )
    logger.info(
        f"{filter_causes.ref_not_in_fasta} reads were filtered because the reference they were mapped to does not have an identifier that matches a reference in {reference_name}.",
    )
    logger.info(
        f"{filter_causes.shorter_than_ref} were filtered because they are shorter than the reference.",
    )


def filter_alignments(input_sam: str, reference_fasta: str, output_sam: str) -> None:  # noqa: C901
    """
    Filter alignments based on the specified criteria and write to output SAM.

    Args:
        input_sam (str): Path to the input SAM file.
        reference_fasta (str): Path to the reference FASTA file.
        output_sam (str): Path to the output filtered SAM file.
    """
    ref_lengths = load_reference_lengths(reference_fasta)

    # Open input BAM
    try:
        samfile = pysam.AlignmentFile(input_sam, "r")
    except FileNotFoundError:
        sys.stderr.write(f"Error: Input SAM file '{input_sam}' not found.\n")
        sys.exit(1)

    # Open output SAM (copy header from input)
    outfile = pysam.AlignmentFile(output_sam, "w", header=samfile.header)

    # initialize a dictionary to store the reasons for each filter
    filtered_tallies = FilterCauses()
    read_count = 0
    retained_tally = 0
    for read in samfile.fetch(until_eof=True):
        # incrememt read count
        read_count += 1

        # Flag to track if the read passes all filters
        is_read_retained = True

        # Skip unmapped reads
        if read.is_unmapped:
            filtered_tallies.unmapped += 1
            continue

        # Check for soft clipping only at the ends in CIGAR
        if read.cigartuples is None:
            filtered_tallies.internal_soft_clips += 1
            continue
        if not is_soft_clipping_only_at_ends(read.cigartuples):
            filtered_tallies.internal_soft_clips += 1
            continue

        # Ensure no substitutions ('X' operations) in CIGAR
        if not has_no_substitutions(read.cigartuples):
            filtered_tallies.subs += 1
            continue

        # Check if alignment starts at position 1 (0-based)
        if read.reference_start != 0:
            filtered_tallies.shorter_than_ref += 1
            continue

        # Calculate alignment length (sum of M, =, X operations)
        aln_length = 0
        for op, length in read.cigartuples:
            if op not in [0, 7, 8]:  # M, =, X
                continue
            aln_length += length

        # Get reference length
        ref_name = read.reference_name
        ref_length = ref_lengths.get(str(ref_name))
        if ref_length is None:
            # Reference name not found in reference lengths
            filtered_tallies.ref_not_in_fasta += 1
            is_read_retained = False
            continue

        # Check if alignment spans the entire reference
        if aln_length < ref_length:
            filtered_tallies.shorter_than_ref += 1
            is_read_retained = False
            continue

        if is_read_retained:
            retained_tally += 1
            outfile.write(read)

    # Close files
    samfile.close()
    outfile.close()

    # Verify that all reads are accounted for
    assert read_count == (
        filtered_tallies.unmapped
        + filtered_tallies.internal_soft_clips
        + filtered_tallies.subs
        + filtered_tallies.shorter_than_ref
        + filtered_tallies.ref_not_in_fasta
        + retained_tally
    ), "Total reads do not match filtered and retained counts"

    # log out how much filtering occurred along with information on why filtering occurred
    log_filtering_stats(
        read_count,
        retained_tally,
        filtered_tallies,
        reference_fasta,
    )


def main() -> None:
    args = parse_arguments()
    filter_alignments(args.input_sam, args.reference_fasta, args.output_sam)


if __name__ == "__main__":
    main()
