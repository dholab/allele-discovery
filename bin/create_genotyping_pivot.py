#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
#     "polars",
#     "pysam",
# ]
# ///

"""

A module for parsing SAM/BAM files to extract and count genotype information from
sequence reads.

This module provides functions for processing SAM/BAM files to parse out sample
identification and genotyping information embedded in sequence read headers. It handles
converting from SAM/BAM formats into a tabular representation suitable for genotype
analysis.

Key functions:
    - parse_record_id(): Parse sample and read count data from FASTA record headers
    - extract_ref_query_pairs(): Extract reference-query pairs from SAM/BAM files
    - count_genotypes(): Count genotype hits per sample from SAM/BAM files
    - pivot_genotypes(): Convert long-format genotype data into a wide sample matrix
"""

from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pysam
from loguru import logger

READ_COUNT_DELIM: str = "_ReadCount-"
CLUSTER_DELIM: str = "_cluster-"


@dataclass
class SampleRefCount:
    """A dataclass for holding sample identification metadata and sequence read metadata.

    Args:
        sample_id: The identifier for a sequencing sample
        ref_id: The identifier for a reference sequence aligning to reads from this sample
        read_count: The number of reads from this sample aligning to the reference
    """

    sample_id: str
    ref_id: str
    read_count: int


def clean_sample_name(putative_id: str) -> str:
    """Clean up the raw ID of a putative sample name.

    Args:
        putative_id: The raw identifier string to be cleaned

    Returns:
        str: The cleaned sample identifier with file extensions removed
    """
    return putative_id.split(".")[0]


def parse_record_id(record_id: str, ref_id: str) -> SampleRefCount:
    """Parse a read-identifying record from an input FASTA file.

    Args:
        record_id: The read identifier string from a FASTA file to parse
        ref_id: The identifier of the reference sequence that this read aligns to

    Returns:
        SampleRefCount: A dataclass containing parsed metadata about the sample id, reference id,
                       and read count for this identifier

    Raises:
        AssertionError: If the record ID is malformed and cannot be properly parsed
        ValueError: If the read count cannot be cast as an integer
    """
    # make sure the expected read count delimiter is present in the input record
    assert (
        READ_COUNT_DELIM in record_id
    ), f"The provided record ID, {record_id} does not contain an '{READ_COUNT_DELIM}' delimiter, and thus can't be parsed properly. Aborting."

    # split out the final element after the read count delimiter, which *should* be the read count
    read_count_split = record_id.split(READ_COUNT_DELIM)
    expected_splits = 2
    assert (
        len(read_count_split) == expected_splits
    ), f"Splitting of {record_id} on {READ_COUNT_DELIM} yielded more than the expected two results, indicating that the cluster FASTA is corrupted:\n{read_count_split}"

    # make sure the read count is an integer
    try:
        int(read_count_split[1])
    except ValueError as error:
        logger.error(
            f"Unable to cast {read_count_split[1]} as an integer, indicating that the input FASTA does not have the expected format: {error}",
        )
    read_count = int(read_count_split[0])

    # parse sample Id from the record
    underscore_split = record_id.split("_")
    assert (
        len(underscore_split) > 1
    ), f"Unable to properly split the Record ID {record_id} by underscores to find the sample ID: {underscore_split}"
    sample_id = clean_sample_name(underscore_split[0])

    return SampleRefCount(
        sample_id=sample_id,
        ref_id=ref_id,
        read_count=read_count,
    )


def flatten_counts(
    sample_ref_lookup: dict[str, dict[str, int]],
) -> list[SampleRefCount]:
    """Flatten a deeply nested dictionary of read counts into a list of sample & reference metadata.

    Args:
        sample_ref_lookup: A dictionary mapping sample IDs to an inner dictionary of reference IDs and read counts

    Returns:
        list[SampleRefCount]: A flattened list of SampleRefCount objects containing the sample/reference metadata
    """
    flat_lookups: list[SampleRefCount] = []
    for sample_id, ref_lookup in sample_ref_lookup.items():
        for ref_id, read_count in ref_lookup.items():
            flat_lookup = SampleRefCount(
                sample_id,
                ref_id,
                read_count,
            )
            flat_lookups.append(flat_lookup)

    return flat_lookups


def extract_ref_query_pairs(sam_file: str) -> list[SampleRefCount]:
    """Extract paired reference and query sequences from a SAM/BAM file.

    Args:
        sam_file: Path to a SAM/BAM file containing aligned reads

    Returns:
        list[SampleRefCount]: A list of SampleRefCount objects containing the sample ID,
                             reference sequence name, and number of supporting reads for each alignment

    Note:
        This function skips unmapped reads and alignments with missing query or reference names.
        It aggregates read counts per sample and reference combination.
    """
    # initialize a dictionary where the first key is the sample id, the
    # second key is the reference name, and the innermost value is the
    # number of reads supporting that call
    sample_ref_lookup: dict[str, dict[str, int]] = {}

    # Open SAM/BAM file
    with pysam.AlignmentFile(sam_file, "r") as sam:
        # Iterate through alignments
        for read in sam:
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Get query sequence
            query_name = read.query_name

            # make sure a string was parsed
            if query_name is None:
                continue

            # get reference sequence
            ref_name = read.reference_name

            # make sure a string was parsed
            if ref_name is None:
                continue

            # parse out metadata
            record_meta = parse_record_id(query_name, ref_name)

            # add to the growing tally
            if record_meta.sample_id not in sample_ref_lookup:
                sample_ref_lookup[record_meta.sample_id] = {
                    ref_name: record_meta.read_count,
                }
                continue

            if ref_name not in sample_ref_lookup[record_meta.sample_id]:
                sample_ref_lookup[record_meta.sample_id][ref_name] = record_meta.read_count
                continue

            sample_ref_lookup[record_meta.sample_id][ref_name] += record_meta.read_count

    return flatten_counts(sample_ref_lookup)


def count_genotypes(sam_files: list[str]) -> pl.DataFrame:
    """Count the number of alignments for each sample-by-reference combination.

    Args:
        sam_files: A list of paths to SAM/BAM alignment files to count hits from

    Returns:
        pl.DataFrame: A dataframe containing the number of reads supporting each sample/reference combination.
                     Columns are 'sample' (sample ID), 'genotype' (reference sequence), and 'support' (read count).
    """
    per_sample_dfs: list[pl.DataFrame] = []
    for sam in sam_files:
        lookup = extract_ref_query_pairs(sam)
        sample_ids = [count.sample_id for count in lookup]
        ref_hits = [count.ref_id for count in lookup]
        counts = [count.read_count for count in lookup]

        if not len(sample_ids) == len(ref_hits) == len(counts):
            logger.error(
                f"Failed to count hits for each sample-reference combination in {sam}. All of the following lists should be the same length:\n{sample_ids}\n{ref_hits}\n{counts}",
            )
            assert len(sample_ids) == len(ref_hits) == len(counts)

        per_sample_dfs.append(
            pl.DataFrame(
                {
                    "sample": sample_ids,
                    "genotype": ref_hits,
                    "support": counts,
                },
            ),
        )

    return pl.concat(per_sample_dfs)


def pivot_genotypes(long_df: pl.DataFrame) -> pl.DataFrame:
    """Pivot a long genotype-by-sample dataframe into an even wider format by sample.

    Args:
        long_df: A polars DataFrame in long format with 'sample', 'genotype' and 'support' columns

    Returns:
        pl.DataFrame: A wider format DataFrame with genotypes as rows and samples as columns
    """
    return (
        long_df.group_by(["sample", "genotype"])  # noqa: PD010
        .agg(pl.count("support").alias("support"))
        .pivot(on="sample", index="genotype", values="support")
    )


def main() -> None:
    sam_files = [str(path) for path in Path.cwd().glob("*.sam")]
    long_df = count_genotypes(sam_files)
    pivot_table = pivot_genotypes(long_df)
    pivot_table.write_excel("genotyping_report.xlsx")


if __name__ == "__main__":
    main()
