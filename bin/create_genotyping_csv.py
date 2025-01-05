#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path

OUTPUT_NAME = sys.argv[1] if len(sys.argv) > 1 else "genotyping.csv"

# Make list from samples
sam_files = Path.cwd().glob("*.sam")


# Define a function to sanitize genotype names
def sanitize_genotype(genotype: str) -> str:
    """
    Replace all non-alphanumeric characters in the genotype string with underscores.

    Parameters:
    - genotype (str): The original genotype string.

    Returns:
    - str: The sanitized genotype string.
    """
    # Replace any character that is not a letter, digit, or underscore with '_'
    return re.sub(r"[^\w]", "_", genotype)


# Create or overwrite the CSV file
with open(OUTPUT_NAME, "w", newline="") as genotyping_csv:
    writer = csv.writer(genotyping_csv)

    # Optionally, write a header row
    writer.writerow(["animal", "genotype"])

    for sam_file in sam_files:
        # Get sample name by removing .sam and any other extensions. This method assumes
        # that the sample name is the first item in a period-delimited file name.
        sam_file_basename = Path(sam_file).name.split(".")[0]
        animal_name = re.sub(r"\.sam$", "", sam_file_basename)

        with open(sam_file) as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                if row[0].startswith("@"):
                    continue

                # Ensure that the row has at least 3 columns
                valid_column_count = 3
                if len(row) < valid_column_count:
                    continue  # Skip malformed rows

                # Get genotype from SAM file (third column, index 2)
                genotype = row[2]

                # Sanitize the genotype string to replace special characters with '_'
                sanitized_genotype = sanitize_genotype(genotype)

                # Write to output CSV
                writer.writerow([animal_name, sanitized_genotype])
