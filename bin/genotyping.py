#!/usr/bin/env python3

import sys

import pandas as pd

# Check if input file is provided
if len(sys.argv) < 2:  # noqa: PLR2004
    print("Usage: python genotyping.py <input_csv> <output_xlsx_name>")  # noqa: T201
    sys.exit(1)


genotype_df = pd.read_csv(sys.argv[1], sep=",", names=["animal", "genotype"])

# add column to hold counts
genotype_df["ct"] = 1

# group by number of times a genotype appears in an animal
df_grouped = genotype_df.groupby(["animal", "genotype"])["ct"].count().reset_index()

# create pivot table
df_pivoted = pd.pivot_table(
    df_grouped,
    index=["genotype"],
    columns=["animal"],
    values=["ct"],
)

# export to Excel
# has index column that should be deleted to avoid confusing numbering
df_pivoted.to_excel(sys.argv[2])
