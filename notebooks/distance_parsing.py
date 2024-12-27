import marimo

__generated_with = "0.8.22"
app = marimo.App(width="medium")


@app.cell
def __():
    import polars as pl
    import re
    from Bio import SeqIO
    return SeqIO, pl, re


@app.cell
def __():
    # normalize how elements of each line are spaced and convert to a proper TSV file
    with open("scratch/distances.txt") as fin, open(
        "scratch/distances_tmp.tsv", "w"
    ) as fout:
        for line in fin:
            cleaned_line = "\t".join(line.split())
            cleaned_line = cleaned_line.replace("0.000000", "1")
            fout.write(f"{cleaned_line}\n")
    return cleaned_line, fin, fout, line


@app.cell
def __(pl):
    # use the first column of the distance matrix to construct a list of column names
    new_column_names = ["alleles"] + pl.scan_csv(
        "scratch/distances_tmp.tsv",
        separator="\t",
        skip_rows=1,
        has_header=False,
    ).select("column_1").collect().to_series().to_list()
    return (new_column_names,)


@app.cell
def __(new_column_names, pl):
    distances = pl.scan_csv(
        "scratch/distances_tmp.tsv",
        separator="\t",
        skip_rows=1,
        has_header=False,
        new_columns=new_column_names,
    )
    return (distances,)


@app.cell
def __(SeqIO):
    identifiers = [
        seq_record.description
        for seq_record in SeqIO.parse("scratch/novel.validated.fasta", "fasta")
    ]
    identifiers
    return (identifiers,)


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
