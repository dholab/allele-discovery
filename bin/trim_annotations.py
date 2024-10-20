import pandas as pd
from Bio import SeqIO

# import GFF file to dataframe
df = pd.read_csv(
    snakemake.input[0],
    sep="\t",
    names=[
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ],
)

# read sequence from FASTA file
with open(snakemake.output[0], "w") as handle:
    for record in SeqIO.parse(snakemake.input[1], "fasta"):
        # filter on single seqid
        df_query = "seqid == " + (str(record.id))
        seqid_df = df.query(df_query)

        # sort by start coordinate
        seqid_df.sort_values(by=["start"], ascending=True)

        # assuming start of CDS is in-frame
        # determine the frame of each subsequent exon

        # initialize counter to store cumulative nucleotide length
        nuc_length = 0

        # initialize offset to 0 for first annotation
        offset = 0

        # set counter for each CDS annotation for a record
        ct = 0

        for idx, row in seqid_df.iterrows():
            # read each GFF row into a variable
            # change these variables when start and stop codon positions change
            seqid = row["seqid"]
            source = row["source"]
            gff_type = row["type"]
            start = row["start"]
            end = row["end"]
            score = row["score"]
            strand = row["strand"]
            phase = row["phase"]
            attributes = row["attributes"]

            # trim first CDS annotation (exon 1 by sorting) so it begins with a methionine
            if ct == 0:
                # translate sequence, using offset from previous annotation
                # GFF is 1-based
                # biopython is 0-based
                translation = record.seq[
                    row["start"] - 1 - offset : row["end"]
                ].translate()
                methionine_start_pos = translation.find("M")
                start = start + (
                    methionine_start_pos * 3
                )  # set new start position in exon 1

            # iterate CDS annotation counter
            ct += 1

            # write new GFF file
            handle.write(
                str(seqid)
                + "\t"
                + str(source)
                + "\t"
                + str(gff_type)
                + "\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t"
                + str(score)
                + "\t"
                + str(strand)
                + "\t"
                + str(phase)
                + "\t"
                + str(attributes)
                + "\n",
            )

            # I tried to handle the stop codon trimming automatically, but could not get accuracy high enough
            # will trim 3' end manually
