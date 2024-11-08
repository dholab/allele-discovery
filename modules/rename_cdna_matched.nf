process RENAME_CDNA_MATCHED_FASTA {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path no_gdna_match
    path matches

    output:
    path "cdna_matches.fasta"

    script:
    """
    #!/usr/bin/env python3

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import csv

    # read input FASTA line-by-line
    for record in SeqIO.parse("${no_gdna_match}", "fasta"):

        # parse file with gdna sequences that match cdna sequences
        with open("${matches}") as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:

                # test if name of sequence in cdna match file matches gdna sequence name
                if row[0] == record.name:
                    # update name of sequence in output file
                    record.description = row[1] + '|' + row[0]

                    # write to file
                    with open("cdna_matches.fasta", "w") as handle:
                        SeqIO.write(record, handle, "fasta")
    """

}
