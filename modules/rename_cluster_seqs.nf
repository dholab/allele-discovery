process RENAME_PUTATIVE_ALLELE_CLUSTERS {

    publishDir params.shared_clusters, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 1

    input:
    path putative_alleles

    output:
    path "putative_alleles.fasta"

    script:
    """
    #!/usr/bin/env python3

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    from datetime import datetime

    # because clustering of shared allele is nondeterminisitc, add timestamp to output FASTA
    timestamp=datetime.now()
    current_time=timestamp.strftime("%Y%m%d%H%M%S")

    # parse input FASTA
    with open("putative_alleles.fasta", "w") as handle:
        for idx, record in enumerate(SeqIO.parse("${putative_alleles}", "fasta")):
            record.id = str(current_time) + '-' + str(idx)
            SeqIO.write(record, handle, "fasta")
    """
}
