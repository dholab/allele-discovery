process VALIDATE_READS {

    /* */

    tag "${label}"

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    cpus 4

    input:
    tuple val(label), path(seq_file)

    output:
    tuple val(label), path("${label}.validated.fastq.gz"), val("success!")

    script:
    if ( file(seq_file).getName().contains(".fastq") || file(seq_file).getName().contains(".fq") )
        """
        seqkit seq --validate-seq ${seq_file} -o ${label}.validated.fastq.gz
        """
    else if ( file(seq_file).getName().endsWith(".bam") )
        """
        samtools quickcheck -v -u ${seq_file} && \
        samtools fastq --threads ${task.cpus} ${seq_file} | bgzip -o ${label}.validated.fastq.gz
        """
    else
        error "Unrecognized file format provided."

}

