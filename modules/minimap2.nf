process MAP_CLUSTERS_TO_FULL_LENGTH_GDNA {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 3

    input:
    each path(gdna_reference_fasta)
    path putative_alleles

    output:
    path "all_mappings.sam"

    script:
    """
    minimap2 \
    -t ${task.cpus} \
    -ax splice -N 10000 \
    ${gdna_reference_fasta} \
    ${putative_alleles} \
    > "all_mappings.sam"
    """

}

process GENOTYPE_AMPLICONS {

    tag "${id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(id), path(trimmed_fastq), path(genotyping_fasta), path(fasta_idx)

    output:
    path "${id}.sam"

    script:
    """
    minimap2 \
    ${genotyping_fasta} \
    ${trimmed_fastq} \
    -ax map-hifi --eqx -t ${task.cpus} \
    > ${id}.sam
    """

}

process GENOTYPE_AMPLICON_CLUSTERS {

    tag "${id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(id), path(trimmed_fastq), path(genotyping_fasta), path(fasta_idx)

    output:
    path "${id}.sam"

    script:
    """
    minimap2 \
    ${genotyping_fasta} \
    ${trimmed_fastq} \
    -ax map-hifi --eqx -t ${task.cpus} \
    > ${id}.sam
    """

}
