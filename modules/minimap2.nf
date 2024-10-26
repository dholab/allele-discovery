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
    ${gdna_reference_fasta} \
    ${putative_alleles} \
    -ax splice -N 10000 \
    > "all_mappings.sam"
    """

}

process GENOTYPE_AMPLICONS {

    tag "${id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple path(genotyping_fasta), path(fasta_idx), path(trimmed_fastq)

    output:
    path "${id}.sam"

    script:
    id = file(trimmed_fastq).getSimpleName()
    """
    minimap2 \
    ${genotyping_fasta} \
    ${trimmed_fastq} \
    -ax map-hifi --eqx -t ${task.cpus} \
    | reformat.sh \
    in=stdin.sam \
    out=${id}.sam \
    ref=${genotyping_fasta} \
    noheader=t \
    threads=${task.cpus} \
    ow=t
    """

}

