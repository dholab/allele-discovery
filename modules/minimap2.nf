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

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    each path(genotyping_fasta)
    path trimmed_fastq

    output:
    path "${id}.sam"

    script:
    id = file(trimmed_fastq).getSimpleName()
    """
    minimap2 \
    ${genotyping_fasta} \
    ${trimmed_fastq} \
    -ax map-hifi --eqx -t 3 \
    | reformat.sh \
    in=stdin.sam \
    out=${id}.sam \
    ref=${genotyping_fasta} \
    noheader=t \
    subfilter=0 \
    threads=1 \
    ow=t \
    -da
    """

}

