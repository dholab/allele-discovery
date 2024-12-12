process FILTER_ALIGNMENTS {

    tag "${file_label}"

    publishDir params.filtered_geno, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path bam
    each path(refseq)

    output:
    path "${file_label}.filtered.bam"

    script:
    file_label = file(bam).getSimpleName()
    """
	samtools faidx ${refseq} && \
    filter_alignments.py \
    --input_bam ${bam} \
    --reference_fasta ${refseq} \
    --output_bam ${file_label}.filtered.bam
    """
}
